%%cython
# distutils: language = c++
import pandas as pd
import numpy as np
from collections import deque
from datetime import datetime
import logging
import time
from syri.bin.func.myUsefulFunctions import *
from multiprocessing import Pool
from functools import partial
from gc import collect

from syri.pyxFiles.function cimport getOverlapWithSynBlocks, getmeblocks
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map as cpp_map
from libcpp.deque cimport deque as cpp_deq
from libcpp.queue cimport queue as cpp_que
from libcpp.vector cimport vector as cpp_vec
from libcpp cimport bool as cpp_bool

cimport numpy as np
cimport cython

cpdef makeBlocksTree(long[:] aStart, long[:] aEnd, long[:] bStart, long[:] bEnd, int threshold, long[:] left, long[:] right):
    """Compute whether two alignments can be part of one translation block. For this:
        the alignments should not be separated by any inPlaceBlock on both ends and
        they should be collinear with respect to each other.

       Returns
       --------
       outOrderedBlocks: pandas DataFrame,
           Dataframe of type Object. Lower half is NA, upper half contains whether two
           alignments can be connected (True) or not (False).
    """

    cdef:
        ssize_t                                     i,j, n = len(aStart)
        cpp_map[long, cpp_deq[long]]                out
        cpp_map[long, cpp_deq[long]].iterator       it
    # cdef np.ndarray allRanges = np.array([range(left[i]+1, right[i]) for i in range(n)])

    for i in range(<Py_ssize_t> n):
        for j in range(<Py_ssize_t> i, <Py_ssize_t> n):
            if (aStart[j] - aStart[i]) > threshold:
                if (aEnd[j] - aEnd[i]) > threshold:
                    if (bStart[j] - bStart[i]) > threshold:
                        if (bEnd[j] - bEnd[i]) > threshold:
                            # Alignments could form a block when they are not separated by inPlaceBlocks. For this, we check whether the alignments have a common intersecting inPlaceAlignments.
                            if right[i]-1 >= left[j]+1 and right[i] <= right[j]:
                                out[i].push_back(j)
                            elif left[i] >= left[j] and left[i]+1 <= right[j]-1:
                                out[i].push_back(j)

    it = out.begin()
    outdict = {}
    while it != out.end():
        outdict[deref(it).first] = [deref(it).second[i] for i in range(<Py_ssize_t> deref(it).second.size())]
        inc(it)
    return(outdict)




def getTransSynOrientation(inPlaceData, transData, threshold, ctx = False):
    """ To get the nearest left and right inPlaceBlocks for all the translated blocks
        in the transData object.
    """
    inPlaceBlocks = inPlaceData.copy()
    if not ctx:
        transRowCount = transData.shape[0]
        transPositions = dict()

        for i in range(transRowCount):
            row = transData.iloc[i]

            upSyn = intersect(np.where(inPlaceBlocks.aStart < (row.aStart - threshold))[0],
                          np.where(inPlaceBlocks.aEnd < (row.aEnd - threshold))[0],
                          np.where(inPlaceBlocks.bStart < (row.bStart - threshold))[0],
                          np.where(inPlaceBlocks.bEnd < (row.bEnd - threshold))[0])

            downSyn = intersect(np.where(inPlaceBlocks.aStart > (row.aStart + threshold))[0],
                                      np.where(inPlaceBlocks.aEnd > (row.aEnd + threshold))[0],
                                      np.where(inPlaceBlocks.bStart > (row.bStart + threshold))[0],
                                      np.where(inPlaceBlocks.bEnd > (row.bEnd + threshold))[0])

            upBlock = max(upSyn) if len(upSyn) > 0 else -1
            downBlock = min(downSyn) if len(downSyn) > 0 else len(inPlaceBlocks)
            transPositions[i] = [upBlock, downBlock]
        transPositions = pd.DataFrame(transPositions).transpose()
        return transPositions




cpdef getProfitableTrans(cpp_map[long, cpp_set[long]] graph, long[:] astart, long[:] aend, long[:] bstart, long[:] bend, np.ndarray achr, np.ndarray bchr, float[:] iden, long[:] alen, long[:] blen, long[:] inastart, long[:] inaend, long[:] inbstart, long[:] inbend, np.ndarray inachr, np.ndarray inbchr, long tUC, float tUP, int isinv = 0):
    """
    Input:
     1) dictionary in which each key corresponds to trans alignment and the corresponding values are the alignments which are colinear with key.
     2) Coordinates of the trans alignments
     3) Coordinates of inplace blocks. Sorted separately for reference and query genome.
     4) tUC, tUP

    Output:
     1) All trans blocks (groups of alignments) with high positive score.
    """
    cdef:
        Py_ssize_t                                  i, j, k
        long                                        id, current
        unsigned long                               cnt
        long                                        nodecnt, edgecnt, len_in = len(inastart)
        long                                        ast, aen, bst, ben
        long                                        ascore, bscore, agap, bgap
        long                                        al, bl, au, bu
        float                                       score
        float[:]                                    weight
        float[:]                                    dist
        long[:]                                     pred
        long[:]                                     topo, indegree, source, target
        long[:]                                     n = np.array(range(len(astart)), dtype='int')
        long[:]                                     achrint, bchrint, inachrint, inbchrint
        cpp_set[long]                               rem
        cpp_que[long]                               q, toporder
        cpp_deq[long]                               path, r_path
        cpp_map[long, cpp_deq[long]]                nodepath
        cpp_map[long, cpp_vec[long]]                almntdata
        cpp_set[long].iterator                      set_it, set_it2
        cpp_deq[long].reverse_iterator              deq_rit
        cpp_map[long, cpp_set[long]].iterator       mapit


    nodecnt = len(n)
    topo = indegree = np.zeros(nodecnt, dtype='int64')

    ## Check that the keys are sorted in the graph
    mapit = graph.begin()
    id = -1
    while mapit != graph.end():
        if id >= deref(mapit).first:
            print('ERROR: unsorted outOrderedBlocks')
        else:
            id = deref(mapit).first
        inc(mapit)


    ## Remove children for which another child is present between parent and itself
    mapit = graph.begin()
    while mapit != graph.end():
        rem.clear()                     ## List of children to be removed
        set_it = deref(mapit).second.begin()
        while set_it != deref(mapit).second.end():
            if rem.count(deref(set_it)) == 0:
                if graph.count(deref(set_it)) == 1:
                    set_it2 = graph[deref(set_it)].begin()
                    while set_it2 !=  graph[deref(set_it)].end():
                        rem.insert(deref(set_it2))
                        inc(set_it2)
            inc(set_it)
        set_it = rem.begin()
        while set_it != rem.end():
            deref(mapit).second.erase(deref(set_it))
            inc(set_it)
        inc(mapit)


    # Get number of edges in the graph
    edgecnt = 0
    mapit = graph.begin()
    while mapit != graph.end():
        edgecnt += <Py_ssize_t> deref(mapit).second.size()
        inc(mapit)

    # Find the topological order in which the nodes should be processed to find paths

    ## Get indegree for all nodes (number of nodes == number of aligments)
    mapit = graph.begin()
    while mapit != graph.end():
        set_it = deref(mapit).second.begin()
        while set_it != deref(mapit).second.end():
            indegree[deref(set_it)]+=1
            inc(set_it)
        inc(mapit)


    ## Push all nodes with indegree=0 in queue
    for i in n:
        if indegree[i] == 0:
            q.push(i)
    cnt = 0

    ## Get topological ordering for the nodes
    while q.size() > 0:
        id = q.front()
        q.pop()
        toporder.push(id)
        if graph.count(id) == 1:
            set_it = graph[id].begin()
            while set_it != graph[id].end():
                indegree[deref(set_it)]-=1
                if indegree[deref(set_it)]==0:
                    q.push(deref(set_it))
                inc(set_it)
        cnt += 1
    if cnt != len(indegree):
        print('ERROR: Cycle found')
    if toporder.size() != len(topo):
        print('ERROR: topological ordering didnt cover all nodes')
    for i in n:
        topo[i] = toporder.front()
        toporder.pop()


    # Get order in which the edges need to be transversed
    source = np.zeros(edgecnt, dtype='int')
    target = np.zeros(edgecnt, dtype='int')
    weight = np.zeros(edgecnt, dtype='float32')
    id = 0
    for i in topo:
        if graph.count(i) == 1:
            set_it = graph[i].begin()
            while set_it != graph[i].end():
                source[id] = i
                target[id] = deref(set_it)
                weight[id] = (aend[i] + bend[i] - astart[i] - bstart[i] + 2) * iden[i]
                id+=1
                inc(set_it)


    ## Convert arary of 'string' chromosome ids to much faster numeric chromosome ids
    if list(np.unique(achr)) == list(np.unique(bchr)) == list(np.unique(inachr)) == list(np.unique(inbchr)):
        _, achrint      = np.unique(achr, return_inverse = True)
        _, bchrint      = np.unique(bchr, return_inverse = True)
        _, inachrint    = np.unique(inachr, return_inverse = True)
        _, inbchrint    = np.unique(inbchr, return_inverse = True)
    else:
        unchrs = np.unique(list(np.unique(achr)) + list(np.unique(bchr)) + list(np.unique(inachr)) + list(np.unique(inbchr)))
        unchrdict = {unchrs[i]:i for i in range(len(unchrs))}
        achrint     = np.array([unchrdict[c] for c in achr], np.int)
        bchrint     = np.array([unchrdict[c] for c in bchr], np.int)
        inachrint   = np.array([unchrdict[c] for c in inachr], np.int)
        inbchrint   = np.array([unchrdict[c] for c in inbchr], np.int)


    ## For each alignment/node calculate the number of bases which are not overlapping with the in-place blocks
    for i in n:
        ascore = 0
        bscore = 0
        ast = astart[i]
        aen = aend[i]
        bst = bstart[i]
        ben = bend[i]

        # print(ast, aen, bst, ben)
        for j in range(len_in):
            if achrint[i] == inachrint[j]:
                if inastart[j] > aen:
                    break
                if inaend[j] < ast:
                    continue
                if inastart[j] < ast:
                    if inaend[j] < aen:
                        ast = inaend[j]+1
                    else:
                        ast = aen+1
                        break
                else:
                    ascore += inastart[j] - ast
                    if inaend[j] < aen:
                        ast = inaend[j]+1
                    else:
                        ast = aen+1
                        break
        ascore += aen - ast + 1

        for j in range(len_in):
            if bchrint[i] == inbchrint[j]:
                if inbstart[j] > ben:
                    break
                if inbend[j] < bst:
                    continue
                if inbstart[j] < bst:
                    if inbend[j] < ben:
                        bst = inbend[j]+1
                    else:
                        bst = ben+1
                        break
                else:
                    bscore += inbstart[j] - bst
                    if inbend[j] < ben:
                        bst = inbend[j]+1
                    else:
                        bst = ben+1
                        break
        bscore += ben - bst + 1
        almntdata[i] = np.array([aend[i]-astart[i]+1, bend[i]-bstart[i]+1, ascore, bscore], dtype=np.int)

    out = deque()
    for i in n:
        # if i%100 == 0:
        #     print(i, str(datetime.now()))
        nodepath.clear()
        pred = np.array([-1]*<Py_ssize_t>len(n), dtype = np.int)
        dist = np.array([np.float32('inf')]*<Py_ssize_t>len(n), dtype = np.float32)
        dist[i] = 0

        # Process vertices in topological order
        id = 0
        for j in n:
            for k in range(id, edgecnt):
                if source[k] != topo[j]:
                    break
                if dist[target[k]] > dist[source[k]] + weight[k]:
                    dist[target[k]] = dist[source[k]] + weight[k]
                    pred[target[k]] = source[k]
                id+=1
        for j in n:
            # Find all connected paths which are profitable
            if dist[topo[j]] != float("inf"):
                current = topo[j]
                path.clear()
                while current!=i:
                    if nodepath.count(current) > 0:
                        deq_it = nodepath[current].rbegin()
                        while deq_it != nodepath[current].rend():
                            path.push_front(deref(deq_it))
                            inc(deq_it)
                        break
                    path.push_front(current)
                    current = pred[current]
                nodepath[topo[j]] = path
                path.push_front(i)                          ## Found the best path between two co-linear nodes

                ## Calculate score of the identified path
                ascore = alen[path[0]]
                bscore = blen[path[0]]
                agap = 0
                bgap = 0

                if path.size() > 1:
                    for k in range(1, <Py_ssize_t> path.size()):
                        ascore += alen[path[k]]
                        bscore += blen[path[k]]
                        agap += 0 if 0 > (astart[path[k]] - aend[path[k-1]]) else astart[path[k]] - aend[path[k-1]]
                        if not isinv:
                            bgap += 0 if 0 > (bstart[path[k]] - bend[path[k-1]]) else bstart[path[k]] - bend[path[k-1]]
                        else:
                            bgap += 0 if 0 > (bstart[path[k-1]] - bend[path[k]]) else bstart[path[k-1]] - bend[path[k]]

                score = min(((ascore - agap)/ascore),((bscore - bgap)/bscore))


                # print([path[id] for id in range(path.size())], score)

                ## Check if the alignments of the path explain sufficient unique alignments and are not excessively overlapped

                # Only process paths for which the alignments explain more than gaps
                if score > 0:
                    al = 0
                    bl = 0
                    au = 0
                    bu = 0
                    for k in range(<Py_ssize_t> path.size()):
                        al += almntdata[path[k]][0]
                        bl += almntdata[path[k]][1]
                        au += almntdata[path[k]][2]
                        bu += almntdata[path[k]][3]

                #Trans block is selected IFF either the unique region on any genome is larger than tUC
                # or length of unique region on a genome is larger than tUP times the length of
                # the overlapping region on that genome
                    if au > tUC or bu > tUC or au > tUP*al or bu > tUP*bl:
                        out.append([path[k] for k in range(<Py_ssize_t> path.size())])
    return out


def blocksdata(outPlaceBlocks, inPlaceBlocks, threshold, tUC, tUP, chromo):
    logger = logging.getLogger('blocksdata.'+ chromo)
    orderedBlocks = outPlaceBlocks[outPlaceBlocks.bDir == 1]
    invertedBlocks = outPlaceBlocks[outPlaceBlocks.bDir == -1]

    logger.debug('Number of directed alignments: '+ str(orderedBlocks.shape[0]) + '. Number of inverted alignments: '+ str(orderedBlocks.shape[0]))

    if len(orderedBlocks) > 0:
        transBlocksNeighbours = getTransSynOrientation(inPlaceBlocks, orderedBlocks, threshold)
        outOrderedBlocks = makeBlocksTree(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, threshold, transBlocksNeighbours[0].values, transBlocksNeighbours[1].values)
        ## need to have ref coords and query coords separately sorted
        transBlocks = getProfitableTrans(outOrderedBlocks,
                                         orderedBlocks.aStart.values,
                                         orderedBlocks.aEnd.values,
                                         orderedBlocks.bStart.values,
                                         orderedBlocks.bEnd.values,
                                         orderedBlocks.aChr.values,
                                         orderedBlocks.bChr.values,
                                         orderedBlocks.iden.values.astype('float32'),
                                         orderedBlocks.aLen.values,
                                         orderedBlocks.bLen.values,
                                         inPlaceBlocks.sort_values(['aChr', 'aStart','aEnd']).aStart.values,
                                         inPlaceBlocks.sort_values(['aChr', 'aStart','aEnd']).aEnd.values,
                                         inPlaceBlocks.sort_values(['bChr', 'bStart','bEnd']).bStart.values,
                                         inPlaceBlocks.sort_values(['bChr', 'bStart','bEnd']).bEnd.values,
                                         inPlaceBlocks.sort_values(['aChr', 'aStart','aEnd']).aChr.values,
                                         inPlaceBlocks.sort_values(['bChr', 'bStart','bEnd']).bChr.values,
                                         tUC,
                                         tUP)
    else:
        transBlocks = []

    logger.debug(str(len(transBlocks)) + ' candidate directed TDs found')
    if len(invertedBlocks) > 0:
        invertedCoords = invertedBlocks.copy()
        invertedCoords.bStart = invertedCoords.bStart + invertedCoords.bEnd
        invertedCoords.bEnd = invertedCoords.bStart - invertedCoords.bEnd
        invertedCoords.bStart = invertedCoords.bStart - invertedCoords.bEnd
        invTransBlocksNeighbours = getTransSynOrientation(inPlaceBlocks, invertedCoords, threshold)

        invertedCoords = invertedBlocks.copy()
        maxCoords = np.max(np.max(invertedCoords[["bStart","bEnd"]]))
        invertedCoords.bStart = maxCoords + 1 - invertedCoords.bStart
        invertedCoords.bEnd = maxCoords + 1 - invertedCoords.bEnd
        outInvertedBlocks = makeBlocksTree(invertedCoords.aStart.values, invertedCoords.aEnd.values, invertedCoords.bStart.values, invertedCoords.bEnd.values, threshold, invTransBlocksNeighbours[0].values, invTransBlocksNeighbours[1].values)

        invertedCoords = invertedBlocks.copy()
        invertedCoords.bStart = invertedCoords.bStart + invertedCoords.bEnd
        invertedCoords.bEnd = invertedCoords.bStart - invertedCoords.bEnd
        invertedCoords.bStart = invertedCoords.bStart - invertedCoords.bEnd
        invTransBlocks = getProfitableTrans(outInvertedBlocks,
                                            invertedCoords.aStart.values,
                                            invertedCoords.aEnd.values,
                                            invertedCoords.bStart.values,
                                            invertedCoords.bEnd.values,
                                            invertedCoords.aChr.values,
                                            invertedCoords.bChr.values,
                                            invertedCoords.iden.values.astype('float32'),
                                            invertedCoords.aLen.values,
                                            invertedCoords.bLen.values,
                                            inPlaceBlocks.sort_values(['aChr', 'aStart','aEnd']).aStart.values,
                                            inPlaceBlocks.sort_values(['aChr', 'aStart','aEnd']).aEnd.values,
                                            inPlaceBlocks.sort_values(['bChr', 'bStart','bEnd']).bStart.values,
                                            inPlaceBlocks.sort_values(['bChr', 'bStart','bEnd']).bEnd.values,
                                            inPlaceBlocks.sort_values(['aChr', 'aStart','aEnd']).aChr.values,
                                            inPlaceBlocks.sort_values(['bChr', 'bStart','bEnd']).bChr.values,
                                            tUC,
                                            tUP,
                                            isinv = 1)
    else:
        invTransBlocks = []
    logger.debug(str(len(transBlocks)) + ' candidate inverted TDs found')

    allTransBlocks, allTransIndexOrder = mergeTransBlocks(transBlocks, orderedBlocks, invTransBlocks, invertedBlocks)
    logger.debug('Finished')
    return (transBlocks, invTransBlocks, allTransBlocks, allTransIndexOrder)

