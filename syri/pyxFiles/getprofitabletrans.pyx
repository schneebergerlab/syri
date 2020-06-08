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

from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map as cpp_map
from libcpp.deque cimport deque as cpp_deq
from libcpp.queue cimport queue as cpp_que
from libcpp.vector cimport vector as cpp_vec
from libcpp cimport bool as cpp_bool

cimport numpy as np
cimport cython
cdef float[:]                                    dist, dist_bkp

cpdef makeBlocksTree_ctx(long[:] astart, long[:] aend, long[:] bstart, long[:] bend, long[:] bdir, np.ndarray achr, np.ndarray bchr, int threshold):
    """Compute whether two alignments can be part of one translation block. For this:
       they should be syntenic with respect to each other.
    """

    cdef:
        Py_ssize_t                                  i,j, n = len(astart)
        cpp_map[long, cpp_deq[long]]                out
        cpp_map[long, cpp_deq[long]].iterator       it
        long[:]                                     achr_int, bchr_int
    _, achr_int = np.unique(achr, return_inverse=True)
    _, bchr_int = np.unique(bchr, return_inverse=True)
    for i in range(n):
        for j in range(i,n):
            if achr_int[i] != achr_int[j]:
                continue
            if bchr_int[i] != bchr_int[j]:
                continue
            if (astart[j] - astart[i]) > threshold:
                if (aend[j] - aend[i]) > threshold:
                    if bdir[i] == 1:
                        if (bstart[j] - bstart[i]) > threshold:
                            if (bend[j] - bend[i]) > threshold:
                                out[i].push_back(j)
                    else:
                        if (bstart[i] - bstart[j]) > threshold:
                            if (bend[i] - bend[j]) > threshold:
                                out[i].push_back(j)

    it = out.begin()
    outdict = {}
    while it != out.end():
        outdict[deref(it).first] = [deref(it).second[i] for i in range(<Py_ssize_t> deref(it).second.size())]
        inc(it)
    return outdict


cpdef getProfitableTrans(cpp_map[long, cpp_set[long]] graph, long[:] astart, long[:] aend, long[:] bstart, long[:] bend, np.ndarray achr, np.ndarray bchr, float[:] iden, long[:] alen, long[:] blen, long[:] inastart, long[:] inaend, long[:] inbstart, long[:] inbend, np.ndarray inachr, np.ndarray inbchr, long tUC, float tUP, int isinv = 0, brk = -1):
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
        float                                       ascore, bscore, agap, bgap
        long                                        al, bl, au, bu
        float                                       score
        float[:]                                    weight
        float[:]                                    dist, dist_bkp
        long[:]                                     pred, pred_bkp
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

    print('Start', str(datetime.now()))


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

    print('Sorted Keys', str(datetime.now()))

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

    print('Removed Children', str(datetime.now()))
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

    print('Got topological order', str(datetime.now()))
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

    print('Starting calculation', str(datetime.now()))

    print(len(n), edgecnt, str(datetime.now()))

    out = deque()
    count = 0
    dist_bkp = np.array([-np.float32('inf')]*<Py_ssize_t>len(n), dtype = np.float32)
    pred_bkp = np.array([-1]*<Py_ssize_t>len(n), dtype = np.int)
    for i in n:
        print(i, str(datetime.now()))
        st = datetime.now()
       
        count+=1
        if count % 1000 == 0:
            print(count, str(datetime.now()))
        if count == brk:
            break

        # if i%100 == 0:
        #     print(i, str(datetime.now()))
        print(i, 'I1', str(datetime.now() - st))
        nodepath.clear()
        pred = pred_bkp.copy()
        dist = dist_bkp.copy()
        dist[i] = 0

        # Process vertices in topological order
        
        print(i, 'I2', str(datetime.now() - st))
        """
        for j in n:
            for k in range(id, edgecnt):
                cnt+=1
                if source[k] != topo[j]:
                    break
                if dist[target[k]] > dist[source[k]] + weight[k]:   
                    dist[target[k]] = dist[source[k]] + weight[k]
                    pred[target[k]] = source[k]
                id+=1
       
        """
        #cnt = 0
        st = datetime.now()
        a = deque()
        a.append(i)
        for k in range(edgecnt):
            if dist[target[k]] < dist[source[k]] + weight[k]:
                a.append(target[k])
                dist[target[k]] = dist[source[k]] + weight[k]
                pred[target[k]] = source[k]
        #print(i, cnt, edgecnt)
        print(i, 'A', str(datetime.now() - st))

        st = datetime.now()
        cnt = 0
        #for j in n:
        for j in set(a):
            # Find all connected paths which are profitable
            #if dist[topo[j]] != float("inf"):
            if True:
                cnt += 1
                #current = topo[j]
                current = j
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
                nodepath[j] = path
                #nodepath[topo[j]] = path
                path.push_front(i)                          ## Found the best path between two co-linear nodes

                #print(i, 'A', str(datetime.now() - st))
                ## Calculate score of the identified path
                ascore = float(alen[path[0]])
                bscore = float(blen[path[0]])
                agap = 0
                bgap = 0

                if path.size() > 1:
                    for k in range(1, <Py_ssize_t> path.size()):
                        ascore += float(alen[path[k]])
                        bscore += float(blen[path[k]])
                        agap += 0 if 0 > (astart[path[k]] - aend[path[k-1]]) else float(astart[path[k]] - aend[path[k-1]])
                        if not isinv:
                            bgap += 0 if 0 > (bstart[path[k]] - bend[path[k-1]]) else float(bstart[path[k]] - bend[path[k-1]])
                        else:
                            bgap += 0 if 0 > (bstart[path[k-1]] - bend[path[k]]) else float(bstart[path[k-1]] - bend[path[k]])

                score = min(((ascore - agap)/ascore),((bscore - bgap)/bscore))

                #print(i, 'B', str(datetime.now() - st))
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
                    #print(i, 'C', str(datetime.now() - st))

                #Trans block is selected IFF either the unique region on any genome is larger than tUC
                # or length of unique region on a genome is larger than tUP times the length of
                # the overlapping region on that genome
                    if au > tUC or bu > tUC or au > tUP*al or bu > tUP*bl:
                        out.append([path[k] for k in range(<Py_ssize_t> path.size())])
                        #print(i, 'D', str(datetime.now() - st))
        print(i, 'D', str(datetime.now() - st))
    print('End', datetime.now())
    return out
	


import os
os.chdir('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/greg_bug/')
import pickle
# with open('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/greg_bug/pickle.data', 'wb') as fout:
#     pickle.dump([orderedBlocks, outOrderedBlocks, invertedBlocks, outInvertedBlocks, annoCoords, threshold, tUC, tUP], fout)

dat = pickle.load(open('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/greg_bug/pickle.data', 'rb'))

%load_ext cython

orderedBlocks = dat[0].copy()
outOrderedBlocks = dat[1].copy()
annoCoords = dat[4]
tUC = dat[6]
tUP = dat[7]

outOrderedBlocks_filtered = {}

outOrderedBlocks = dat[1]
for i in sorted(list(outOrderedBlocks.keys())):
    #print(i)
    j = 0
    while j < len(outOrderedBlocks[i]):
        try:
            outOrderedBlocks[i] = sorted(set(outOrderedBlocks[i]) - set(outOrderedBlocks[outOrderedBlocks[i][j]]))
        except KeyError:
            pass
        j+=1
    g = []
    for j in outOrderedBlocks[i]:
        if orderedBlocks.iloc[j]['aStart'] - orderedBlocks.iloc[i]['aEnd'] < 1000000:
            if orderedBlocks.iloc[j]['bStart'] - orderedBlocks.iloc[i]['bEnd'] < 1000000:
                g.append(j)
    outOrderedBlocks[i]  = g

isinv = 0
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
                                     annoCoords.sort_values(['aChr', 'aStart','aEnd']).aStart.values,
                                     annoCoords.sort_values(['aChr', 'aStart','aEnd']).aEnd.values,
                                     annoCoords.sort_values(['bChr', 'bStart','bEnd']).bStart.values,
                                     annoCoords.sort_values(['bChr', 'bStart','bEnd']).bEnd.values,
                                     annoCoords.sort_values(['aChr', 'aStart','aEnd']).aChr.values,
                                     annoCoords.sort_values(['bChr', 'bStart','bEnd']).bChr.values,
                                     tUC,
                                     tUP,
                                     isinv,
                                     brk = 1000)