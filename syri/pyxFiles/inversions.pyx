# distutils: language = c++
import numpy as np
from igraph import Graph
from collections import deque
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vec
from libcpp.deque cimport deque as cpp_deq
from libcpp cimport bool as bool_t
from syri.bin.func.myUsefulFunctions import *
# from scipy.stats import *
import pandas as pd
from gc import collect
import logging
from datetime import datetime
from syri.pyxFiles.synsearchFunctions import apply_TS, alignmentBlock
from syri.pyxFiles.function cimport getConnectivityGraph

cimport numpy as np

np.random.seed(1)

cpdef invPath(cpp_map[long, cpp_vec[long]] invpos, long[:, :] neighbour, float[:] profit, long[:] aStart, long[:] aEnd, long[:] bStart, long[:] bEnd, long threshold):
    cdef:
        cpp_deq[long]               st, end, stb, endb, path
        long[:]                     parents
        float[:]                    totscore
        Py_ssize_t                  i, j, maxid
        int                         lp = invpos.size()
        float                       max

    for i in range(lp):
        st.push_back(aStart[invpos[i].front()])
        end.push_back(aEnd[invpos[i].back()])
        stb.push_back(bEnd[invpos[i].back()])
        endb.push_back(bStart[invpos[i].front()])

    totscore = profit
    parents = np.array([-1]*lp, dtype = 'int')

    for i in range(lp):
        for j in range(i, lp):
            if st[j] > end[i]-threshold:
                if stb[j] > endb[i] -threshold:
                    if profit[j] + totscore[i] > totscore[j]:
                        totscore[j] = profit[j] + totscore[i]
                        parents[j] = i

    maxid = -1
    max = -1

    for i in range(lp):
        if totscore[i] > i:
            max = totscore[i]
            maxid = i

    path.push_front(maxid)
    while parents[maxid] != -1:
        path.push_front(parents[maxid])
        maxid = parents[maxid]
    return [path[i] for i in range(path.size())]

cpdef getProfitable(invblocks, long[:] aStart, long[:] aEnd, long[:] bStart, long[:] bEnd, float[:] iDen, cpp_map[int, cpp_vec[long]] neighbourSyn, float[:] synBlockScore, long[:] aStartSyn, long[:] aEndSyn, long tUC, float tUP,  brk = -1):
    cdef:
        long                            i, j, k, l, current
        long                            n_topo, n_edges, n_syn
        long                            leftSyn, rightSyn, leftEnd, rightEnd, overlapLength
        float                           w, revenue, cost, profit
        Py_ssize_t                      index
        long[:]                         n = np.array(range(len(invblocks)), dtype=np.int)
        long[:]                         topo
        long[:]                         source, target
        long[:]                         pred
        float[:]                        dist
        float[:]                        weight
        long[:,:]                       nsynmap = np.zeros((neighbourSyn.size(), 2), dtype='int64')                  # neighbours of inverted alignments
        cpp_map[long, cpp_deq[long]]    nodepath
        cpp_deq[long]                   path, r_path
        cpp_deq[long]                   startA, endA, startB, endB
        cpp_deq[float]                  iden
        bool_t                          isMore

    n_syn = len(synBlockScore)
    invG = getConnectivityGraph(invblocks)
    out = deque()

    # get neighbours of inverted alignments
    for i in range(neighbourSyn.size()):
        nsynmap[i, 0] = neighbourSyn[i][0]
        nsynmap[i, 1] = neighbourSyn[i][1]

    ## Get topological ordering of the graph
    indegree = invG.vs.degree('IN')
    q = deque()
    toporder = deque()
    for i in n:
        if indegree[i] == 0:
            q.append(i)

    cnt = 0
    while len(q) > 0:
        u = q.popleft()
        toporder.append(u)
        for i in invG.neighbors(u,'OUT'):
            indegree[i] -= 1
            if indegree[i] == 0:
                q.append(i)
        cnt += 1

    if cnt != len(indegree):
        print('Cycle found')

    topo = np.array(toporder, dtype = np.int)
    n_topo = len(topo)

    # Get order in which the edges need to be transversed
    source = np.zeros_like(invG.es['source'], dtype=int)
    target = np.zeros_like(invG.es['source'], dtype=int)
    weight = np.zeros_like(invG.es['source'], dtype=np.float32)

    index = 0
    garb = invG.get_adjlist('OUT')
    for i in range(n_topo):
        if not len(garb[topo[i]]) > 0:
            continue
        w = invG.es[invG.get_eid(topo[i], garb[topo[i]][0])]['weight']
        for j in range(len(garb[topo[i]])):
            source[index] = topo[i]
            target[index] = garb[topo[i]][j]
            weight[index] = w
            index += 1
    n_edges = len(source)

    # Find shortest path to all other nodes from each node

    for i in n:
        if i%500 == 0:
            print(i, str(datetime.now()))
        nodepath.clear()
        pred = np.array([-1]*len(n), dtype = np.int)
        dist = np.array([np.float32('inf')]*len(n), dtype = np.float32)
        dist[i] = 0

        # Process vertices in topological order
        index = 0
        for j in range(n_topo):
            for k in range(index, n_edges):
                if source[k] != topo[j]:
                    break
                if dist[target[k]] > dist[source[k]] + weight[k]:
                    dist[target[k]] = dist[source[k]] + weight[k]
                    pred[target[k]] = source[k]
                index+=1

        for j in range(n_topo):
            # Find all connected paths which are profitable
            if dist[topo[j]] != float("inf"):
                current = topo[j]
                path.clear()
                while current!=i:
                    if nodepath.count(current) > 0:
                        for index in range(nodepath[current].size()):
                            path.push_back(nodepath[current][index])
                        # path.extend(nodepath[current].copy())
                        break
                    path.push_back(current)
                    current = pred[current]
                nodepath[topo[j]] = path
                path.push_back(i)
                r_path.clear()

                current = path.size()
                for index in range(path.size()):
                    r_path.push_back(path[current-index-1])

                # calculate revenue of the identified path
                if r_path.size() == 1:
                    revenue = iDen[r_path[0]]*(aEnd[r_path[0]] - aStart[r_path[0]] + 1 + bStart[r_path[0]] - bEnd[r_path[0]] + 1)
                else:
                    revenue = 0
                    # Initiate by adding coordinates of first alignment
                    startA.push_back(aStart[r_path[0]])
                    endA.push_back(aEnd[r_path[0]])
                    startB.push_back(bEnd[r_path[0]])
                    endB.push_back(bStart[r_path[0]])
                    iden.push_back(iDen[r_path[0]])

                    # Add remaining alignments' coordinates iteratively
                    for k in range(1, current):
                        l = r_path[k]
                        isMore = True if iDen[k] > iden.back() else False
                        if aStart[l] < endA.back():
                            # In case of overlapping bases, choose score of the alignment with higher identity
                            if isMore:
                                endA.pop_back()
                                endA.push_back(aStart[l])
                                startA.push_back(aStart[l])
                                endA.push_back(aEnd[l])
                            else:
                                startA.push_back(endA.back())
                                endA.push_back(aEnd[l])
                        else:
                            startA.push_back(aStart[l])
                            endA.push_back(aEnd[l])

                        if bStart[l] > startB.back():
                            # In case of overlapping bases, choose score of the alignment with higher identity
                            if isMore:
                                startB.pop_back()
                                startB.push_back(bStart[l])
                                startB.push_back(bEnd[l])
                                endB.push_back(bStart[l])
                            else:
                                endB.push_back(startB.back())
                                startB.push_back(bEnd[l])
                        else:
                            startB.push_back(bEnd[l])
                            endB.push_back(bStart[l])
                        iden.push_back(iDen[l])

                    if startA.size() == endA.size() == startB.size() == endB.size() == iden.size():
                        for k in range(iden.size()):
                            revenue += iden[k]*((endA[k] - startA[k] + 1) + (endB[k] - startB[k] + 1))
                        startA.clear()
                        endA.clear()
                        startB.clear()
                        endB.clear()
                        iden.clear()
                    else:
                        print('ERROR in calculating revenue')

                # Calculate cost of the identified path

                # Get left and right syntenic neighbours
                leftSyn = nsynmap[r_path.front()][0] if nsynmap[r_path.front()][0] < nsynmap[r_path.back()][0] else nsynmap[r_path.back()][0]
                rightSyn = nsynmap[r_path.front()][1] if nsynmap[r_path.front()][1] > nsynmap[r_path.back()][1] else nsynmap[r_path.back()][1]

                #cost of removing all intersecting neighbours
                cost = 0
                for k in range(leftSyn+1, rightSyn):
                    cost += synBlockScore[k]

                # Check whether inversion is overlapping. If yes, then check for uniqueness. If not sufficiently unique, then set high cost.
                leftEnd = aEndSyn[leftSyn] if leftSyn > -1 else 0
                rightEnd = aStartSyn[rightSyn] if rightSyn < n_syn else aEnd[r_path.back()]
                if rightEnd - leftEnd <= tUC:
                    overlapLength = (leftEnd - aStart[r_path.front()]) + (aEnd[r_path.back()] - rightEnd)
                    if overlapLength/(rightEnd - leftEnd) < tUP:
                        cost = 10000000000000

                # Select those candidate inversions for which the score of
                #  adding them would be at least 10% better than the score
                #  of syntenic regions needed to be removed
                if revenue > 1.1*cost:
                    out.append(([r_path[k] for k in range(current)], revenue - cost, leftSyn, rightSyn))
        if i == brk:
            return out
    return out


cpdef getInvBlocks(invTree, invertedCoordsOri):
    cdef int nrow, i, child
    nrow = invTree.shape[0]
    invBlocks = [alignmentBlock(i, np.where(invTree.iloc[i,] == True)[0], invertedCoordsOri.iloc[i]) for i in range(nrow)]

    for block in invBlocks:
        i = 0
        while(i < len(block.children)):
            block.children = list(set(block.children) - set(invBlocks[block.children[i]].children))
            i+=1
        block.children.sort()

        for child in block.children:
            invBlocks[child].addParent(block.id)
    return(invBlocks)


cpdef dict getNeighbourSyn(np.ndarray aStartInv, np.ndarray aEndInv, np.ndarray bStartInv, np.ndarray bEndInv, np.ndarray indexInv, np.ndarray bDirInv, np.ndarray aStartSyn, np.ndarray aEndSyn, np.ndarray bStartSyn, np.ndarray bEndSyn, np.ndarray indexSyn, np.ndarray bDirSyn, int threshold):
    cdef:
        cdef Py_ssize_t i, j, index
        dict neighbourSyn = dict()
        int upBlock, downBlock
        list upSyn, downSyn
    for i in range(len(indexInv)):
        index = indexInv[i]
        upSyn = np.where(indexSyn < index)[0].tolist()
        downSyn = np.where(indexSyn > index)[0].tolist()

        upBlock  = -1
        downBlock = len(indexSyn)
        for j in upSyn[::-1]:
            if bDirSyn[j] == bDirInv[i]:
                if (aStartInv[i] - aStartSyn[j]) > threshold and (aEndInv[i] - aEndSyn[j]) > threshold and (bStartInv[i] - bStartSyn[j]) > threshold and (bEndInv[i] - bEndSyn[j]) > threshold:
                    upBlock = j
                    break
            else:
                if (aStartInv[i] - aStartSyn[j]) > threshold and (aEndInv[i] - aEndSyn[j]) > threshold and (bEndInv[i] - bStartSyn[j]) > threshold and (bStartInv[i] - bEndSyn[j]) > threshold:
                    upBlock = j
                    break

        for j in downSyn:
            if bDirSyn[j] == bDirInv[i]:
                if (aStartSyn[j] - aStartInv[i]) > threshold and (aEndSyn[j] - aEndInv[i]) > threshold and (bStartSyn[j] - bStartInv[i]) > threshold and (bEndSyn[j] - bEndInv[i]) > threshold:
                    downBlock = j
                    break
            else:
                if (aStartSyn[j] - aStartInv[i]) > threshold and (aEndSyn[j] - aEndInv[i]) > threshold and (bStartSyn[j] - bEndInv[i]) > threshold and (bEndSyn[j] - bStartInv[i]) > threshold:
                    downBlock = j
                    break
        neighbourSyn[i] = [upBlock, downBlock]
    return(neighbourSyn)



def getInversions(coords,chromo, threshold, synData, tUC, tUP):
    logger = logging.getLogger("getinversion."+chromo)

    class inversion:
        def __init__(self, i):
            self.profit = i[1]
            self.neighbours = [i[2], i[3]]
            self.invPos = i[0]

    invertedCoordsOri = coords.loc[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == -1)]

    if len(invertedCoordsOri) == 0:
        return(invertedCoordsOri, [],[],invertedCoordsOri,[],[])

    invertedCoords = invertedCoordsOri.copy()
    maxCoords = np.max(np.max(invertedCoords[["bStart","bEnd"]]))

    invertedCoords.bStart = maxCoords + 1 - invertedCoords.bStart
    invertedCoords.bEnd = maxCoords + 1 - invertedCoords.bEnd

    nrow = pd.Series(range(invertedCoords.shape[0]))

    if len(invertedCoordsOri) > 0:
        invTree = pd.DataFrame(apply_TS(invertedCoords.aStart.values,invertedCoords.aEnd.values,invertedCoords.bStart.values,invertedCoords.bEnd.values, threshold), index = range(len(invertedCoords)), columns = invertedCoords.index.values)
    else:
        invTree = pd.DataFrame([], index = range(len(invertedCoords)), columns = invertedCoords.index.values)

    logger.debug("found inv Tree " + chromo)

    #######################################################################
    ###### Create list of inverted alignments
    #######################################################################
    invTree = pd.DataFrame(apply_TS(invertedCoords.aStart.values,invertedCoords.aEnd.values,invertedCoords.bStart.values,invertedCoords.bEnd.values, threshold), index = range(len(invertedCoords)), columns = invertedCoords.index.values)
    invblocks = getInvBlocks(invTree, invertedCoordsOri)
    logger.debug("found inv blocks " + chromo)

    #########################################################################
    ###### Finding profitable inversions (group of inverted blocks)
    #########################################################################

    neighbourSyn = getNeighbourSyn(invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, invertedCoordsOri.index.values, invertedCoordsOri.bDir.values, synData.aStart.values, synData.aEnd.values, synData.bStart.values, synData.bEnd.values, synData.index.values, synData.bDir.values, threshold)

    logger.debug("found neighbours " + chromo)

    synBlockScore = [(i.aLen + i.bLen)*i.iden for index, i in synData.iterrows()]

    ##invPos are 0-indexed positions of inverted alignments in the invertedCoordsOri object
    # profitable = [inversion(cost[i][j], revenue[i][j],
    #                          getNeighbours(neighbourSyn, shortest[i][j]),shortest[i][j])
    #                          for i in range(len(profit)) for j in range(len(profit[i]))\
    #                              if profit[i][j] > (0.1*cost[i][j])]     ##Select only those inversions for which the profit is more than  10% of the cost

    profitable = [inversion(i) for i in getProfitable(invblocks, invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, invertedCoordsOri.iden.values.astype('float32'), neighbourSyn, np.array(synBlockScore, dtype = 'float32'), synData.aStart.values, synData.aEnd.values, tUC, tUP)]

    logger.debug("found profitable " + chromo)

    del(invblocks, invTree, neighbourSyn, synBlockScore)
    collect()
    #####################################################################
    #### Find optimal set of inversions from all profitable inversions
    #####################################################################
    if len(profitable) > 0:
        bestInvPath = invPath({i:profitable[i].invPos for i in range(len(profitable))}, np.array([i.neighbours for i in profitable]), np.array([i.profit for i in profitable], dtype='float32'), invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, threshold)
    else:
        bestInvPath = []

    logger.debug("found bestInvPath " + chromo)

    invBlocksIndex = unlist([profitable[_i].invPos for _i in bestInvPath])
    invData = invertedCoordsOri.iloc[invBlocksIndex]

    badSyn = []
    synInInv = []
    for _i in bestInvPath:
        invNeighbour = profitable[_i].neighbours
#        synInInv = list(range(invNeighbour[0]+1, invNeighbour[1]))
        invPos = profitable[_i].invPos
        invCoord = [invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]]
        for _j in range(invNeighbour[0]+1, invNeighbour[1]):
            sd = synData.iloc[_j][["aStart","aEnd","bStart","bEnd"]]
            if (invCoord[0] - sd[0] < threshold) and (sd[1] - invCoord[1] < threshold) and (invCoord[2] - sd[2] < threshold) and (sd[3] - invCoord[2] < threshold):
                synInInv.append(_j)
            else:
                badSyn.append(_j)
    return(invertedCoordsOri, profitable, bestInvPath,invData, synInInv, badSyn)

# import numpy as np
# from syri.bin.func.myUsefulFunctions import *
# import sys
# import time
# from igraph import *
# from collections import Counter, deque, defaultdict
# from scipy.stats import *
# from datetime import datetime, date
# import pandas as pd
# from multiprocessing import Pool
# from functools import partial
# import os
# from gc import collect
# from Bio.SeqIO import parse
# import logging
# import psutil
# from syri.pyxFiles.synsearchFunctions import apply_TS, alignmentBlock
# from syri.pyxFiles.function cimport getAllLongestPaths, getConnectivityGraph
#
#
# cimport numpy as np
# cimport cython
#
# np.random.seed(1)
#
#
#
# cpdef getInvBlocks(invTree, invertedCoordsOri):
#     cdef int nrow, i, child
#     nrow = invTree.shape[0]
#     invBlocks = [alignmentBlock(i, np.where(invTree.iloc[i,] == True)[0], invertedCoordsOri.iloc[i]) for i in range(nrow)]
#
#     for block in invBlocks:
#         i = 0
#         while(i < len(block.children)):
#             block.children = list(set(block.children) - set(invBlocks[block.children[i]].children))
#             i+=1
#         block.children.sort()
#
#         for child in block.children:
#             invBlocks[child].addParent(block.id)
#     return(invBlocks)
#
#
# cpdef list getShortest(invBlocks):
#     cdef:
#         shortest = deque()
#         int i
#         list j = list(range(len(invBlocks)))
#     invG = getConnectivityGraph(invBlocks)
#     source = np.array(invG.es['source'], dtype = np.int32)
#     target = np.array(invG.es['target'], dtype = np.int32)
#     weight = np.array(invG.es['weight'], dtype = np.float32)
#     n = len(invG.vs.indices)
#     for i in j:
#         if i%100 == 0:
#             print(i, str(datetime.now()))
#         shortest.append(getAllLongestPaths(n,i,j,source,target,weight))
#     short = [list(s) for s in shortest]
#     return short
#
#
# cpdef list getRevenue(invBlocks, shortest, np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.ndarray iDen):
#     cdef:
#         list revenue,i, values, startA, endA, startB, endB, iden
#         np.ndarray[np.int32_t] j
#         np.int32_t k
#         Py_ssize_t l
#     revenue = []
#     for i in shortest:
#         values = []
#         for j in i:
#             if len(j) == 1:
#                 values.append(invBlocks[j[0]].score)
#             else:
#                 score = 0
#                 startA = [aStart[j[0]]]
#                 endA = [aEnd[j[0]]]
#                 startB = [bEnd[j[0]]]
#                 endB = [bStart[j[0]]]
#                 iden = [iDen[j[0]]]
#                 for k in j[1:]:
#                     isMore = True if iDen[k] > iden[-1] else False
#                     if aStart[k] < endA[-1]:
#                         if isMore:
#                             endA[-1] = aStart[k]
#                             startA.append(aStart[k])
#                             endA.append(aEnd[k])
#                         else:
#                             startA.append(endA[-1])
#                             endA.append(aEnd[k])
#                     else:
#                         startA.append(aStart[k])
#                         endA.append(aEnd[k])
#
#                     if bStart[k] > startB[-1]:
#                         if isMore:
#                             startB[-1] = bStart[k]
#                             startB.append(bEnd[k])
#                             endB.append(bStart[k])
#                         else:
#                             endB.append(startB[-1])
#                             startB.append(bEnd[k])
#                     else:
#                         startB.append(bEnd[k])
#                         endB.append(bStart[k])
#                     iden.append(iDen[k])
#                 if len(startA) == len(endA) == len(startB) == len(endB) == len(iden):
#                     for l in range(len(iden)):
#                         score += iden[l]*((endA[l] - startA[l]) + (endB[l] - startB[l]))
#                 values.append(score)
#         revenue = revenue + [values]
#     return(revenue)
#
#
# cpdef dict getNeighbourSyn(np.ndarray aStartInv, np.ndarray aEndInv, np.ndarray bStartInv, np.ndarray bEndInv, np.ndarray indexInv, np.ndarray bDirInv, np.ndarray aStartSyn, np.ndarray aEndSyn, np.ndarray bStartSyn, np.ndarray bEndSyn, np.ndarray indexSyn, np.ndarray bDirSyn, int threshold):
#
#     cdef:
#         cdef Py_ssize_t i, j, index
#         dict neighbourSyn = dict()
#         int upBlock, downBlock
#         list upSyn, downSyn
#     for i in range(len(indexInv)):
#         index = indexInv[i]
#         upSyn = np.where(indexSyn < index)[0].tolist()
#         downSyn = np.where(indexSyn > index)[0].tolist()
#
#         upBlock  = -1
#         downBlock = len(indexSyn)
#         for j in upSyn[::-1]:
#             if bDirSyn[j] == bDirInv[i]:
#                 if (aStartInv[i] - aStartSyn[j]) > threshold and (aEndInv[i] - aEndSyn[j]) > threshold and (bStartInv[i] - bStartSyn[j]) > threshold and (bEndInv[i] - bEndSyn[j]) > threshold:
#                     upBlock = j
#                     break
#             else:
#                 if (aStartInv[i] - aStartSyn[j]) > threshold and (aEndInv[i] - aEndSyn[j]) > threshold and (bEndInv[i] - bStartSyn[j]) > threshold and (bStartInv[i] - bEndSyn[j]) > threshold:
#                     upBlock = j
#                     break
#
#
#         for j in downSyn:
#             if bDirSyn[j] == bDirInv[i]:
#                 if (aStartSyn[j] - aStartInv[i]) > threshold and (aEndSyn[j] - aEndInv[i]) > threshold and (bStartSyn[j] - bStartInv[i]) > threshold and (bEndSyn[j] - bEndInv[i]) > threshold:
#                     downBlock = j
#                     break
#             else:
#                 if (aStartSyn[j] - aStartInv[i]) > threshold and (aEndSyn[j] - aEndInv[i]) > threshold and (bStartSyn[j] - bEndInv[i]) > threshold and (bEndSyn[j] - bStartInv[i]) > threshold:
#                     downBlock = j
#                     break
#         neighbourSyn[i] = [upBlock, downBlock]
#     return(neighbourSyn)
#
#
# cpdef list getCost(list synPath, list shortest, dict neighbourSyn, list synBlockScore, synData, invertedCoordsOri):
#     cdef:
#         list cost, i, values
#         int leftSyn, rightSyn, leftEnd, rightEnd, overlapLength
#         double syncost
#         np.ndarray[np.int32_t] j
#     cost = []
#     synLength = len(synPath)
#     for i in shortest:
#         values = []
#         for j in i:
#             leftSyn, rightSyn = getNeighbours(neighbourSyn, j)
#             synCost = sum([synBlockScore[synIndex] for synIndex in range(leftSyn+1,rightSyn)])
#             leftEnd = synData.iat[leftSyn, 1] if leftSyn > -1 else 0
#             rightEnd = synData.iat[rightSyn,0] if rightSyn < synLength else invertedCoordsOri.iat[j[-1],1]
#             if rightEnd - leftEnd > 1000:
#                 values.append(synCost)
#             else:
#                 overlapLength = (leftEnd - invertedCoordsOri.iat[j[0], 0]) + (invertedCoordsOri.iat[j[-1],1] - rightEnd)
#                 if overlapLength > ((rightEnd - leftEnd)/2):
#                     values.append(synCost + 10000000000000)
#                 else:
#                     values.append(synCost)
#         cost = cost + [values]
#     return(cost)
#
# def getNeighbours(neighbourSyn, j):
#     return(min(neighbourSyn[j[0]]+neighbourSyn[j[-1]]), max(neighbourSyn[j[0]]+neighbourSyn[j[-1]]))
#
#
# def getInversions(coords,chromo, threshold, synData, synPath):
#     logger = logging.getLogger("getinversion."+chromo)
#
#     class inversion:
#         def __init__(self, cost, revenue, neighbours, invPos):
#             self.cost = cost
#             self.revenue = revenue
#             self.profit = revenue - cost
#             self.neighbours = list(neighbours)
#             self.invPos = invPos
#
#     invertedCoordsOri = coords.loc[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == -1)]
#
#     if len(invertedCoordsOri) == 0:
#         return(invertedCoordsOri, [],[],invertedCoordsOri,[],[])
#
#     invertedCoords = invertedCoordsOri.copy()
#     maxCoords = np.max(np.max(invertedCoords[["bStart","bEnd"]]))
#
#     invertedCoords.bStart = maxCoords + 1 - invertedCoords.bStart
#     invertedCoords.bEnd = maxCoords + 1 - invertedCoords.bEnd
#
#     nrow = pd.Series(range(invertedCoords.shape[0]))
#
#     if len(invertedCoordsOri) > 0:
#         invTree = pd.DataFrame(apply_TS(invertedCoords.aStart.values,invertedCoords.aEnd.values,invertedCoords.bStart.values,invertedCoords.bEnd.values, threshold), index = range(len(invertedCoords)), columns = invertedCoords.index.values)
#     else:
#         invTree = pd.DataFrame([], index = range(len(invertedCoords)), columns = invertedCoords.index.values)
#
#     logger.debug("found inv Tree " + chromo)
#
#     #######################################################################
#     ###### Create list of inverted alignments
#     #######################################################################
#
#     invblocks = getInvBlocks(invTree, invertedCoordsOri)
#     logger.debug("found inv blocks " + chromo)
#
#     #########################################################################
#     ###### Finding profitable inversions (group of inverted blocks)
#     #########################################################################
#
#     shortest = getShortest(invblocks)
#     logger.debug("found shortest " + chromo )
#
# #    revenue = getRevenue(invBlocks, shortest, invertedCoordsOri)
#
#     revenue = getRevenue(invblocks, shortest, invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, invertedCoordsOri.iden.values)
#     logger.debug("found revenue " + chromo)
#
#     ## Get syntenic neighbouring blocks of inversions
#
#
# #    neighbourSyn = getNeighbourSyn(invertedCoordsOri, synData, threshold)
#
#     neighbourSyn = getNeighbourSyn(invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, invertedCoordsOri.index.values, invertedCoordsOri.bDir.values, synData.aStart.values, synData.aEnd.values, synData.bStart.values, synData.bEnd.values, synData.index.values, synData.bDir.values, threshold)
#
#     logger.debug("found neighbours " + chromo)
#
#     synBlockScore = [(i.aLen + i.bLen)*i.iden for index, i in synData.iterrows()]
#
#     ## Calculate cost adding an inversion, i.e sum of all synblocks which need to be removed to accomodate teh synblocks
#     cost = getCost(synPath, shortest, neighbourSyn, synBlockScore, synData, invertedCoordsOri)
#     logger.debug("found cost " + chromo)
#
#     ## Calculate profit (or loss) associated with the addition of an inversion
#     profit = []
#     for i in range(len(revenue)):
#         profit = profit + [[revenue[i][j] - cost[i][j] for j in range(len(revenue[i]))]]
#     logger.debug("found profit " + chromo)
#
#     ## Create list of all profitable inversions
#
#     ##invPos are 0-indexed positions of inverted alignments in the invertedCoordsOri object
#     profitable = [inversion(cost[i][j], revenue[i][j],
#                              getNeighbours(neighbourSyn, shortest[i][j]),shortest[i][j])
#                              for i in range(len(profit)) for j in range(len(profit[i]))\
#                                  if profit[i][j] > (0.1*cost[i][j])]     ##Select only those inversions for which the profit is more than  10% of the cost
#     logger.debug("found profitable " + chromo)
#
#     del(invblocks, revenue, neighbourSyn, shortest, synBlockScore)
#     collect()
#     #####################################################################
#     #### Find optimal set of inversions from all profitable inversions
#     #####################################################################
#     profitInvs = [p.profit for p in profitable]
#
#     if len(profitInvs) > 0:
#         lp = len(profitable)
#         iAStart = deque()
#         iAEnd = deque()
#         iBStart = deque()
#         iBEnd = deque()
#         for i in profitable:
#             iAStart.append(invertedCoordsOri.iat[i.invPos[0], 0])
#             iAEnd.append(invertedCoordsOri.iat[i.invPos[-1], 1])
#             iBStart.append(invertedCoordsOri.iat[i.invPos[-1], 3])
#             iBEnd.append(invertedCoordsOri.iat[i.invPos[0], 2])
#
#         iAStart = np.array(iAStart)
#         iAEnd = np.array(iAEnd)
#         iBStart = np.array(iBStart)
#         iBEnd = np.array(iBEnd)
#
#         scores = np.array([i.profit for i in profitable], dtype= int)
#         parents = np.array([-1]*lp, dtype = int)
#         totscore = scores.copy()
#         for i in range(lp):
#             nonOverlapA = np.where(iAStart > (iAEnd[i] - threshold))[0]
#             nonOverlapB = np.where(iBStart > (iBEnd[i] - threshold))[0]
#             childNodes = np.intersect1d(nonOverlapA, nonOverlapB, assume_unique=True) #.astype("uint32") + 1               ## two inversions can co-exist only if the overlap between them is less than threshold on both genomes
#             chIndex =  np.where(scores[childNodes] + totscore[i] > totscore[childNodes])[0]
#             totscore[childNodes[chIndex]] = scores[childNodes[chIndex]] + totscore[i]
#             parents[childNodes[chIndex]] = i
#
#         maxid = totscore.argmax()
#         bestInvPath = deque([maxid])
#         while parents[i] != -1:
#             bestInvPath.append(parents[i])
#             i = parents[i]
#         bestInvPath = list(bestInvPath)[::-1]
#
#     else:
#         bestInvPath = []
#
#
#     logger.debug("found bestInvPath " + chromo)
#
#     invBlocksIndex = unlist([profitable[_i].invPos for _i in bestInvPath])
#     invData = invertedCoordsOri.iloc[invBlocksIndex]
#
#     badSyn = []
#     synInInv = []
#     for _i in bestInvPath:
#         invNeighbour = profitable[_i].neighbours
# #        synInInv = list(range(invNeighbour[0]+1, invNeighbour[1]))
#         invPos = profitable[_i].invPos
#         invCoord = [invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]]
#         for _j in range(invNeighbour[0]+1, invNeighbour[1]):
#             sd = synData.iloc[_j][["aStart","aEnd","bStart","bEnd"]]
#             if (invCoord[0] - sd[0] < threshold) and (sd[1] - invCoord[1] < threshold) and (invCoord[2] - sd[2] < threshold) and (sd[3] - invCoord[2] < threshold):
#                 synInInv.append(_j)
#             else:
#                 badSyn.append(_j)
#
#     return(invertedCoordsOri, profitable, bestInvPath,invData, synInInv, badSyn)
