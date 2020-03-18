%%cython
# distutils: language = c++
import numpy as np
from igraph import Graph
from collections import deque, Counter
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
#from syri.pyxFiles.function cimport getConnectivityGraph
from cython.operator cimport dereference as deref, preincrement as inc


cimport numpy as np

np.random.seed(1)


cpdef inline getConnectivityGraph(blocksList):
    outOG = Graph().as_directed()
    outOG.add_vertices(len(blocksList))
    if len(blocksList) == 0:
        return outOG

    ## Add edges and edge weight
    edgeList = deque()
    esWeight = deque()
    sourceList = deque()
    targetList = deque()
    for i in blocksList:
        if len(i.children) > 0:
            edgeList.extend(list(zip([i.id]*len(i.children), i.children)))
            esWeight.extend([-i.score]*len(i.children))
            sourceList.extend([i.id]*len(i.children))
            targetList.extend(i.children)
    outOG.add_edges(list(edgeList))
    outOG.es["weight"] = list(esWeight)
    outOG.es["source"] = list(sourceList)
    outOG.es["target"] = list(targetList)
    return outOG
    

cpdef getProfitable_test(invblocks, long[:] aStart, long[:] aEnd, long[:] bStart, long[:] bEnd, float[:] iDen, cpp_map[int, cpp_vec[long]] neighbourSyn, float[:] synBlockScore, long[:] aStartSyn, long[:] aEndSyn, long tUC, float tUP, long threshold, brk = -1):
    cdef:
        long                            i, j, k, l, current, count
        long                            n_topo, n_edges, n_syn
        long                            leftSyn, rightSyn, leftEnd, rightEnd, overlapLength
        long                            maxid
        int                             lp, bothuni
        float                           maxscore, synscore
        float                           w, revenue, cost
        Py_ssize_t                      index
        long[:]                         n = np.array(range(len(invblocks)), dtype=np.int)
        long[:]                         topo
        long[:]                         source, target
        long[:]                         pred
        long[:]                         st_list, end_list, stb_list, endb_list, parents
        float[:]                        profit_list, totscore
        float[:]                        dist
        float[:]                        weight
        long[:,:]                       nsynmap = np.zeros((neighbourSyn.size(), 2), dtype='int64')                  # neighbours of inverted alignments
        cpp_map[long, cpp_deq[long]]    nodepath
        cpp_deq[long]                   path, r_path
        cpp_deq[long]                   startA, endA, startB, endB
        cpp_deq[float]                  iden
        bool_t                          isMore
        # Variable for invPath identification
        cpp_deq[long]                   st, end, stb, endb
        cpp_deq[float]                  profit
        cpp_deq[long]                   stsyn, endsyn
    
    print('starting')
    n_syn = len(synBlockScore)
    invG = getConnectivityGraph(invblocks)
    out = deque()

    # get neighbours of inverted alignments
    for i in range(<Py_ssize_t> neighbourSyn.size()):
        nsynmap[i, 0] = neighbourSyn[i][0]
        nsynmap[i, 1] = neighbourSyn[i][1]
    print('1')
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
    print('2')
    # Get order in which the edges need to be transversed
    source = np.zeros_like(invG.es['source'], dtype=int)
    target = np.zeros_like(invG.es['source'], dtype=int)
    weight = np.zeros_like(invG.es['source'], dtype=np.float32)

    index = 0
    garb = invG.get_adjlist('OUT')
    for i in range(n_topo):
        if len(garb[topo[i]]) == 0:
            continue
        w = invG.es[invG.get_eid(topo[i], garb[topo[i]][0])]['weight']
        for j in range(len(garb[topo[i]])):
            source[index] = topo[i]
            target[index] = garb[topo[i]][j]
            weight[index] = w
            index += 1
    n_edges = len(source)
    
    print('3')
    outlist = invG.get_adjlist('OUT')
    inlist = invG.get_adjlist('IN')
    unistart = deque()
    uniend = deque()
    for i in n:
        if len(outlist[i]) == 1:
            if len(inlist[outlist[i][0]]) == 1 and inlist[outlist[i][0]][0]==i:
                stsyn.clear()
                endsyn.clear()
                bothuni = 0
                
                for j in range(nsynmap[i][0]+1,nsynmap[i][1]):
                    if j <= nsynmap[outlist[i][0]][0] or j >= nsynmap[outlist[i][0]][1]:
                        stsyn.push_back(j)
                for j in range(nsynmap[outlist[i][0]][0]+1, nsynmap[outlist[i][0]][1]):
                    if j <= nsynmap[i][0] or j >= nsynmap[i][1]:
                        endsyn.push_back(j)
                
                if stsyn.size()==0 and endsyn.size()==0:
                    print(i, outlist[i], inlist[outlist[i][0]] )
                    unistart.append(i)
                    uniend.append(outlist[i][0])
                
                ## get start node score
                synscore = np.sum([synBlockScore[j] for j in range(<Py_ssize_t> stsyn.size())])
                if (iDen[i]*(aEnd[i] - aStart[i] + 1 + bStart[i] - bEnd[i] + 1)) > synscore:
                    bothuni += 1
                else:
                    continue
                ## get end node score
                synscore = np.sum([synBlockScore[j] for j in range(<Py_ssize_t> endsyn.size())])
                if (iDen[outlist[i][0]]*(aEnd[outlist[i][0]] - aStart[outlist[i][0]] + 1 + bStart[outlist[i][0]] - bEnd[outlist[i][0]] + 1)) > synscore:
                    bothuni += 1
                else:
                    continue
                    
                if bothuni == 2:
                    #print(i, outlist[i], inlist[outlist[i][0]] )
                    unistart.append(i)
                    uniend.append(outlist[i][0])
        
    print('unistart size: ', len(unistart))
    print('uniend size: ', len(uniend))
    unistart = set(list(unistart))
    uniend = set(list(uniend))

    # Find shortest path to all other nodes from each node
    for i in n:
        if i in uniend:
            continue
        nodepath.clear()
        pred = np.array([-1]* <Py_ssize_t> len(n), dtype = np.int)
        dist = np.array([np.float32('inf')]*  <Py_ssize_t> len(n), dtype = np.float32)
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
            if j in unistart:
                continue
            # Find all connected paths which are profitable
            if dist[topo[j]] != float("inf"):
                current = topo[j]
                path.clear()
                while current!=i:
                    if nodepath.count(current) > 0:
                        for index in range(<Py_ssize_t> nodepath[current].size()):
                            path.push_back(nodepath[current][index])
                        # path.extend(nodepath[current].copy())
                        break
                    path.push_back(current)
                    current = pred[current]
                nodepath[topo[j]] = path
                path.push_back(i)
                r_path.clear()

                current = path.size()
                for index in range(<Py_ssize_t> path.size()):
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
                        for k in range(<Py_ssize_t> iden.size()):
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
                    if (rightEnd - leftEnd)/(aEnd[r_path.back()] - aStart[r_path.front()]) < tUP:
                        cost = 10000000000000

                #  Select those candidate inversions for which the score of
                #  adding them would be at least 10% better than the score
                #  of syntenic regions needed to be removed
                if revenue > 1.1*cost:
                    st.push_back(aStart[r_path.front()])
                    end.push_back(aEnd[r_path.back()])
                    stb.push_back(bEnd[r_path.back()])
                    endb.push_back(bStart[r_path.front()])
                    profit.push_back(revenue - cost)
        if i == brk:
            break
    print(st.size())
    
    path.clear()
    lp = st.size()
    totscore = np.array([profit[i] for i in range(<Py_ssize_t> profit.size())], dtype=np.float32)
    st_list  = np.array([st[i] for i in range(<Py_ssize_t> st.size())], dtype=np.int)
    end_list  = np.array([end[i] for i in range(<Py_ssize_t> end.size())], dtype=np.int)
    stb_list  = np.array([stb[i] for i in range(<Py_ssize_t> stb.size())], dtype=np.int)
    endb_list  = np.array([endb[i] for i in range(<Py_ssize_t> endb.size())], dtype=np.int)

    parents = np.array([-1]*lp, dtype = 'int')

    for i in range(lp):
        if i%100==0:
            print(i, str(datetime.now()))
        for j in range(lp-1,i,-1):
            if st_list[j] > end_list[i]-threshold:
                if stb_list[j] > endb_list[i] -threshold:
                    if profit[j] + totscore[i] > totscore[j]:
                        totscore[j] = profit[j] + totscore[i]
                        parents[j] = i
            else:
                break
    maxid = -1
    maxscore = -1
    for i in range(lp):
        if totscore[i] > maxscore:
            maxscore = totscore[i]
            maxid = i
    if maxid == -1:
        return out
    path.push_front(maxid)
    while parents[maxid] != -1:
        path.push_front(parents[maxid])
        maxid = parents[maxid]
    goodinvs = set([path[i] for i in range(<Py_ssize_t> path.size())])
    
    #print(goodinvs)
    count = 0
    for i in n:
        if i in uniend:
            continue
        nodepath.clear()
        pred = np.array([-1]* <Py_ssize_t> len(n), dtype = np.int)
        dist = np.array([np.float32('inf')]*  <Py_ssize_t> len(n), dtype = np.float32)
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
            if j in unistart:
                continue
            # Find all connected paths which are profitable
            if dist[topo[j]] != float("inf"):
                current = topo[j]
                path.clear()
                while current!=i:
                    if nodepath.count(current) > 0:
                        for index in range(<Py_ssize_t> nodepath[current].size()):
                            path.push_back(nodepath[current][index])
                        # path.extend(nodepath[current].copy())
                        break
                    path.push_back(current)
                    current = pred[current]
                nodepath[topo[j]] = path
                path.push_back(i)
                r_path.clear()

                current = path.size()
                for index in range(<Py_ssize_t> path.size()):
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
                        for k in range(<Py_ssize_t> iden.size()):
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
                    if (rightEnd - leftEnd)/(aEnd[r_path.back()] - aStart[r_path.front()]) < tUP:
                        cost = 10000000000000

                #  Select those candidate inversions for which the score of
                #  adding them would be at least 10% better than the score
                #  of syntenic regions needed to be removed
                if revenue > 1.1*cost:
                    if count in goodinvs:
                        out.append(([r_path[k] for k in range(current)], revenue - cost, leftSyn, rightSyn))
                    count += 1
        if i == brk:
            return out
    return out


import pickle
with open("/srv/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/tests/trash/syri_test_genomes/analysis3/profitable_arguments.pickle", 'rb') as fin:
    input = pickle.load(fin)
    invblocks = input[0]
    invertedCoordsOri = input[1]
    neighbourSyn = input[2]
    synBlockScore = input[3]
    synData = input[4]
    tUC = input[5]
    tUP = input[6]
    threshold = input[7]

print('running')
getProfitable_test(invblocks, invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, invertedCoordsOri.iden.values.astype('float32'), neighbourSyn, np.array(synBlockScore, dtype = 'float32'), synData.aStart.values, synData.aEnd.values, tUC, tUP, threshold)