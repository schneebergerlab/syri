import numpy as np
from syri.bin.func.myUsefulFunctions import *
import sys
import time
from igraph import *
from collections import Counter, deque, defaultdict
from scipy.stats import *
from datetime import datetime, date
import pandas as pd
from multiprocessing import Pool
from functools import partial
import os
from gc import collect
from Bio.SeqIO import parse
import logging
import psutil
from syri.pyxFiles.synsearchFunctions import apply_TS, alignmentBlock
from syri.pyxFiles.function cimport getAllLongestPaths, getConnectivityGraph


cimport numpy as np
cimport cython

np.random.seed(1)



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


cpdef list getShortest(invBlocks):
    cdef:
        shortest = deque()
        int i
        list j = list(range(len(invBlocks)))
    invG = getConnectivityGraph(invBlocks)
    source = np.array(invG.es['source'], dtype = np.int32)
    target = np.array(invG.es['target'], dtype = np.int32)
    weight = np.array(invG.es['weight'], dtype = np.float32)
    n = len(invG.vs.indices)
    for i in j:
        shortest.append(getAllLongestPaths(n,i,j,source,target,weight))
    short = [list(s) for s in shortest]
    return short


cpdef list getRevenue(invBlocks, shortest, np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.ndarray iDen):
    cdef:
        list revenue,i, values, startA, endA, startB, endB, iden
        np.ndarray[np.int32_t] j
        np.int32_t k
        Py_ssize_t l
    revenue = []
    for i in shortest:
        values = []
        for j in i:
            if len(j) == 1:
                values.append(invBlocks[j[0]].score)
            else:
                score = 0
                startA = [aStart[j[0]]]
                endA = [aEnd[j[0]]]
                startB = [bEnd[j[0]]]
                endB = [bStart[j[0]]]
                iden = [iDen[j[0]]]
                for k in j[1:]:
                    isMore = True if iDen[k] > iden[-1] else False
                    if aStart[k] < endA[-1]:
                        if isMore:
                            endA[-1] = aStart[k]
                            startA.append(aStart[k])
                            endA.append(aEnd[k])
                        else:
                            startA.append(endA[-1])
                            endA.append(aEnd[k])
                    else:
                        startA.append(aStart[k])
                        endA.append(aEnd[k])

                    if bStart[k] > startB[-1]:
                        if isMore:
                            startB[-1] = bStart[k]
                            startB.append(bEnd[k])
                            endB.append(bStart[k])
                        else:
                            endB.append(startB[-1])
                            startB.append(bEnd[k])
                    else:
                        startB.append(bEnd[k])
                        endB.append(bStart[k])
                    iden.append(iDen[k])
                if len(startA) == len(endA) == len(startB) == len(endB) == len(iden):
                    for l in range(len(iden)):
                        score += iden[l]*((endA[l] - startA[l]) + (endB[l] - startB[l]))
                values.append(score)
        revenue = revenue + [values]
    return(revenue)


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


cpdef list getCost(list synPath, list shortest, dict neighbourSyn, list synBlockScore, synData, invertedCoordsOri):
    cdef:
        list cost, i, values
        int leftSyn, rightSyn, leftEnd, rightEnd, overlapLength
        double syncost
        np.ndarray[np.int32_t] j
    cost = []
    synLength = len(synPath)
    for i in shortest:
        values = []
        for j in i:
            leftSyn, rightSyn = getNeighbours(neighbourSyn, j)
            synCost = sum([synBlockScore[synIndex] for synIndex in range(leftSyn+1,rightSyn)])
            leftEnd = synData.iat[leftSyn, 1] if leftSyn > -1 else 0
            rightEnd = synData.iat[rightSyn,0] if rightSyn < synLength else invertedCoordsOri.iat[j[-1],1]
            if rightEnd - leftEnd > 1000:
                values.append(synCost)
            else:
                overlapLength = (leftEnd - invertedCoordsOri.iat[j[0], 0]) + (invertedCoordsOri.iat[j[-1],1] - rightEnd)
                if overlapLength > ((rightEnd - leftEnd)/2):
                    values.append(synCost + 10000000000000)
                else:
                    values.append(synCost)
        cost = cost + [values]
    return(cost)

def getNeighbours(neighbourSyn, j):
    return(min(neighbourSyn[j[0]]+neighbourSyn[j[-1]]), max(neighbourSyn[j[0]]+neighbourSyn[j[-1]]))


def getInversions(coords,chromo, threshold, synData, synPath):
    logger = logging.getLogger("getinversion."+chromo)

    class inversion:
        def __init__(self, cost, revenue, neighbours, invPos):
            self.cost = cost
            self.revenue = revenue
            self.profit = revenue - cost
            self.neighbours = list(neighbours)
            self.invPos = invPos

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

    invblocks = getInvBlocks(invTree, invertedCoordsOri)
    logger.debug("found inv blocks " + chromo)

    #########################################################################
    ###### Finding profitable inversions (group of inverted blocks)
    #########################################################################

    shortest = getShortest(invblocks)
    logger.debug("found shortest " + chromo )

#    revenue = getRevenue(invBlocks, shortest, invertedCoordsOri)

    revenue = getRevenue(invblocks, shortest, invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, invertedCoordsOri.iden.values)
    logger.debug("found revenue " + chromo)

    ## Get syntenic neighbouring blocks of inversions


#    neighbourSyn = getNeighbourSyn(invertedCoordsOri, synData, threshold)

    neighbourSyn = getNeighbourSyn(invertedCoordsOri.aStart.values, invertedCoordsOri.aEnd.values, invertedCoordsOri.bStart.values, invertedCoordsOri.bEnd.values, invertedCoordsOri.index.values, invertedCoordsOri.bDir.values, synData.aStart.values, synData.aEnd.values, synData.bStart.values, synData.bEnd.values, synData.index.values, synData.bDir.values, threshold)

    logger.debug("found neighbours " + chromo)

    synBlockScore = [(i.aLen + i.bLen)*i.iden for index, i in synData.iterrows()]

    ## Calculate cost adding an inversion, i.e sum of all synblocks which need to be removed to accomodate teh synblocks
    cost = getCost(synPath, shortest, neighbourSyn, synBlockScore, synData, invertedCoordsOri)
    logger.debug("found cost " + chromo)

    ## Calculate profit (or loss) associated with the addition of an inversion
    profit = []
    for i in range(len(revenue)):
        profit = profit + [[revenue[i][j] - cost[i][j] for j in range(len(revenue[i]))]]
    logger.debug("found profit " + chromo)

    ## Create list of all profitable inversions

    ##invPos are 0-indexed positions of inverted alignments in the invertedCoordsOri object
    profitable = [inversion(cost[i][j], revenue[i][j],
                             getNeighbours(neighbourSyn, shortest[i][j]),shortest[i][j])
                             for i in range(len(profit)) for j in range(len(profit[i]))\
                                 if profit[i][j] > (0.1*cost[i][j])]     ##Select only those inversions for which the profit is more than  10% of the cost
    logger.debug("found profitable " + chromo)

    del(invblocks, revenue, neighbourSyn, shortest, synBlockScore)
    collect()
    #####################################################################
    #### Find optimal set of inversions from all profitable inversions
    #####################################################################
    profitInvs = [p.profit for p in profitable]

    if len(profitInvs) > 0:
        lp = len(profitable)
        iAStart = deque()
        iAEnd = deque()
        iBStart = deque()
        iBEnd = deque()
        for i in profitable:
            iAStart.append(invertedCoordsOri.iat[i.invPos[0], 0])
            iAEnd.append(invertedCoordsOri.iat[i.invPos[-1], 1])
            iBStart.append(invertedCoordsOri.iat[i.invPos[-1], 3])
            iBEnd.append(invertedCoordsOri.iat[i.invPos[0], 2])

        iAStart = np.array(iAStart)
        iAEnd = np.array(iAEnd)
        iBStart = np.array(iBStart)
        iBEnd = np.array(iBEnd)

        scores = np.array([i.profit for i in profitable], dtype= int)
        parents = np.array([-1]*lp, dtype = int)
        totscore = scores.copy()
        for i in range(lp):
            nonOverlapA = np.where(iAStart > (iAEnd[i] - threshold))[0]
            nonOverlapB = np.where(iBStart > (iBEnd[i] - threshold))[0]
            childNodes = np.intersect1d(nonOverlapA, nonOverlapB, assume_unique=True) #.astype("uint32") + 1               ## two inversions can co-exist only if the overlap between them is less than threshold on both genomes
            chIndex =  np.where(scores[childNodes] + totscore[i] > totscore[childNodes])[0]
            totscore[childNodes[chIndex]] = scores[childNodes[chIndex]] + totscore[i]
            parents[childNodes[chIndex]] = i

        maxid = totscore.argmax()
        bestInvPath = deque([maxid])
        while parents[i] != -1:
            bestInvPath.append(parents[i])
            i = parents[i]
        bestInvPath = list(bestInvPath)[::-1]

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
