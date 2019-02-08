import numpy as np
from syri.bin.func.myUsefulFunctions import *
import sys
import time
from igraph import *
from collections import deque
from datetime import datetime, date
# from scipy.stats import *
import pandas as pd
from multiprocessing import Pool
from functools import partial
from gc import collect
import logging
from syri.pyxFiles.function cimport getAllLongestPaths, getConnectivityGraph
from syri.pyxFiles.function cimport getmeblocks, getOverlapWithSynBlocks
from syri.pyxFiles.synsearchFunctions import alignmentBlock


cimport numpy as np
cimport cython

np.random.seed(1)

class transGroups:
    def __init__(self, leftEnd, rightEnd, index, threshold):
        self.leftEnd = leftEnd
        self.rightEnd = rightEnd
        self.member = [index]
        self.threshold = threshold

    def checkOverlap(self, leftEnd, rightEnd):
        if leftEnd < self.leftEnd:
            print("Blocks must be sorted")
            sys.exit()
        elif rightEnd < self.rightEnd + self.threshold:
            return True
        elif (rightEnd - self.rightEnd) < 0.5*(self.rightEnd - leftEnd):
            return True
        elif (leftEnd - self.leftEnd) < 0.5*(self.rightEnd - leftEnd):
            return True
        return False

    def addMember(self, leftEnd, rightEnd, index):
        self.leftEnd = min(self.leftEnd, leftEnd)
        self.rightEnd = max(self.rightEnd, rightEnd)
        self.member.append(index)

def findOverlappingSynBlocks(inPlaceBlocks, aStart, aEnd, bStart, bEnd):
    aBlocks = list(np.where((inPlaceBlocks.aStart.values < aEnd) & (inPlaceBlocks.aEnd.values > aStart) == True)[0])
    bBlocks = list(np.where((inPlaceBlocks.bStart.values < bEnd) & (inPlaceBlocks.bEnd.values > bStart) == True)[0])
    return(aBlocks, bBlocks)


cpdef np.ndarray getTranslocationScore(np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.ndarray aLen, np.ndarray bLen, np.ndarray translocations):
    """Function to score the proposed translocation block based on the number of
        basepairs it explains and the gaps between alignments of the block
    """
    cdef Py_ssize_t i,j,k, n = len(translocations)
    cdef int l
    cdef np.float blockScore, aScore, bScore, aGap, bGap
    cdef np.ndarray transScores = np.array([-1]*n, dtype = object), blocksScores
    for i in range(n):
        l = len(translocations[i])
        blocksScores = np.array([-1]*l, dtype = object)
        for j in range(l):
            aScore = np.float(aLen[translocations[i][j][0]])
            bScore = np.float(bLen[translocations[i][j][0]])
            aGap = np.float(0)
            bGap = np.float(0)
            if len(translocations[i][j]) > 1:
                for k in range(1, len(translocations[i][j])):
                    aScore += np.float(aLen[translocations[i][j][k]])
                    bScore += np.float(bLen[translocations[i][j][k]])
                    aGap += np.float(max(0, aStart[translocations[i][j][k]] - aEnd[translocations[i][j][k-1]]))
                    bGap += np.float(max(0, bStart[translocations[i][j][k]] - bEnd[translocations[i][j][k-1]]))
                blockScore = min(((aScore - aGap)/aScore),((bScore - bGap)/bScore))
                blocksScores[j] = blockScore
            else:
                blocksScores[j] = 1
        transScores[i] = blocksScores
    return transScores
#
#

cpdef np.ndarray getTranslocationScore_ctx(np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.ndarray aLen, np.ndarray bLen, np.ndarray bDir, np.ndarray translocations):
    """Function to score the proposed translocation block based on the number of
        basepairs it explains and the gaps between alignments of the block
    """
    cdef Py_ssize_t i,j,k, n = len(translocations)
    cdef int l
    cdef np.float blockScore, aScore, bScore, aGap, bGap
    cdef np.ndarray transScores = np.array([-1]*n, dtype = object), blocksScores
    for i in range(n):
        l = len(translocations[i])
        blocksScores = np.array([-1]*l, dtype = object)
        for j in range(l):
            aScore = np.float(aLen[translocations[i][j][0]])
            bScore = np.float(bLen[translocations[i][j][0]])
            aGap = np.float(0)
            bGap = np.float(0)
            if len(translocations[i][j]) > 1:
                for k in range(1, len(translocations[i][j])):
                    aScore += np.float(aLen[translocations[i][j][k]])
                    bScore += np.float(bLen[translocations[i][j][k]])
                    aGap += np.float(max(0, aStart[translocations[i][j][k]] - aEnd[translocations[i][j][k-1]]))
                    if bDir[k] == 1:
                        bGap += np.float(max(0, bStart[translocations[i][j][k]] - bEnd[translocations[i][j][k-1]]))
                    else:
                        bGap += np.float(max(0, bStart[translocations[i][j][k-1]] - bEnd[translocations[i][j][k]]))
                blockScore = min(((aScore - aGap)/aScore),((bScore - bGap)/bScore))
                blocksScores[j] = blockScore
            else:
                blocksScores[j] = 1
        transScores[i] = blocksScores
    return transScores

def getDupGenome(dupData, allTransBlocksData, transClasses):
    dupGenomes = []
    for row in dupData.itertuples(index = True):
        found = False
        tempTransBlock = allTransBlocksData[row.Index]
        if not tempTransBlock.aUni:
            dupGenomes.append("B")
            continue
        elif not tempTransBlock.bUni:
            dupGenomes.append("A")
            continue
        elif tempTransBlock.genomeAUni:
            dupGenomes.append("A")
            continue
        elif tempTransBlock.genomeBUni:
            dupGenomes.append("B")
            continue
        for i in tempTransBlock.meAlist:
            if i in transClasses["translocation"] or i in transClasses["invTranslocation"]:
                found = True
                dupGenomes.append("B")
                break
        if not found:
            dupGenomes.append("A")
    dupData["dupGenomes"] = pd.Series(dupGenomes, index = dupData.index)
    return(dupData)


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

def readAnnoCoords(cwdPath, uniChromo, prefix):
    annoCoords = pd.DataFrame(columns=["aStart","aEnd","bStart","bEnd","aChr","bChr"])
    synData = []
    fin = open(cwdPath+prefix+"synOut.txt","r")
    for line in fin:
        line = line.strip().split("\t")
        if line[0] == "#":
            chromo = line[1]
            continue
        synData.append(list(map(int,line[:4]))+[chromo,chromo])
    fin.close()
    synData = pd.DataFrame(synData,columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"])
    annoCoords = annoCoords.append(synData)

    for i in ["invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt"]:
        data = []
        fin = open(cwdPath+prefix+i,"r")
        for line in fin:
            line = line.strip().split("\t")
            if line[0] == "#":
                data.append(list(map(int,getValues(line,[2,3,6,7]))) + [line[1],line[5]])
        fin.close()
        data = pd.DataFrame(data, columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"], dtype=object)
        annoCoords = annoCoords.append(data)

    annoCoords[["aStart","aEnd","bStart","bEnd"]] = annoCoords[["aStart","aEnd","bStart","bEnd"]].astype("int64")
    annoCoords.sort_values(by = ["bChr","bStart","bEnd","aChr","aStart","aEnd"],inplace = True)
    annoCoords["bIndex"] = range(len(annoCoords))
    annoCoords.sort_values(by = ["aChr","aStart","aEnd","bChr","bStart","bEnd"],inplace = True)
    annoCoords.index = range(len(annoCoords))
    annoCoords["aIndex"] = range(len(annoCoords))
    return(annoCoords)


def findOrderedTranslocations(outOrderedBlocks, orderedBlocks, inPlaceBlocks, threshold, tUC, tUP, ctx = False):
    logger = logging.getLogger("findOrderedTranslocations")
    if not isinstance(ctx, bool):
        print("CTX status must be a boolean")
        sys.exit()
    def makeBlocksList(blocksTree, blocksData):
        _blocksList = [alignmentBlock(_i, np.where(blocksTree.iloc[_i] == True)[0],blocksData.iloc[_i]) for _i in range(blocksTree.shape[0])]
        for _block in _blocksList:
            _i = 0
            while _i < len(_block.children):
                _block.children = list(set(_block.children) - set(_blocksList[_block.children[_i]].children))
                _i+=1
            _block.children.sort()

            for _child in _block.children:
                _blocksList[_child].addParent(_block.id)
        return _blocksList

    def getTransBlocks(transScores, shortestOutOG, orderedBlocks, inPlaceBlocks, threshold, tUC,tUP, ctx):
        """This method filters possible translocation blocks to select those which have a positive gap based score
           (output of `getTransLocationsScore`) and those which dont overlap significantly with the inPlaceBlocks.

           Parameters
           ----------
           transScores: list,
               output of getTransLocationScores, scores of all shortest blocks
           shortTrans: list,
               list of best blocks between everypair of blocks
           orderedBlocks: DataFrame,
               all translocated alignment blocks
           inPlaceBlocks: DataFrame,
               all syntenic alignment blocks
           threshold: int,
               cut-off value.
           Returns
           --------
           outBlocks: list,
               selected shortTrans blocks
        """

        positiveTransScores = [np.where(np.array(i) >= 0)[0] for i in transScores]
        transBlocks = [getValues(shortestOutOG[i], positiveTransScores[i]) for i in range(len(shortestOutOG))]
        transBlocks = [i for j in transBlocks for i in j]
        outBlocks = []
        if not isinstance(ctx,bool):
            print("CTX status must be a boolean")
            sys.exit()
        if not ctx:
            allAlmnt = pd.unique(unlist(transBlocks))
            almntData = {}
            for almnt in allAlmnt:
                almntData[almnt] = {}
                aScore = 0
                bScore = 0
                aStart = orderedBlocks.iat[almnt,0]
                aEnd = orderedBlocks.iat[almnt,1]
                if orderedBlocks.iat[almnt,8] == 1:
                    bStart = orderedBlocks.iat[almnt,2]
                    bEnd = orderedBlocks.iat[almnt,3]
                else:
                    bStart = orderedBlocks.iat[almnt,3]
                    bEnd = orderedBlocks.iat[almnt,2]

                aBlocks, bBlocks = findOverlappingSynBlocks(inPlaceBlocks, aStart, aEnd, bStart, bEnd)

                for aBlock in aBlocks:
                    if inPlaceBlocks.iat[aBlock,0] - aStart < threshold and aEnd - inPlaceBlocks.iat[aBlock,1] < threshold:
                        aStart = aEnd
                        break
                    elif inPlaceBlocks.iat[aBlock,0] < aStart and inPlaceBlocks.iat[aBlock,1] < aEnd:
                        aStart = inPlaceBlocks.iat[aBlock,1]
                    else:
                        aScore += inPlaceBlocks.iat[aBlock,0] - aStart
                        if inPlaceBlocks.iat[aBlock,1] < aEnd:
                            aStart = inPlaceBlocks.iat[aBlock,1]
                        else:
                            aStart = aEnd
                            break
                aScore += aEnd - aStart

                for bBlock in bBlocks:
                    bBlockEnd = inPlaceBlocks.iat[bBlock,3]
                    bBlockStart = inPlaceBlocks.iat[bBlock,2]
                    if bBlockStart - bStart < threshold and bEnd - bBlockEnd < threshold:
                        bStart = bEnd
                        break
                    elif bBlockStart < bStart and bBlockEnd < bEnd:
                        bStart = bBlockEnd
                    else:
                        bScore += bBlockStart - bStart
                        if bBlockEnd < bEnd:
                            bStart = bBlockEnd
                        else:
                            bStart = bEnd
                            break
                bScore += bEnd - bStart

                almntData[almnt]['aLen'] = orderedBlocks.iat[almnt,4]
                almntData[almnt]['bLen'] = orderedBlocks.iat[almnt,5]
                almntData[almnt]['aScore'] = aScore
                almntData[almnt]['bScore'] = bScore

            for block in transBlocks:
                blockAlength = 0
                blockBlength = 0
                blockAUni = 0
                blockBUni = 0
                for almnt in block:
                    blockAlength += almntData[almnt]['aLen']
                    blockBlength += almntData[almnt]['bLen']
                    blockAUni += almntData[almnt]['aScore']
                    blockBUni += almntData[almnt]['bScore']

        #Trans block is selected IFF either the unique region on any genome is larger than 1kb
        # or length of unique region on a genome is larger than 0.5 times the length of
        # the overlapping region on that genome
                if blockAUni > tUC or blockBUni > tUC or blockAUni > tUP*blockAlength or blockBUni > tUP*blockBlength:
                    outBlocks.append(block)
            return(outBlocks)
        ##########
        ## With CTX
        ##########
        if ctx:
            allAlmnt = np.unique(unlist(transBlocks))
            almntData = {}
            for almnt in allAlmnt:
                almntData[almnt] = {}
                aScore = 0
                bScore = 0
                index = orderedBlocks.index.values[almnt]
                aStart = orderedBlocks.at[index,"aStart"]
                aEnd = orderedBlocks.at[index,"aEnd"]
                aChr = orderedBlocks.at[index,"aChr"]
                bStart = orderedBlocks.at[index,"bStart"]
                bEnd = orderedBlocks.at[index,"bEnd"]
                bChr = orderedBlocks.at[index,"bChr"]

                if bEnd < bStart:
                    print("CTX Input: bStart must be less than bEnd")
                    sys.exit()

                aBlocks = list(intersect(np.where(inPlaceBlocks.aStart.values <  aEnd)[0],
                                              np.where(inPlaceBlocks.aEnd.values >  aStart)[0],
                                              np.where(inPlaceBlocks.aChr == aChr)[0]))
                aBlocks = getValues(inPlaceBlocks.index.values,aBlocks)

                for aBlock in aBlocks:
                    if inPlaceBlocks.at[aBlock,"aStart"] - aStart < threshold and aEnd - inPlaceBlocks.at[aBlock,"aEnd"] < threshold:
                        aStart = aEnd
                        break
                    elif inPlaceBlocks.at[aBlock,"aStart"] < aStart and inPlaceBlocks.at[aBlock,"aEnd"] < aEnd:
                        aStart = inPlaceBlocks.at[aBlock,"aEnd"]
                    else:
                        aScore += inPlaceBlocks.at[aBlock,"aStart"] - aStart
                        if inPlaceBlocks.at[aBlock,"aEnd"] < aEnd:
                            aStart = inPlaceBlocks.at[aBlock,"aEnd"]
                        else:
                            aStart = aEnd
                            break
                aScore += aEnd - aStart
                bBlocks = list(intersect(np.where(inPlaceBlocks.bStart.values <  bEnd)[0],
                                              np.where(inPlaceBlocks.bEnd.values >  bStart)[0],
                                              np.where(inPlaceBlocks.bChr == bChr)[0]))
                bBlocks = getValues(inPlaceBlocks.index.values, bBlocks)

                for bBlock in bBlocks:
                    bBlockStart = inPlaceBlocks.at[bBlock,"bStart"]
                    bBlockEnd = inPlaceBlocks.at[bBlock,"bEnd"]
                    if bBlockStart -bStart < threshold and bEnd - bBlockEnd < threshold:
                        bStart = bEnd
                        break
                    elif bBlockStart < bStart and bBlockEnd < bEnd:
                        bStart = bBlockEnd
                    else:
                        bScore += bBlockStart - bStart
                        if bBlockEnd < bEnd:
                            bStart = bBlockEnd
                        else:
                            bStart = bEnd
                            break
                bScore += bEnd - bStart

                almntData[almnt]['aScore'] = aScore
                almntData[almnt]['bScore'] = bScore
                almntData[almnt]['aLen'] = orderedBlocks.at[index,"aLen"]
                almntData[almnt]['bLen'] = orderedBlocks.at[index,"bLen"]

            for block in transBlocks:
                blockAlength = 0
                blockBlength = 0
                blockAUni = 0
                blockBUni = 0
                for almnt in block:
                    blockAlength += almntData[almnt]['aLen']
                    blockBlength += almntData[almnt]['bLen']
                    blockAUni += almntData[almnt]['aScore']
                    blockBUni += almntData[almnt]['bScore']

        #Trans block is selected IFF either the unique region on any genome is larger than tUC
        # or length of unique region on a genome is larger than tUP times the length of
        # the overlapping region on that genome
                if blockAUni > tUC or blockBUni > tUC or blockAUni > tUP*blockAlength or blockBUni > tUP*blockBlength:
                    outBlocks.append(block)
            return(outBlocks)


    orderedBlocksList = makeBlocksList(outOrderedBlocks, orderedBlocks)
    outOG = getConnectivityGraph(orderedBlocksList)
    shortestOutOG = []
    source = np.array(outOG.es['source'], dtype = np.int32)
    target = np.array(outOG.es['target'], dtype = np.int32)
    weight = np.array(outOG.es['weight'], dtype = np.float32)
    for i in range(len(orderedBlocksList)):
        eNode = [i]
        eNode.extend(list(np.where(outOrderedBlocks.iloc[i] == True)[0]))
        shortestOutOG.append(getAllLongestPaths(len(outOG.vs.indices),i,eNode,source, target, weight, "weight"))
    shortestOutOG = np.array(shortestOutOG)
    logger.debug("starting getTranslocationScore")
    if not ctx:
        transScores = getTranslocationScore(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, orderedBlocks.aLen.values, orderedBlocks.bLen.values, shortestOutOG)
    elif ctx:
        transScores = getTranslocationScore_ctx(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, orderedBlocks.aLen.values, orderedBlocks.bLen.values, orderedBlocks.bDir.values, shortestOutOG)
        logger.debug("finished getTranslocationScore")

#    transScores = getTranslocationScore(shortestOutOG, orderedBlocks, ctx)
    transBlocks = getTransBlocks(transScores, shortestOutOG, orderedBlocks, inPlaceBlocks, threshold, tUC, tUP, ctx)
    logger.debug("finished getTransBlocks")

    return transBlocks




cpdef np.ndarray[np.npy_bool, ndim=2] makeBlocksTree_ctx(np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.ndarray bDir, np.ndarray aChr, np.ndarray bChr, np.ndarray index, np.int threshold):
    """Compute whether two alignments can be part of one translation block. For this:
        the alignments should not be separated by any inPlaceBlock on both ends and
        they should be syntenic with respect to each other.
       
       Returns
       --------
       outOrderedBlocks: pandas DataFrame,
           Dataframe of type Object. Lower half is NA, upper half contains whether two
           alignments can be connected (True) or not (False).
    """
    assert(aStart.dtype==np.int and aEnd.dtype==np.int and bStart.dtype==np.int and bEnd.dtype==np.int and bDir.dtype==np.int and aChr.dtype==np.object and bChr.dtype==np.object and index.dtype==np.int)
    cdef Py_ssize_t i,j, n = len(aStart)
    assert(n == len(aEnd) == len(bStart) == len(bEnd) == len(index) == len(bDir) == len(aChr) == len(bChr))

    cdef np.ndarray outOrderedBlocks =  np.array([[np.False_]*n]*n, dtype=np.bool)

    for i in range(n):
        for j in range(i,n):
            if bDir[i] != bDir[j]:
                sys.exit("ERROR: bDir not matching")
            elif aChr[i] != aChr[j] or bChr[i] != bChr[j]:
                continue
            elif bDir[i] == 1:
                if (aStart[j] - aStart[i]) > threshold and (aEnd[j] - aEnd[i]) > threshold and (bStart[j] - bStart[i]) > threshold and (bEnd[j] - bEnd[i]) > threshold:
                    outOrderedBlocks[i][j] = True          #True
                else:
                    continue
            elif bDir[i] == -1:
                if (aStart[j] - aStart[i]) > threshold and (aEnd[j] - aEnd[i]) > threshold and (bStart[i] - bStart[j]) > threshold and (bEnd[i] - bEnd[j]) > threshold:
                    outOrderedBlocks[i][j] = True
                else:
                    continue
            else:
                sys.exit("ERROR: ILLEGAL BDIR VALUE")
    return(outOrderedBlocks)

def getBlocks(orderedBlocks, annoCoords, threshold, tUC, tUP):
    if len(orderedBlocks) == 0:
        return([])
    outOrderedBlocks = pd.DataFrame(makeBlocksTree_ctx(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, orderedBlocks.bDir.values, orderedBlocks.aChr.values, orderedBlocks.bChr.values, orderedBlocks.index.values, threshold))
    transBlocks = findOrderedTranslocations(outOrderedBlocks, orderedBlocks, annoCoords, threshold, tUC, tUP, ctx = True)
    return(transBlocks)


def mergeTransBlocks(transBlocks, orderedBlocks, invTransBlocks, invertedBlocks, ctx = False):
    if not isinstance(ctx,bool):
        print("CTX status must be a boolean")
        sys.exit()
    if not ctx:
        transBlocksData = []
        for i in transBlocks:
            aStart = orderedBlocks.iat[i[0],0]
            aEnd = orderedBlocks.iat[i[-1],1]
            if orderedBlocks.iat[i[0],8] == 1:
                bStart = orderedBlocks.iat[i[0],2]
                bEnd = orderedBlocks.iat[i[-1],3]
                transDir = 1
            else:
                bStart = orderedBlocks.iat[i[-1],3]
                bEnd = orderedBlocks.iat[i[0],2]
                transDir = -1
            transBlocksData.append([aStart, aEnd, bStart, bEnd, transDir])

        for i in invTransBlocks:
            aStart = invertedBlocks.iat[i[0],0]
            aEnd = invertedBlocks.iat[i[-1],1]
            if invertedBlocks.iat[i[0],8] == 1:
                bStart = invertedBlocks.iat[i[0],2]
                bEnd = invertedBlocks.iat[i[-1],3]
                transDir = 1
            else:
                bStart = invertedBlocks.iat[i[-1],3]
                bEnd = invertedBlocks.iat[i[0],2]
                transDir = -1
            transBlocksData.append([aStart, aEnd, bStart, bEnd, transDir])
        transBlocksData = pd.DataFrame(transBlocksData, columns =  ["aStart","aEnd","bStart","bEnd","dir"])
        transBlocksData.index = list(range(len(transBlocks))) + list(range(len(invTransBlocks)))
        transBlocksData.sort_values(["aStart","aEnd","bStart","bEnd"], inplace = True)
        orderedIndex = transBlocksData.index.values
        transBlocksData.index = range(transBlocksData.shape[0])
        return(transBlocksData, orderedIndex)

    if ctx:
        transBlocksData = []
        for i in transBlocks:
            indices = getValues(orderedBlocks.index.values,i)
            aStart = orderedBlocks.at[indices[0],"aStart"]
            aEnd = orderedBlocks.at[indices[-1],"aEnd"]
            bStart = orderedBlocks.at[indices[0],"bStart"]
            bEnd = orderedBlocks.at[indices[-1],"bEnd"]
            aDir = 1
            bDir = 1
            aChr = orderedBlocks.at[indices[0],"aChr"]
            bChr = orderedBlocks.at[indices[0],"bChr"]
            transBlocksData.append([aStart, aEnd, bStart, bEnd, aDir, bDir, aChr, bChr])

        for i in invTransBlocks:
            indices = getValues(invertedBlocks.index.values,i)
            aStart = invertedBlocks.at[indices[0],"aStart"]
            aEnd = invertedBlocks.at[indices[-1],"aEnd"]
            bStart = invertedBlocks.at[indices[-1],"bStart"]
            bEnd = invertedBlocks.at[indices[0],"bEnd"]
            aDir = 1
            bDir = -1
            aChr = invertedBlocks.at[indices[0],"aChr"]
            bChr = invertedBlocks.at[indices[0],"bChr"]
            transBlocksData.append([aStart, aEnd, bStart, bEnd, aDir, bDir, aChr, bChr])

        transBlocksData = pd.DataFrame(transBlocksData, columns =  ["aStart","aEnd","bStart","bEnd","aDir","bDir", "aChr","bChr"])
        transBlocksData.index = list(range(len(transBlocks))) + list(range(len(invTransBlocks)))
        transBlocksData.sort_values(["aChr","aStart","aEnd","bChr","bStart","bEnd"], inplace = True)
        orderedIndex = transBlocksData.index.values
        transBlocksData.index = range(transBlocksData.shape[0])
        return transBlocksData, orderedIndex

def makeTransGroupList(transBlocksData, startC, endC, threshold):
    transBlocksTable = transBlocksData.sort_values([startC,endC])
    indices = transBlocksTable.index.values
    if len(transBlocksData) > 0:
        genomeGroups = [transGroups(transBlocksTable.at[indices[0],startC],\
                                    transBlocksTable.at[indices[0],endC], indices[0], threshold)]
        for i in indices[1:]:
            if transBlocksTable.at[i, startC] > genomeGroups[-1].rightEnd:
                genomeGroups.append(transGroups(transBlocksTable.at[i,startC],\
                                                transBlocksTable.at[i,endC], i, threshold))
            elif genomeGroups[-1].checkOverlap(transBlocksTable.at[i,startC],\
                             transBlocksTable.at[i,endC]):
                genomeGroups[-1].addMember(transBlocksTable.at[i,startC],\
                            transBlocksTable.at[i,endC], i)
            else:
                genomeGroups.append(transGroups(transBlocksTable.at[i,startC],\
                                                transBlocksTable.at[i,endC], i, threshold))
        return genomeGroups
    else:
        return []


def getTransCluster(transGroupIndices, transGenomeAGroups, transGenomeBGroups):
    assert(list(transGroupIndices.keys()) == list(range(len(transGroupIndices))))
    nodeStack = np.zeros(len(transGroupIndices), dtype='uint8')
    visitedTransBlock = np.zeros(len(transGroupIndices), dtype='uint8')
    visitedIndices = deque()
    transCluster = []
    addedAGroups = []
    addedBGroups = []
    count = 0
    for key,value in transGroupIndices.items():
        if visitedTransBlock[key] == 0:
            visitedTransBlock[key]=1
            visitedIndices.append(key)
            newGroup = [key]
            node1 = value[0]
            node2 = value[1]
            addedAGroups.append(value[0])
            addedBGroups.append(value[1])
            nodeStack[transGenomeAGroups[node1].member] = 1
            nodeStack[transGenomeBGroups[node2].member] = 1
            nodeStack[visitedIndices] = 0
#            nodeStack[np.nonzero(visitedTransBlock)[0]] = 0
            while 1 in nodeStack:
                count=count+1
                if count%20000 ==0:
                    print(count, str(datetime.now()))

                newKey = np.where(nodeStack == 1)[0][0]
                if visitedTransBlock[newKey]== 0:
                    visitedTransBlock[newKey] = 1
                    visitedIndices.append(newKey)
                    newGroup.append(newKey)
                    aInd = transGroupIndices[newKey][0]
                    bInd = transGroupIndices[newKey][1]
                    nodeStack[newKey] = 0

                    if aInd not in addedAGroups:
                        nodeStack[transGenomeAGroups[aInd].member] = 1
                        addedAGroups.append(aInd)
                        nodeStack[visitedIndices] = 0
                    if bInd not in addedBGroups:
                        nodeStack[transGenomeBGroups[bInd].member] = 1
                        addedBGroups.append(bInd)
                        nodeStack[visitedIndices] = 0
            newGroup.sort()
            transCluster.append(list(newGroup))
    return(transCluster)



class transBlock:
    def __init__(self, i):
        self.aStart = None
        self.aEnd = None
        self.bStart = None
        self.bEnd = None
        self.dir = None
        self.transBlocksID = i
        self.transClusterIndex = None
        self.status = 0
        self.overlappingInPlaceBlocks = []
        self.aUni = False
        self.bUni = False
        self.transGroupIndices = None
        self.genomeAUni = False
        self.genomeBUni = False

    def checkoverlaps(self, agroups, bgroups):
        a = self.getoverlappingregions(agroups, "a")
        b = self.getoverlappingregions(bgroups, "b")
        self.genomeAUni = True if len(a) == 0 else False
        self.genomeBUni = True if len(b) == 0 else False

    def getoverlappingregions(self, groups, genome):
        if genome=="a":
            reg= groups[self.transGroupIndices[0]].member.copy()
        elif genome =="b":
            reg = groups[self.transGroupIndices[1]].member.copy()
        reg.remove(self.transBlocksID)
        return reg

    def addOrderedData(self, orderedData):
        self.orderedData = orderedData

    def checkOverlapWithSynBlocks(self,inPlaceBlocks, threshold):
        aBlocks, bBlocks = findOverlappingSynBlocks(inPlaceBlocks, self.aStart, self.aEnd, self.bStart, self.bEnd)

        blockAUni = 0
        blockBUni = 0

        start = self.aStart
        end = self.aEnd

        for aBlock in aBlocks:
            if inPlaceBlocks.iat[aBlock,0] - start < threshold and\
                end - inPlaceBlocks.iat[aBlock,1] < threshold:
                start = end
                break
            elif inPlaceBlocks.iat[aBlock,0] < start and inPlaceBlocks.iat[aBlock,1] < end:
                start = inPlaceBlocks.iat[aBlock,1]
            else:
                blockAUni += inPlaceBlocks.iat[aBlock,0] - start
                if inPlaceBlocks.iat[aBlock,1] < end:
                    start = inPlaceBlocks.iat[aBlock,1]
                else:
                    start = end
                    break
        blockAUni += end - start

        start = self.bStart
        end = self.bEnd
        for bBlock in bBlocks:
            bBlockStart = inPlaceBlocks.iat[bBlock,2]
            bBlockEnd = inPlaceBlocks.iat[bBlock,3]

            if bBlockStart - start < threshold and\
            end - bBlockEnd < threshold:
                start = end
                break
            elif bBlockStart < start and bBlockEnd < end:
                start = bBlockEnd
            else:
                blockBUni += bBlockStart - start
                if bBlockEnd< end:
                    start = bBlockEnd
                else:
                    start = end
                    break
        blockBUni += end - start

        self.overlappingInPlaceBlocks.extend([aBlocks, bBlocks])
        self.aUni = True if blockAUni > 1000 or blockAUni > 0.5*(self.aEnd-self.aStart) else False
        self.bUni = True if blockBUni > 1000 or blockBUni > 0.5*(self.bEnd-self.bStart) else False

    def checkOverlapWithSynBlocks_A(self,inPlaceBlocks, threshold):
        aBlocks = list(np.where((inPlaceBlocks.aStart.values < self.aEnd) & (inPlaceBlocks.aEnd.values > self.aStart) == True)[0])


        blockAUni = 0

        astart = self.aStart
        end = self.aEnd
        for aBlock in aBlocks:
            if inPlaceBlocks.iat[aBlock,0] - astart < threshold and\
                end - inPlaceBlocks.iat[aBlock,1] < threshold:
                astart = end
                break
            elif inPlaceBlocks.iat[aBlock,0] < astart and inPlaceBlocks.iat[aBlock, 1] < end:
                astart = inPlaceBlocks.iat[aBlock, 1]
            else:
                blockAUni += inPlaceBlocks.iat[aBlock,0] - astart
                if inPlaceBlocks.iat[aBlock,1] < end:
                    astart = inPlaceBlocks.iat[aBlock, 1]
                else:
                    astart = end
                    break
        blockAUni += end - astart

        self.overlappingInPlaceBlocks.append(aBlocks)
        self.aUni = True if blockAUni > 1000 or blockAUni > 0.5*(self.aEnd-self.aStart) else False

    def checkOverlapWithSynBlocks_B(self,inPlaceBlocks, threshold):
        bBlocks = list(np.where((inPlaceBlocks.bStart.values < self.bEnd) & (inPlaceBlocks.bEnd.values > self.bStart) == True)[0])

        blockBUni = 0
        start = self.bStart
        end = self.bEnd
        for bBlock in bBlocks:
            bBlockStart = inPlaceBlocks.iat[bBlock,2]
            bBlockEnd = inPlaceBlocks.iat[bBlock,3]

            if bBlockStart - start < threshold and\
            end - bBlockEnd < threshold:
                start = end
                break
            elif bBlockStart < start and bBlockEnd < end:
                start = bBlockEnd
            else:
                blockBUni += bBlockStart - start
                if bBlockEnd< end:
                    start = bBlockEnd
                else:
                    start = end
                    break
        blockBUni += end - start

        self.overlappingInPlaceBlocks.append(bBlocks)
        self.bUni = True if blockBUni > 1000 or blockBUni > 0.5*(self.bEnd-self.bStart) else False

    def addMEBlock(self, blockID):
        """List of Blocks which prohibit the entry of current block in the
        optimal solution"""
        self.meTo = np.array(blockID, dtype = "uint32")
        # try:
        #     self.meTo .extend(blockID) if type(blockID) == list else self.meTo.append(blockID)
        # except AttributeError:
        #     self.meTo = blockID if type(blockID) == list else [blockID]
    def setMEList(self, meAlist, meBlist):
        """Lists of a-overlap and b-overlap blocks. If at least 1 block has
        been selected from both lists then this block would become redundant"""
        self.meAlist = np.array(meAlist, dtype="uint32")
        self.meBlist = np.array(meBlist, dtype="uint32")

    def setStatus(self,stat):
        """stat = 1 ==> transBlock is important/necessary/unique"""
        self.status = stat



def getBestClusterSubset(cluster, transBlocksData, bRT):
    seedBlocks = [i for i in cluster if transBlocksData[i].status == 1]
    if len(cluster) < 50:
        output = bruteSubsetSelector(cluster, transBlocksData, seedBlocks, bRT)
        if output == "Failed":
            output = greedySubsetSelector(cluster, transBlocksData, seedBlocks)
    else:
        output = greedySubsetSelector(cluster, transBlocksData, seedBlocks)
    return output

def bruteSubsetSelector(cluster, transBlocksData, seedBlocks, bRT):
    posComb = [seedBlocks]
    skipList = [seedBlocks]
    for i in cluster:
        startTime = time.time()
        if hasattr(transBlocksData[i], "meTo"):
            newPosComb = []
            newSkipList = []
            for j in range(len(posComb)):
                if not any(a in posComb[j] for a in transBlocksData[i].meTo) and i not in skipList[j]:
                    newPosComb.append(posComb[j] + [i])
                    skipIndices = []
                    for k in posComb[j]:
                        if hasattr(transBlocksData[k],"meAlist"):
                            if i in transBlocksData[k].meAlist:
                                skipIndices.extend(transBlocksData[k].meBlist)
                            if i in transBlocksData[k].meBlist:
                                skipIndices.extend(transBlocksData[k].meAlist)
                        skipIndices.extend(transBlocksData[i].meTo)
                    newSkipList.append(skipList[j] + skipIndices)
            posComb.extend(newPosComb)
            skipList.extend(newSkipList)
        elif hasattr(transBlocksData[i], "meAlist"):
            newPosComb = []
            newSkipList = []
            for j in range(len(posComb)):
                check1 = not any(a in posComb[j] for a in transBlocksData[i].meAlist)
                check2 = not any(a in posComb[j] for a in transBlocksData[i].meBlist)
                if ( check1 or check2) and i not in skipList[j]:
                    newPosComb.append(posComb[j] + [i])
                    skipIndices = []
                    for k in posComb[j]:
                        if hasattr(transBlocksData[k],"meAlist"):
                            if i in transBlocksData[k].meAlist:
                                skipIndices.extend(transBlocksData[k].meBlist)
                            if i in transBlocksData[k].meBlist:
                                skipIndices.extend(transBlocksData[k].meAlist)
                        if k in transBlocksData[i].meAlist:
                            skipIndices.extend(transBlocksData[i].meBlist)
                        elif k in transBlocksData[i].meBlist:
                            skipIndices.extend(transBlocksData[i].meAlist)
                    for meElement in transBlocksData[i].meAlist:
                        if meElement in transBlocksData[i].meBlist:
                            skipIndices.append(meElement)
                    newSkipList.append(skipList[j] + skipIndices)
            posComb.extend(newPosComb)
            skipList.extend(newSkipList)
        else:
            newPosComb = []
            newSkipList = []
            for j in range(len(posComb)):
                if i not in skipList[j]:
                    newPosComb.append(posComb[j]+[i])
                    skipIndices = []
                    for k in posComb[j]:
                        if hasattr(transBlocksData[k],"meAlist"):
                            if i in transBlocksData[k].meAlist:
                                skipIndices.extend(transBlocksData[k].meBlist)
                            if i in transBlocksData[k].meBlist:
                                skipIndices.extend(transBlocksData[k].meAlist)
                    newSkipList.append(skipList[j]+skipIndices)
            posComb.extend(newPosComb)
            skipList.extend(newSkipList)
        timeTaken = time.time() - startTime
        remainingIterations = len(cluster) - cluster.index(i)

        if (timeTaken*(1.5**remainingIterations) > bRT):
            print("Cluster is too big for Brute Force\nTime taken for last iteration ",
                  timeTaken, " iterations remaining ",remainingIterations)
            return "Failed"


    if [] in posComb:
        posComb.remove([])
    ## Find the best set of alignments from a cluster applicable only in case where the number
    ## of all possible combinations is small
    bestScore = getScore(posComb[0], transBlocksData)
    bestComb = posComb[0]
    for i in range(1, len(posComb)):
        outBlocks = posComb[i]
        bestScore, bestComb = updateBestComb(bestScore, bestComb, outBlocks, transBlocksData)
    return( bestScore, bestComb)


def greedySubsetSelector(cluster, transBlocksData, seedBlocks, iterCount = 100):
    np.random.seed(1)
    bestScore = 0
    bestComb = []
    for i in range(iterCount):
        tempCluster = np.zeros(len(transBlocksData), dtype="uint8")
        outBlocks = np.zeros(len(transBlocksData), dtype="uint8")
        skipList = np.zeros(len(transBlocksData), dtype="uint8")
        tempCluster[cluster] = 1
        outBlocks[seedBlocks] = 1
        length = tempCluster.sum()
        tempCluster[seedBlocks] = 0
        transBlocksScore = {}
        for i in np.nonzero(tempCluster)[0]:
            transBlocksScore[i] = (transBlocksData[i].aEnd - transBlocksData[i].aStart) + (transBlocksData[i].bEnd - transBlocksData[i].bStart)
        while tempCluster.sum() > 0:
            while tempCluster.sum() != length:
                length = tempCluster.sum()
                for i in np.nonzero(tempCluster)[0]:
                    if hasattr(transBlocksData[i],"meTo"):
                        if outBlocks[transBlocksData[i].meTo].sum() > 0:
                            tempCluster[i] = 0
                            skipList[i]=1
                    elif hasattr(transBlocksData[i], "meAlist"):
                        if len(np.where(outBlocks[transBlocksData[i].meAlist] == 1)[0]) > 0 and len(np.where(outBlocks[transBlocksData[i].meBlist] == 1)[0]) > 0:
                            tempCluster[i] = 0
                            skipList[i] = 1
                for i in np.nonzero(tempCluster)[0]:
                    if hasattr(transBlocksData[i],"meTo"):
                        if skipList[transBlocksData[i].meTo].sum() == len(transBlocksData[i].meTo):
                            tempCluster[i] = 0
                            outBlocks[i]=1
                    elif hasattr(transBlocksData[i], "meAlist"):
                        if skipList[transBlocksData[i].meAlist].sum() == len(transBlocksData[i].meAlist) and skipList[transBlocksData[i].meBlist].sum() == len(transBlocksData[i].meBlist):
                            tempCluster[i] = 0
                            outBlocks[i] = 1

            if tempCluster.sum() > 0:
                topBlocks = sorted(np.nonzero(tempCluster)[0], key = lambda x: transBlocksScore[x], reverse = True)[:20]
                totalScore = sum(transBlocksScore[i] for i in topBlocks)
                prob = [transBlocksScore[i]/totalScore for i in topBlocks]
                newBlock = int(np.random.choice(topBlocks, size = 1, p = prob))
                outBlocks[newBlock] = 1
                tempCluster[newBlock] = 0
                if hasattr(transBlocksData[newBlock],"meTo"):
                    tempCluster[transBlocksData[newBlock].meTo] = 0
                    skipList[transBlocksData[newBlock].meTo] = 1
                elif hasattr(transBlocksData[newBlock],"meAlist"):
                    if outBlocks[transBlocksData[newBlock].meAlist].sum() > 0:
                        tempCluster[transBlocksData[newBlock].meBlist] = 0
                        skipList[transBlocksData[newBlock].meBlist] = 1
                    elif outBlocks[transBlocksData[newBlock].meAlist].sum() > 0:
                        tempCluster[transBlocksData[newBlock].meAlist] = 0
                        skipList[transBlocksData[newBlock].meBlist] = 1
                    for meElement in transBlocksData[newBlock].meAlist:
                        if meElement in transBlocksData[newBlock].meBlist:
                            tempCluster[meElement] = 0
                            skipList[meElement] = 1
        bestScore, bestComb = updateBestComb(bestScore, bestComb, np.nonzero(outBlocks)[0], transBlocksData)
    return(bestScore, bestComb)


def updateBestComb(bestScore, bestComb, outBlocks, transBlocksData):
    score = getScore(outBlocks, transBlocksData)
    if (score - bestScore > 1000) or (score > bestScore and len(outBlocks) <= len(bestComb)) or (bestScore - score < 1000 and len(outBlocks) < len(bestComb)):
        bestScore = score
        bestComb = outBlocks
    return(bestScore, bestComb)

def getScore(outBlocks, transBlocksData):
    aIndices = np.array([[transBlocksData[j].aStart, transBlocksData[j].aEnd] for j in outBlocks if transBlocksData[j].aUni])
    bIndices = np.array([[transBlocksData[j].bStart, transBlocksData[j].bEnd] for j in outBlocks if transBlocksData[j].bUni])
    aScore = count_uniq_elems(aIndices) if len(aIndices) > 0 else 0
    bScore = count_uniq_elems(bIndices) if len(bIndices) > 0 else 0
    return(aScore + bScore)


def count_uniq_elems(coordinates):
    a = coordinates[coordinates[:,0].argsort()]
    subs = a[1:,0] - a[:-1,1]                               ## overlapping regions
    overf = (a[:-1,1] - a[1:,1])
    return (a[:,1] - a[:,0]).sum() + subs[subs < 0].sum() + overf[overf > 0].sum()


def getTransClasses(clusterSolutionBlocks, transData, agroups, bgroups):
    def settl(j):
        if transData[j].dir == 1:
            transClasses["translocation"].append(j)
        elif transData[j].dir == -1:
            transClasses["invTranslocation"].append(j)
        else:
            print("ERROR ERROR ERROR", j)

    def setdup(j):
        if transData[j].dir == 1:
             transClasses["duplication"].append(j)
        elif transData[j].dir == -1:
            transClasses["invDuplication"].append(j)
        else:
            print("ERROR ERROR ERROR", j)

    transClasses = {"translocation":[],
                    "invTranslocation":[],
                    "duplication":[],
                    "invDuplication":[]}

    for i in clusterSolutionBlocks:
        for j in i:
            if not transData[j].aUni and not transData[j].bUni:
                print("ERROR ERROR ERROR", j)
            elif transData[j].status == 1:
                if not transData[j].aUni or not transData[j].bUni:
                    setdup(j)
                elif transData[j].aUni and transData[j].bUni:
                    if transData[j].genomeAUni and transData[j].genomeBUni:
                        settl(j)
                    elif not transData[j].genomeAUni:
                        istrans = 1
                        for k in transData[j].getoverlappingregions(groups = agroups,genome="a"):
                            if k in i:
                                if getScore([k], transData) >= getScore([j],transData):
                                    istrans = 0
                                    break
                        if istrans:
                            settl(j)
                        else:
                            setdup(j)
                    elif not transData[j].genomeBUni:
                        istrans = 1
                        for k in transData[j].getoverlappingregions(groups = bgroups,genome="b"):
                            if k in i:
                                if getScore([k], transData) >= getScore([j],transData):
                                    istrans = 0
                                    break
                        if istrans:
                            settl(j)
                        else:
                            setdup(j)
            elif not transData[j].aUni or not transData[j].bUni:
                setdup(j)
            elif transData[j].aUni and transData[j].bUni:
                if hasattr(transData[j],"meTo"):
                    if len(np.intersect1d(transData[j].meTo, i)) > 0:
                        setdup(j)
                    else:
                        settl(j)
                elif hasattr(transData[j],"meAlist"):
                    if len(np.intersect1d(transData[j].meAlist, i)) > 0 or\
                    len(np.intersect1d(transData[j].meBlist, i)) > 0:
                        setdup(j)
                    else:
                        settl(j)
                else:
                     print("ERROR ERROR ERROR", j)
    return transClasses


def getCTX(coords, cwdPath, uniChromo, threshold, bRT, prefix, tUC, tUP, nCores):
    logger = logging.getLogger("getCTX")
    logger.info("Identifying cross-chromosomal translocation and duplication for chromosome" + str(datetime.now()))

    def getDupCTX(indices, allTransBlocksData, transClasses):
        dupGenomes = {}
        for index in indices:
            if index in transClasses["translocation"] or index in transClasses["invTranslocation"]:
                dupGenomes[index] = ""
                continue
            found = False
            tempTransBlock = allTransBlocksData[index]
            if not tempTransBlock.aUni:
                dupGenomes[index] = "B"
                continue
            elif not tempTransBlock.bUni:
                dupGenomes[index] = "A"
                continue
            elif tempTransBlock.genomeAUni:
                dupGenomes[index] = "A"
                continue
            elif tempTransBlock.genomeBUni:
                dupGenomes[index] = "B"
                continue
            for i in tempTransBlock.meAlist:
                if i in transClasses["translocation"] or i in transClasses["invTranslocation"]:
                    found = True
                    dupGenomes[index] = "B"
                    break
            if not found:
                dupGenomes[index] = "A"
        return(dupGenomes)

    def printCTX(cwdPath, clusterSolutionBlocks, ctxBlocksData, orderedBlocks, invertedBlocks ,transBlocks, invTransBlocks, ctxTransIndexOrder, ctxTransBlocks, genomeagroups, genomebgroups):
        transClasses = getTransClasses(clusterSolutionBlocks, ctxBlocksData, genomeagroups, genomebgroups)
        indices = sorted(unlist(list(transClasses.values())))
        keys = [key for index in indices for key in list(transClasses.keys()) if index in transClasses[key]]
        blocksClasses = dict(zip(indices,keys))
        dupGenomes = getDupCTX(indices, ctxBlocksData, transClasses)

        fout = open(cwdPath+prefix+"ctxOut.txt","w")
        for index in indices:
            if ctxBlocksData[index].dir == 1:
                alignIndices = transBlocks[ctxTransIndexOrder[index]]
                fout.write("#\t" + "\t".join(map(str,[ctxTransBlocks.iloc[index]["aChr"], ctxTransBlocks.iloc[index]["aStart"], ctxTransBlocks.iloc[index]["aEnd"], "-", ctxTransBlocks.iloc[index]["bChr"], ctxTransBlocks.iloc[index]["bStart"],ctxTransBlocks.iloc[index]["bEnd"]])) + "\t" + blocksClasses[index]+ "\t" +  dupGenomes[index]+"\n")
                for i in alignIndices:
                    fout.write("\t".join(map(str,orderedBlocks.iloc[i,0:4]))+"\n")
            elif ctxBlocksData[index].dir == -1:
                alignIndices = invTransBlocks[ctxTransIndexOrder[index]]
                fout.write("#\t" + "\t".join(map(str,[ctxTransBlocks.iloc[index]["aChr"], ctxTransBlocks.iloc[index]["aStart"], ctxTransBlocks.iloc[index]["aEnd"], "-", ctxTransBlocks.iloc[index]["bChr"], ctxTransBlocks.iloc[index]["bStart"],ctxTransBlocks.iloc[index]["bEnd"]]))+ "\t" + blocksClasses[index] + "\t" +  dupGenomes[index]+ "\n")
                for i in alignIndices:
                    fout.write("\t".join(map(str,invertedBlocks.iloc[i,[0,1,3,2]])) + "\n")
        fout.close()


    logger.debug("Reading Coords" + str(datetime.now()))

    annoCoords = readAnnoCoords(cwdPath, uniChromo, prefix)
    ctxData = coords.loc[coords['aChr'] != coords['bChr']].copy()
    ctxData.index = range(len(ctxData))
    invCTXIndex = ctxData.index[ctxData.bDir == -1]
    ctxData.loc[invCTXIndex,"bStart"] = ctxData.loc[invCTXIndex].bStart + ctxData.loc[invCTXIndex].bEnd
    ctxData.loc[invCTXIndex, "bEnd"] = ctxData.loc[invCTXIndex].bStart - ctxData.loc[invCTXIndex].bEnd
    ctxData.loc[invCTXIndex, "bStart"] = ctxData.loc[invCTXIndex].bStart - ctxData.loc[invCTXIndex].bEnd
    ctxData.sort_values(by= ["aChr","aStart","aEnd","bChr","bStart","bEnd"], inplace = True)
    ctxData["aIndex"] = range(ctxData.shape[0])
    ctxData.sort_values(by= ["bChr","bStart","bEnd","aChr","aStart","aEnd"], inplace = True)
    ctxData["bIndex"] = range(ctxData.shape[0])
    ctxData.sort_values("aIndex", inplace = True)

    logger.debug("CTX identification: ctxdata size" + str(ctxData.shape))

    orderedBlocks = ctxData[ctxData.bDir == 1]
    invertedBlocks = ctxData[ctxData.bDir == -1]

    ## Create connectivity tree for directed blocks

    logger.debug("Making Tree")

    nCores = nCores if nCores < 2 else 2

    with Pool(processes = nCores) as pool:
        blks = pool.map(partial(getBlocks, annoCoords=annoCoords, threshold=threshold, tUC=tUC, tUP=tUP), [orderedBlocks,invertedBlocks])
    transBlocks = blks[0]
    invTransBlocks = blks[1]
    del(blks)
    collect()
    logger.debug("finding Blocks")
    logger.debug("Preparing for cluster analysis")

    ctxTransBlocks, ctxTransIndexOrder = mergeTransBlocks(transBlocks, orderedBlocks, invTransBlocks, invertedBlocks, ctx = True)

    ctxTransGenomeAGroups = []
    for chromo in uniChromo:
        ctxTransGenomeAGroups += makeTransGroupList(ctxTransBlocks.loc[ctxTransBlocks.aChr == chromo, ["aStart","aEnd","bStart","bEnd"]], "aStart","aEnd",threshold)


    ctxTransGenomeBGroups = []
    for chromo in uniChromo:
        ctxTransGenomeBGroups += makeTransGroupList(ctxTransBlocks.loc[ctxTransBlocks.bChr == chromo, ["aStart","aEnd","bStart","bEnd"]], "bStart","bEnd",threshold)


    ctxGroupIndices = {}
    for i in range(len(ctxTransGenomeAGroups)):
        for block in ctxTransGenomeAGroups[i].member:
            ctxGroupIndices[block] = [i]

    for i in range(len(ctxTransGenomeBGroups)):
        for block in ctxTransGenomeBGroups[i].member:
            ctxGroupIndices[block].append(i)

    ctxCluster = getTransCluster(ctxGroupIndices, ctxTransGenomeAGroups, ctxTransGenomeBGroups)

    ctxClusterIndices = dict()
    for i in range(len(ctxCluster)):
        ctxClusterIndices.update(dict.fromkeys(ctxCluster[i], i))

    if len(ctxTransBlocks) > 0:
        auni = getOverlapWithSynBlocks(np.array(ctxTransBlocks.aStart), np.array(ctxTransBlocks.aEnd),
                                     np.array(ctxTransBlocks.aChr), np.array(annoCoords.aStart),
                                     np.array(annoCoords.aEnd), np.array(annoCoords.aChr), 50,
                                     ctxTransBlocks.shape[0], tUC, tUP)

        sortedInPlace = annoCoords.sort_values(["bStart", "bEnd"])
        buni = getOverlapWithSynBlocks(np.array(ctxTransBlocks.bStart), np.array(ctxTransBlocks.bEnd),
                                   np.array(ctxTransBlocks.bChr), np.array(sortedInPlace.bStart),
                                   np.array(sortedInPlace.bEnd), np.array(sortedInPlace.bChr), 50,
                                     ctxTransBlocks.shape[0], tUC, tUP)

    genomeGroupLengths = ([len(i.member) for i in ctxTransGenomeAGroups], [ len(i.member) for i in ctxTransGenomeBGroups])

    ctxBlocksData = [transBlock(i) for i in range(ctxTransBlocks.shape[0])]

    count=0
    for row in ctxTransBlocks.itertuples(index = False):
        ctxBlocksData[count].aStart = row.aStart
        ctxBlocksData[count].aEnd = row.aEnd
        ctxBlocksData[count].bStart = row.bStart
        ctxBlocksData[count].bEnd = row.bEnd
        ctxBlocksData[count].dir = row.bDir
        ctxBlocksData[count].transClusterIndex = ctxClusterIndices[count]
        ctxBlocksData[count].transGroupIndices = ctxGroupIndices[count]
        ctxBlocksData[count].aUni = auni[count]
        ctxBlocksData[count].bUni = buni[count]
        if genomeGroupLengths[0][ctxBlocksData[count].transGroupIndices[0]] == 1:
            ctxBlocksData[count].genomeAUni = True
        if genomeGroupLengths[1][ctxBlocksData[count].transGroupIndices[1]] == 1:
            ctxBlocksData[count].genomeBUni = True
        if (ctxBlocksData[count].aUni and ctxBlocksData[count].genomeAUni) or (ctxBlocksData[count].bUni and ctxBlocksData[count].genomeBUni):
            ctxBlocksData[count].setStatus(1)
        count+=1

    aUni = np.array([ctxBlocksData[i].aUni for i in range(ctxTransBlocks.shape[0])], dtype="int")
    bUni = np.array([ctxBlocksData[i].bUni for i in range(ctxTransBlocks.shape[0])], dtype="int")
    status = np.array([ctxBlocksData[i].status for i in range(ctxTransBlocks.shape[0])], dtype="int")
    aIndex = np.array([ctxBlocksData[i].transGroupIndices[0] for i in range(ctxTransBlocks.shape[0])], dtype="int")
    bIndex = np.array([ctxBlocksData[i].transGroupIndices[1] for i in range(ctxTransBlocks.shape[0])], dtype="int")
    aGroups = {i: np.array(ctxTransGenomeAGroups[i].member, dtype="int") for i in range(len(ctxTransGenomeAGroups))}
    bGroups = {i: np.array(ctxTransGenomeBGroups[i].member, dtype="int") for i in range(len(ctxTransGenomeBGroups))}

    if len(ctxTransBlocks) > 0:
        out = getmeblocks(np.array(ctxTransBlocks.aStart), np.array(ctxTransBlocks.aEnd), np.array(ctxTransBlocks.bStart),
                          np.array(ctxTransBlocks.bEnd), 50, ctxTransBlocks.shape[0], aUni, bUni, status, aIndex, bIndex,
                          aGroups, bGroups)

        for i in range(len(out[0])):
            if out[0][i]:
                ctxCluster[ctxClusterIndices[i]].remove(i)

        for i in out[1].keys():
            if len(out[1][i]) > 0:
                ctxBlocksData[i].addMEBlock(list(out[1][i]))

        for i in out[2].keys():
            ctxBlocksData[i].setMEList(list(out[2][i][0]),list(out[2][i][1]))

        del(aUni, bUni, status, aIndex, bIndex, aGroups, bGroups, out)
        collect()

    logger.debug("Finding clusters")

    clusterSolutions = []
    for i in range(len(ctxCluster)):
        tempCluster = ctxCluster[i].copy()
        if len(tempCluster) == 0:
            continue
        else:
            clusterSolutions.append(getBestClusterSubset(tempCluster, ctxBlocksData, bRT))

    clusterSolutionBlocks = [i[1] for i in clusterSolutions]

    printCTX(cwdPath, clusterSolutionBlocks, ctxBlocksData, orderedBlocks, invertedBlocks ,transBlocks, invTransBlocks, ctxTransIndexOrder, ctxTransBlocks, ctxTransGenomeAGroups, ctxTransGenomeBGroups)
    return 0
