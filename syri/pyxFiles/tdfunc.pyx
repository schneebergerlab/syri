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

    def addMEBlock(self, blockID):
        """List of Blocks which prohibit the entry of current block in the
        optimal solution"""
        self.meTo = np.array(blockID, dtype = "uint32")

    def setMEList(self, meAlist, meBlist):
        """Lists of a-overlap and b-overlap blocks. If at least 1 block has
        been selected from both lists then this block would become redundant"""
        self.meAlist = np.array(meAlist, dtype="uint32")
        self.meBlist = np.array(meBlist, dtype="uint32")

    def setStatus(self,stat):
        """stat = 1 ==> transBlock is important/necessary/unique"""
        self.status = stat


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


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline getmeblocks2(unsigned long[:] ast, unsigned long[:] aen, unsigned long[:] bst, unsigned long[:] ben, cpp_vec[long] amem, cpp_vec[long] bmem, int i, int threshold, int meclass):
    cdef:
        unsigned short int[:]                   meb, meb_a, meb_b
        Py_ssize_t                              index
        long                                    j
    if meclass == 1:
        meb = np.zeros(bmem.size(), dtype=np.uint16)              ## vector of mutually exclusive block
        for index in range(<Py_ssize_t>bmem.size()):
            j = bmem[index]
            if ben[j] < bst[i]:
                continue
            if bst[j] > ben[i]:
                break
            if bst[j] - threshold < bst[i] and ben[j] + threshold > ben[i]:
                if j!=i:
                    meb[index]=1
        return np.array([bmem[index] for index in range(<Py_ssize_t>bmem.size()) if meb[index] == 1], np.uint)
    elif meclass == 2:
        meb = np.zeros(amem.size(), dtype=np.uint16)               ## vector of mutually exclusive block
        for index in range(<Py_ssize_t>amem.size()):
            j = amem[index]
            if aen[j] < ast[i]:
                continue
            if ast[j] > aen[i]:
                break
            if ast[j] - threshold < ast[i] and aen[j]+threshold > aen[i]:
                if j!=i:
                    meb[index] = 1
        return np.array([amem[index] for index in range(<Py_ssize_t>amem.size()) if meb[index] ==1], np.uint)
    elif meclass == 3:
        meb_a = np.zeros(amem.size(), dtype=np.uint16)             ## vector of mutually exclusive block on A genome
        for index in range(<Py_ssize_t>amem.size()):
            j = amem[index]
            if aen[j] < ast[i]:
                continue
            if ast[j] > aen[i]:
                break
            if ast[j] - threshold < ast[i] and aen[j]+threshold > aen[i]:
                if j!=i:
                    meb_a[index] = 1

        meb_b = np.zeros(bmem.size(), dtype=np.uint16)             ## vector of mutually exclusive block on B genome
        for index in range(<Py_ssize_t> bmem.size()):
            j = bmem[index]
            if ben[j] < bst[i]:
                continue
            if bst[j] > ben[i]:
                break
            if bst[j] - threshold < bst[i] and ben[j] + threshold > ben[i]:
                if j!=i:
                    meb_b[index] = 1
        return np.array([amem[index] for index in range(<Py_ssize_t>amem.size()) if meb_a[index]==1], np.uint), np.array([bmem[index] for index in range(<Py_ssize_t>bmem.size()) if meb_b[index] == 1], np.uint)



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


def count_uniq_elems(coordinates):
    a = coordinates[coordinates[:,0].argsort()]
    subs = a[1:,0] - a[:-1,1]                               ## overlapping regions
    overf = (a[:-1,1] - a[1:,1])
    return (a[:,1] - a[:,0]).sum() + subs[subs < 0].sum() + overf[overf > 0].sum()


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


def bruteSubsetSelector(cluster, transBlocksData, seedBlocks, bRT):
    logger = logging.getLogger('Brute-force TD identification')
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
            logger.info("Cluster is too big for Brute Force\nTime taken for last iteration " +
                  str(timeTaken) + ". iterations remaining " + str(remainingIterations))
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


def getCTX(coords, cwdPath, uniChromo, threshold, bRT, prefix, tUC, tUP, nCores, tdgl):
    logger = logging.getLogger("getCTX")
    logger.info("Identifying cross-chromosomal translocation and duplication for chromosome" + str(datetime.now()))

    def getDupCTX(indices, allTransBlocksData, transClasses, astart, aend, bstart, bend, aindex, bindex, agroup, bgroup, threshold, meclass):
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
            if meclass[index] != 3:
                print('Error in dup class identification ', index)
                continue
            mealist, meblist = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[index]],  bgroup[bindex[index]], index, threshold, meclass=3)
            for i in mealist:
                if i in transClasses["translocation"] or i in transClasses["invTranslocation"]:
                    found = True
                    dupGenomes[index] = "B"
                    break
            if not found:
                dupGenomes[index] = "A"
        return(dupGenomes)

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

    nCorestemp = nCores if nCores < 2 else 2

    with Pool(processes = nCorestemp) as pool:
        blks = pool.starmap(partial(getBlocks, annoCoords=annoCoords, threshold=threshold, tUC=tUC, tUP=tUP, tdgl=tdgl), [[orderedBlocks, 0], [invertedBlocks, 1]])

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

    logger.debug("Getting clusters")
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

    logger.debug("Finding ME Blocks")
    aUni = np.array([ctxBlocksData[i].aUni for i in range(ctxTransBlocks.shape[0])], dtype="int")
    bUni = np.array([ctxBlocksData[i].bUni for i in range(ctxTransBlocks.shape[0])], dtype="int")
    status = np.array([ctxBlocksData[i].status for i in range(ctxTransBlocks.shape[0])], dtype="int")
    aIndex = np.array([ctxBlocksData[i].transGroupIndices[0] for i in range(ctxTransBlocks.shape[0])], dtype="int")
    bIndex = np.array([ctxBlocksData[i].transGroupIndices[1] for i in range(ctxTransBlocks.shape[0])], dtype="int")

    ## get sorted values. sorted based on genome coordinate in ctxTransBlocks
    aGroups = {}
    for i in range(len(ctxTransGenomeAGroups)):
        aGroups[i] = ctxTransBlocks.iloc[ctxTransGenomeAGroups[i].member].sort_values(['aStart', 'aEnd']).index.values
    bGroups = {}
    for i in range(len(ctxTransGenomeBGroups)):
        bGroups[i] = ctxTransBlocks.iloc[ctxTransGenomeBGroups[i].member].sort_values(['bStart', 'bEnd']).index.values
    clstrsize = np.array([len(ctxCluster[i.transClusterIndex]) for i in ctxBlocksData], np.int)

    if len(ctxTransBlocks) > 0:
        out = getmeblocks(np.array(ctxTransBlocks.aStart),
                          np.array(ctxTransBlocks.aEnd),
                          np.array(ctxTransBlocks.bStart),
                          np.array(ctxTransBlocks.bEnd),
                          threshold,
                          aUni,
                          bUni,
                          status,
                          aIndex,
                          bIndex,
                          aGroups,
                          bGroups,
                          clstrsize)

        for i in range(len(out[0])):
            if out[0][i]:
                ctxCluster[ctxClusterIndices[i]].remove(i)

        for i in out[1].keys():
            if len(out[1][i]) > 0:
                ctxBlocksData[i].addMEBlock(list(out[1][i]))

        for i in out[2].keys():
            ctxBlocksData[i].setMEList(list(out[2][i][0]),list(out[2][i][1]))
        #
        # del(aUni, bUni, status, aIndex, bIndex, aGroups, bGroups, out)
        # collect()

    logger.debug("Finding best subset of clusters")

    with Pool(processes=nCores) as pool:
        clusterSolutions = pool.map(partial(getBestClusterSubset,
                                             transBlocksData=ctxBlocksData,
                                             chromo='CTX',
                                             bRT=bRT,
                                             aGroups=aGroups,
                                             bGroups=bGroups,
                                             threshold=threshold), ctxCluster)

    clusterSolutionBlocks = [i[1] for i in clusterSolutions if i != None]

    garb = deque()
    for i in range(len(ctxBlocksData)):
        if not ctxBlocksData[i].aUni and not ctxBlocksData[i].bUni:
            garb.append(0)
        elif ctxBlocksData[i].status == 1:
            garb.append(0)
        elif not ctxBlocksData[i].aUni:
            garb.append(1)
        elif not ctxBlocksData[i].bUni:
            garb.append(2)
        else:
            garb.append(3)
    meclass = np.array(list(garb), np.uint16)

    transClasses = getTransClasses(clusterSolutionBlocks,
                                   ctxBlocksData,
                                   ctxTransGenomeAGroups,
                                   ctxTransGenomeBGroups,
                                   ctxTransBlocks.aStart.values.astype(np.uint),
                                   ctxTransBlocks.aEnd.values.astype(np.uint),
                                   ctxTransBlocks.bStart.values.astype(np.uint),
                                   ctxTransBlocks.bEnd.values.astype(np.uint),
                                   aIndex,
                                   bIndex,
                                   aGroups,
                                   bGroups,
                                   threshold,
                                   meclass)

    indices = sorted(unlist(list(transClasses.values())))
    keys = [key for index in indices for key in list(transClasses.keys()) if index in transClasses[key]]
    blocksClasses = dict(zip(indices,keys))
    dupGenomes = getDupCTX(indices,
                           ctxBlocksData,
                           transClasses,
                           ctxTransBlocks.aStart.values.astype(np.uint),
                           ctxTransBlocks.aEnd.values.astype(np.uint),
                           ctxTransBlocks.bStart.values.astype(np.uint),
                           ctxTransBlocks.bEnd.values.astype(np.uint),
                           aIndex,
                           bIndex,
                           aGroups,
                           bGroups,
                           threshold,
                           meclass)

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

    return 0

cpdef makeBlocksTree(long[:] aStart, long[:] aEnd, long[:] bStart, long[:] bEnd, int threshold, long[:] left, long[:] right, int tdgl):
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

    for i in range(<Py_ssize_t> n):
        for j in range(<Py_ssize_t> i, <Py_ssize_t> n):
            if (aStart[j] - aEnd[i]) < tdgl:        # Select only alignments with small gaps
                if (aStart[j] - aStart[i]) > threshold:
                    if (aEnd[j] - aEnd[i]) > threshold:
                        if (bStart[j] - bEnd[i]) < tdgl:        # Select only alignments with small gaps
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

cpdef makeBlocksTree_ctx(long[:] astart, long[:] aend, long[:] bstart, long[:] bend, long[:] bdir, np.ndarray achr, np.ndarray bchr, int threshold, int tdgl):
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
    for i in range(<Py_ssize_t> n):
        for j in range(<Py_ssize_t> i, <Py_ssize_t> n):
            if achr_int[i] != achr_int[j]:
                continue
            if bchr_int[i] != bchr_int[j]:
                continue
            if (astart[j] - aend[i]) < tdgl:
                if (astart[j] - astart[i]) > threshold:
                    if (aend[j] - aend[i]) > threshold:
                        if bdir[i] == 1:
                            if (bstart[j] - bend[i]) < tdgl:
                                if (bstart[j] - bstart[i]) > threshold:
                                    if (bend[j] - bend[i]) > threshold:
                                        out[i].push_back(j)
                        else:
                            if (bend[i] - bstart[j]) < tdgl:
                                if (bstart[i] - bstart[j]) > threshold:
                                    if (bend[i] - bend[j]) > threshold:
                                        out[i].push_back(j)
    it = out.begin()
    outdict = {}
    while it != out.end():
        outdict[deref(it).first] = [deref(it).second[i] for i in range(<Py_ssize_t> deref(it).second.size())]
        inc(it)
    return outdict

# cpdef makeBlocksTree(long[:] aStart, long[:] aEnd, long[:] bStart, long[:] bEnd, int threshold, long[:] left, long[:] right):
#     """Compute whether two alignments can be part of one translation block. For this:
#         the alignments should not be separated by any inPlaceBlock on both ends and
#         they should be collinear with respect to each other.
#
#        Returns
#        --------
#        outOrderedBlocks: pandas DataFrame,
#            Dataframe of type Object. Lower half is NA, upper half contains whether two
#            alignments can be connected (True) or not (False).
#     """
#
#     cdef:
#         ssize_t                                     i,j, n = len(aStart)
#         cpp_map[long, cpp_deq[long]]                out
#         cpp_map[long, cpp_deq[long]].iterator       it
#
#     for i in range(<Py_ssize_t> n):
#         for j in range(<Py_ssize_t> i, <Py_ssize_t> n):
#             if (aStart[j] - aStart[i]) > threshold:
#                 if (aEnd[j] - aEnd[i]) > threshold:
#                     if (bStart[j] - bStart[i]) > threshold:
#                         if (bEnd[j] - bEnd[i]) > threshold:
#                             # Alignments could form a block when they are not separated by inPlaceBlocks. For this, we check whether the alignments have a common intersecting inPlaceAlignments.
#                             if right[i]-1 >= left[j]+1 and right[i] <= right[j]:
#                                 out[i].push_back(j)
#                             elif left[i] >= left[j] and left[i]+1 <= right[j]-1:
#                                 out[i].push_back(j)
#
#     it = out.begin()
#     outdict = {}
#     while it != out.end():
#         outdict[deref(it).first] = [deref(it).second[i] for i in range(<Py_ssize_t> deref(it).second.size())]
#         inc(it)
#     return(outdict)
#
#
# cpdef makeBlocksTree_ctx(long[:] astart, long[:] aend, long[:] bstart, long[:] bend, long[:] bdir, np.ndarray achr, np.ndarray bchr, int threshold):
#     """Compute whether two alignments can be part of one translation block. For this:
#        they should be syntenic with respect to each other.
#     """
#     cdef:
#         Py_ssize_t                                  i,j, n = len(astart)
#         cpp_map[long, cpp_deq[long]]                out
#         cpp_map[long, cpp_deq[long]].iterator       it
#         long[:]                                     achr_int, bchr_int
#     _, achr_int = np.unique(achr, return_inverse=True)
#     _, bchr_int = np.unique(bchr, return_inverse=True)
#     for i in range(n):
#         for j in range(i,n):
#             if achr_int[i] != achr_int[j]:
#                 continue
#             if bchr_int[i] != bchr_int[j]:
#                 continue
#             if (astart[j] - astart[i]) > threshold:
#                 if (aend[j] - aend[i]) > threshold:
#                     if bdir[i] == 1:
#                         if (bstart[j] - bstart[i]) > threshold:
#                             if (bend[j] - bend[i]) > threshold:
#                                 out[i].push_back(j)
#                     else:
#                         if (bstart[i] - bstart[j]) > threshold:
#                             if (bend[i] - bend[j]) > threshold:
#                                 out[i].push_back(j)
#
#     it = out.begin()
#     outdict = {}
#     while it != out.end():
#         outdict[deref(it).first] = [deref(it).second[i] for i in range(<Py_ssize_t> deref(it).second.size())]
#         inc(it)
#     return outdict
#

def getBlocks(orderedBlocks, isinv, annoCoords, threshold, tUC, tUP, tdgl):
    if len(orderedBlocks) == 0:
        return([])

    outOrderedBlocks = makeBlocksTree_ctx(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, orderedBlocks.bDir.values, orderedBlocks.aChr.values, orderedBlocks.bChr.values, threshold, tdgl)

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
                                     isinv)
    return(transBlocks)



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
        # print(i, str(datetime.now()))
        # st = datetime.now()
        #
        count+=1
        if count % 1000 == 0:
            print(count, str(datetime.now()))
        if count == brk:
            break

        # if i%100 == 0:
        #     print(i, str(datetime.now()))
        # print(i, 'I1', str(datetime.now() - st))
        nodepath.clear()
        pred = pred_bkp.copy()
        dist = dist_bkp.copy()
        dist[i] = 0

        # Process vertices in topological order

        # print(i, 'I2', str(datetime.now() - st))
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
        # print(i, 'A', str(datetime.now() - st))

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
                # print(i, 'D', str(datetime.now() - st))
    print('End', datetime.now())
    return out

# cpdef getProfitableTrans(cpp_map[long, cpp_set[long]] graph, long[:] astart, long[:] aend, long[:] bstart, long[:] bend, np.ndarray achr, np.ndarray bchr, float[:] iden, long[:] alen, long[:] blen, long[:] inastart, long[:] inaend, long[:] inbstart, long[:] inbend, np.ndarray inachr, np.ndarray inbchr, long tUC, float tUP, int isinv = 0):
#     """
#     Input:
#      1) dictionary in which each key corresponds to trans alignment and the corresponding values are the alignments which are colinear with key.
#      2) Coordinates of the trans alignments
#      3) Coordinates of inplace blocks. Sorted separately for reference and query genome.
#      4) tUC, tUP
#
#     Output:
#      1) All trans blocks (groups of alignments) with high positive score.
#     """
#     cdef:
#         Py_ssize_t                                  i, j, k
#         long                                        id, current
#         unsigned long                               cnt
#         long                                        nodecnt, edgecnt, len_in = len(inastart)
#         long                                        ast, aen, bst, ben
#         float                                       ascore, bscore, agap, bgap
#         long                                        al, bl, au, bu
#         float                                       score
#         float[:]                                    weight
#         float[:]                                    dist
#         long[:]                                     pred
#         long[:]                                     topo, indegree, source, target
#         long[:]                                     n = np.array(range(len(astart)), dtype='int')
#         long[:]                                     achrint, bchrint, inachrint, inbchrint
#         cpp_set[long]                               rem
#         cpp_que[long]                               q, toporder
#         cpp_deq[long]                               path, r_path
#         cpp_map[long, cpp_deq[long]]                nodepath
#         cpp_map[long, cpp_vec[long]]                almntdata
#         cpp_set[long].iterator                      set_it, set_it2
#         cpp_deq[long].reverse_iterator              deq_rit
#         cpp_map[long, cpp_set[long]].iterator       mapit
#
#
#     nodecnt = len(n)
#     topo = indegree = np.zeros(nodecnt, dtype='int64')
#
#     ## Check that the keys are sorted in the graph
#     mapit = graph.begin()
#     id = -1
#     while mapit != graph.end():
#         if id >= deref(mapit).first:
#             print('ERROR: unsorted outOrderedBlocks')
#         else:
#             id = deref(mapit).first
#         inc(mapit)
#
#
#     ## Remove children for which another child is present between parent and itself
#     mapit = graph.begin()
#     while mapit != graph.end():
#         rem.clear()                     ## List of children to be removed
#         set_it = deref(mapit).second.begin()
#         while set_it != deref(mapit).second.end():
#             if rem.count(deref(set_it)) == 0:
#                 if graph.count(deref(set_it)) == 1:
#                     set_it2 = graph[deref(set_it)].begin()
#                     while set_it2 !=  graph[deref(set_it)].end():
#                         rem.insert(deref(set_it2))
#                         inc(set_it2)
#             inc(set_it)
#         set_it = rem.begin()
#         while set_it != rem.end():
#             deref(mapit).second.erase(deref(set_it))
#             inc(set_it)
#         inc(mapit)
#
#
#     # Get number of edges in the graph
#     edgecnt = 0
#     mapit = graph.begin()
#     while mapit != graph.end():
#         edgecnt += <Py_ssize_t> deref(mapit).second.size()
#         inc(mapit)
#
#
#     # Find the topological order in which the nodes should be processed to find paths
#
#     ## Get indegree for all nodes (number of nodes == number of aligments)
#     mapit = graph.begin()
#     while mapit != graph.end():
#         set_it = deref(mapit).second.begin()
#         while set_it != deref(mapit).second.end():
#             indegree[deref(set_it)]+=1
#             inc(set_it)
#         inc(mapit)
#
#
#     ## Push all nodes with indegree=0 in queue
#     for i in n:
#         if indegree[i] == 0:
#             q.push(i)
#     cnt = 0
#
#     ## Get topological ordering for the nodes
#     while q.size() > 0:
#         id = q.front()
#         q.pop()
#         toporder.push(id)
#         if graph.count(id) == 1:
#             set_it = graph[id].begin()
#             while set_it != graph[id].end():
#                 indegree[deref(set_it)]-=1
#                 if indegree[deref(set_it)]==0:
#                     q.push(deref(set_it))
#                 inc(set_it)
#         cnt += 1
#     if cnt != len(indegree):
#         print('ERROR: Cycle found')
#     if toporder.size() != len(topo):
#         print('ERROR: topological ordering didnt cover all nodes')
#     for i in n:
#         topo[i] = toporder.front()
#         toporder.pop()
#
#
#     # Get order in which the edges need to be transversed
#     source = np.zeros(edgecnt, dtype='int')
#     target = np.zeros(edgecnt, dtype='int')
#     weight = np.zeros(edgecnt, dtype='float32')
#     id = 0
#     for i in topo:
#         if graph.count(i) == 1:
#             set_it = graph[i].begin()
#             while set_it != graph[i].end():
#                 source[id] = i
#                 target[id] = deref(set_it)
#                 weight[id] = (aend[i] + bend[i] - astart[i] - bstart[i] + 2) * iden[i]
#                 id+=1
#                 inc(set_it)
#
#
#     ## Convert arary of 'string' chromosome ids to much faster numeric chromosome ids
#     if list(np.unique(achr)) == list(np.unique(bchr)) == list(np.unique(inachr)) == list(np.unique(inbchr)):
#         _, achrint      = np.unique(achr, return_inverse = True)
#         _, bchrint      = np.unique(bchr, return_inverse = True)
#         _, inachrint    = np.unique(inachr, return_inverse = True)
#         _, inbchrint    = np.unique(inbchr, return_inverse = True)
#     else:
#         unchrs = np.unique(list(np.unique(achr)) + list(np.unique(bchr)) + list(np.unique(inachr)) + list(np.unique(inbchr)))
#         unchrdict = {unchrs[i]:i for i in range(len(unchrs))}
#         achrint     = np.array([unchrdict[c] for c in achr], np.int)
#         bchrint     = np.array([unchrdict[c] for c in bchr], np.int)
#         inachrint   = np.array([unchrdict[c] for c in inachr], np.int)
#         inbchrint   = np.array([unchrdict[c] for c in inbchr], np.int)
#
#
#     ## For each alignment/node calculate the number of bases which are not overlapping with the in-place blocks
#     for i in n:
#         ascore = 0
#         bscore = 0
#         ast = astart[i]
#         aen = aend[i]
#         bst = bstart[i]
#         ben = bend[i]
#
#         # print(ast, aen, bst, ben)
#         for j in range(len_in):
#             if achrint[i] == inachrint[j]:
#                 if inastart[j] > aen:
#                     break
#                 if inaend[j] < ast:
#                     continue
#                 if inastart[j] < ast:
#                     if inaend[j] < aen:
#                         ast = inaend[j]+1
#                     else:
#                         ast = aen+1
#                         break
#                 else:
#                     ascore += inastart[j] - ast
#                     if inaend[j] < aen:
#                         ast = inaend[j]+1
#                     else:
#                         ast = aen+1
#                         break
#         ascore += aen - ast + 1
#
#         for j in range(len_in):
#             if bchrint[i] == inbchrint[j]:
#                 if inbstart[j] > ben:
#                     break
#                 if inbend[j] < bst:
#                     continue
#                 if inbstart[j] < bst:
#                     if inbend[j] < ben:
#                         bst = inbend[j]+1
#                     else:
#                         bst = ben+1
#                         break
#                 else:
#                     bscore += inbstart[j] - bst
#                     if inbend[j] < ben:
#                         bst = inbend[j]+1
#                     else:
#                         bst = ben+1
#                         break
#         bscore += ben - bst + 1
#         almntdata[i] = np.array([aend[i]-astart[i]+1, bend[i]-bstart[i]+1, ascore, bscore], dtype=np.int)
#
#     out = deque()
#     for i in n:
#         # if i%100 == 0:
#         #     print(i, str(datetime.now()))
#         nodepath.clear()
#         pred = np.array([-1]*<Py_ssize_t>len(n), dtype = np.int)
#         dist = np.array([np.float32('inf')]*<Py_ssize_t>len(n), dtype = np.float32)
#         dist[i] = 0
#
#         # Process vertices in topological order
#         id = 0
#         for j in n:
#             for k in range(id, edgecnt):
#                 if source[k] != topo[j]:
#                     break
#                 if dist[target[k]] > dist[source[k]] + weight[k]:
#                     dist[target[k]] = dist[source[k]] + weight[k]
#                     pred[target[k]] = source[k]
#                 id+=1
#         for j in n:
#             # Find all connected paths which are profitable
#             if dist[topo[j]] != float("inf"):
#                 current = topo[j]
#                 path.clear()
#                 while current!=i:
#                     if nodepath.count(current) > 0:
#                         deq_it = nodepath[current].rbegin()
#                         while deq_it != nodepath[current].rend():
#                             path.push_front(deref(deq_it))
#                             inc(deq_it)
#                         break
#                     path.push_front(current)
#                     current = pred[current]
#                 nodepath[topo[j]] = path
#                 path.push_front(i)                          ## Found the best path between two co-linear nodes
#
#
#                 ## Calculate score of the identified path
#                 ascore = float(alen[path[0]])
#                 bscore = float(blen[path[0]])
#                 agap = 0
#                 bgap = 0
#
#                 if path.size() > 1:
#                     for k in range(1, <Py_ssize_t> path.size()):
#                         ascore += float(alen[path[k]])
#                         bscore += float(blen[path[k]])
#                         agap += 0 if 0 > (astart[path[k]] - aend[path[k-1]]) else float(astart[path[k]] - aend[path[k-1]])
#                         if not isinv:
#                             bgap += 0 if 0 > (bstart[path[k]] - bend[path[k-1]]) else float(bstart[path[k]] - bend[path[k-1]])
#                         else:
#                             bgap += 0 if 0 > (bstart[path[k-1]] - bend[path[k]]) else float(bstart[path[k-1]] - bend[path[k]])
#
#                 score = min(((ascore - agap)/ascore),((bscore - bgap)/bscore))
#
#
#                 # print([path[id] for id in range(path.size())], score)
#                 ## Check if the alignments of the path explain sufficient unique alignments and are not excessively overlapped
#
#                 # Only process paths for which the alignments explain more than gaps
#
#                 if score > 0:
#                     al = 0
#                     bl = 0
#                     au = 0
#                     bu = 0
#                     for k in range(<Py_ssize_t> path.size()):
#                         al += almntdata[path[k]][0]
#                         bl += almntdata[path[k]][1]
#                         au += almntdata[path[k]][2]
#                         bu += almntdata[path[k]][3]
#
#                 #Trans block is selected IFF either the unique region on any genome is larger than tUC
#                 # or length of unique region on a genome is larger than tUP times the length of
#                 # the overlapping region on that genome
#                     if au > tUC or bu > tUC or au > tUP*al or bu > tUP*bl:
#                         out.append([path[k] for k in range(<Py_ssize_t> path.size())])
#     return out


def blocksdata(outPlaceBlocks, inPlaceBlocks, threshold, tUC, tUP, chromo, tdgl):
    logger = logging.getLogger('blocksdata.'+ chromo)
    orderedBlocks = outPlaceBlocks[outPlaceBlocks.bDir == 1]
    invertedBlocks = outPlaceBlocks[outPlaceBlocks.bDir == -1]

    logger.debug('Number of directed alignments: '+ str(orderedBlocks.shape[0]) + '. Number of inverted alignments: '+ str(orderedBlocks.shape[0]))

    if len(orderedBlocks) > 0:
        transBlocksNeighbours = getTransSynOrientation(inPlaceBlocks, orderedBlocks, threshold)
        outOrderedBlocks = makeBlocksTree(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, threshold, transBlocksNeighbours[0].values, transBlocksNeighbours[1].values, tdgl)
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
        outInvertedBlocks = makeBlocksTree(invertedCoords.aStart.values, invertedCoords.aEnd.values, invertedCoords.bStart.values, invertedCoords.bEnd.values, threshold, invTransBlocksNeighbours[0].values, invTransBlocksNeighbours[1].values, tdgl)

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



@cython.boundscheck(False)
@cython.wraparound(False)
cdef greedySubsetSelector2(long[:] cluster, transBlocksData, long[:] seedBlocks, cpp_map[long, cpp_vec[long]] agroup, cpp_map[long, cpp_vec[long]] bgroup, int threshold):
    np.random.seed(1)

    cdef:
        ssize_t                                         i, j, k
        cpp_bool                                        fnd, outchanged, changed
        long                                            n = len(transBlocksData)
        long                                            ncls = len(cluster)
        long                                            length=0, ntmp=0            # number of temp cluster still need to be classified
        unsigned long                                   newblock
        unsigned long[:]                                astart, aend, bstart, bend, aindex, bindex, clstrsize
        unsigned short int[:]                           meclass
        unsigned long[:]                                meto, mealist, meblist
        unsigned short int                              auni, buni, status
        unsigned short int[:]                           tempcluster, outblocks, skiplist
        unsigned short int[:]                           intrlist

    bestScore = 0
    bestComb = []

    astart = np.array([transBlocksData[i].aStart for i in range(n)], dtype=np.uint)
    aend = np.array([transBlocksData[i].aEnd for i in range(n)], dtype=np.uint)
    bstart = np.array([transBlocksData[i].bStart for i in range(n)], dtype=np.uint)
    bend = np.array([transBlocksData[i].bEnd for i in range(n)], dtype=np.uint)
    aindex = np.array([transBlocksData[i].transGroupIndices[0] for i in range(n)], dtype=np.uint)
    bindex = np.array([transBlocksData[i].transGroupIndices[1] for i in range(n)], dtype=np.uint)

    tempcluster = np.zeros(len(transBlocksData), dtype=np.uint16)
    outblocks = np.zeros(len(transBlocksData), dtype=np.uint16)
    skiplist = np.zeros(len(transBlocksData), dtype=np.uint16)
    intrlist = np.zeros(len(transBlocksData), dtype=np.uint16)

    for i in range(ncls):
        tempcluster[cluster[i]] = 1
        length+=1

    ntmp = length

    for i in range(len(seedBlocks)):
        outblocks[seedBlocks[i]] = 1
        tempcluster[seedBlocks[i]] = 0
        ntmp-=1

    transBlocksScore = {}
    for i in range(n):
        if tempcluster[i] == 1:
            transBlocksScore[i] = aend[i] - astart[i] + bend[i] - bstart[i] + 2


    garb = deque()
    for i in range(n):
        if not transBlocksData[i].aUni and not transBlocksData[i].bUni:
            garb.append(0)
        elif transBlocksData[i].status == 1:
            garb.append(0)
        elif not transBlocksData[i].aUni:
            garb.append(1)
        elif not transBlocksData[i].bUni:
            garb.append(2)
        else:
            garb.append(3)
    meclass = np.array(list(garb), np.uint16)

    while ntmp > 0:
        outchanged = True
        changed = True
        while ntmp != length:
            length = ntmp

            if changed or outchanged:
                changed = False
            ## Check whether any of the remaining elements are mutually exclusive to the already selected candidates.
            ## Remove any such mutually exclusive element
                for i in range(n):
                    if tempcluster[i] == 0:
                        continue
                    ## when the tempcluster element is overlapping with inplace blocks at one of the genomes
                    if meclass[i] == 1 or meclass[i] == 2:
                        meto = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]], bgroup[bindex[i]], i, threshold, meclass[i])
                        for j in meto:
                            if outblocks[j] == 1:
                                tempcluster[i] = 0
                                ntmp-=1
                                skiplist[i]=1
                                changed = True
                                break
                    ## when the tempcluster element does not overlap with inplace blocks
                    elif meclass[i] == 3:
                        mealist, meblist = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]],  bgroup[bindex[i]], i, threshold, meclass[i])
                        ## remove candidate when at least ME candidate on both genomes has already been selected
                        for j in mealist:
                            if outblocks[j] == 1:
                                for k in meblist:
                                    if outblocks[k] == 1:
                                        tempcluster[i] = 0
                                        ntmp-=1
                                        skiplist[i]=1
                                        changed = True
                                        break
                                break

            if changed or outchanged:
                outchanged = False
                changed = False
            ## For a given candidate, if all of its ME candidates have already been added to the skiplist, then that candidate will be selected as part of output
                for i in range(n):
                    if tempcluster[i] ==0:
                        continue
                    ## when the tempcluster element is overlapping with inplace blocks at one of the genomes
                    if meclass[i] == 1 or meclass[i] == 2:
                        meto = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]], bgroup[bindex[i]], i, threshold, meclass[i])
                        fnd = True
                        for j in meto:
                            if skiplist[j] == 0:
                                fnd = False
                                break
                        if fnd:
                            tempcluster[i]=0
                            ntmp-=1
                            outblocks[i]=1
                            changed = True
                    ## when the tempcluster element does not overlap with inplace blocks
                    elif meclass[i] == 3:
                        mealist, meblist = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]], bgroup[bindex[i]], i, threshold, meclass[i])
                        fnd = True
                        for j in mealist:
                            if skiplist[j]==0:
                                fnd = False
                                break
                        if fnd:
                            for j in meblist:
                                if skiplist[j] ==0:
                                    fnd = False
                                    break
                        if fnd:
                            tempcluster[i] = 0
                            ntmp-=1
                            outblocks[i] = 1
                            changed = True

        if ntmp>0:
            ## select one of the twenty highest scoring candidate randomly to break the deadlock
            topblocks = deque()
            for i in range(n):
                if tempcluster[i] == 1:
                    topblocks.append(i)
            topblocks = sorted(list(topblocks), key = lambda x: transBlocksScore[x], reverse = True)[:20]
            totalscore = sum(transBlocksScore[i] for i in topblocks)
            prob = [transBlocksScore[i]/totalscore for i in topblocks]
            i = int(np.random.choice(topblocks, size = 1, p = prob))

            ## update lists to include the selected candidate
            tempcluster[i] = 0
            ntmp-=1
            outblocks[i] = 1

            ## Remove ME candidates
            if meclass[i] == 1 or meclass[i] == 2:
                meto = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]], bgroup[bindex[i]], i, threshold, meclass[i])
                for j in meto:
                    if tempcluster[j] == 1:
                        tempcluster[j] = 0
                        ntmp-=1
                    skiplist[j] = 1

            elif meclass[i] == 3:
                mealist, meblist = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]], bgroup[bindex[i]], i, threshold, meclass[i])
                fnd = False
                for j in mealist:
                    if outblocks[j] == 1:
                        fnd = True
                        for k in meblist:
                            if tempcluster[k] == 1:
                                tempcluster[k] = 0
                                ntmp -=1
                            skiplist[k] = 1
                        break
                if not fnd:
                    for j in meblist:
                        if outblocks[j] == 1:
                            for k in mealist:
                                if tempcluster[k] == 1:
                                    tempcluster[k] = 0
                                    ntmp -=1
                                skiplist[k] = 1
                            break

                ## Also remove candidates which are ME to the selected candidate on both genomes
                for j in range(n):
                    intrlist[j] = 0

                for j in mealist:
                    intrlist[j] = 1

                for j in meblist:
                    if intrlist[j] != 0:
                        if tempcluster[j] == 1:
                            tempcluster[j] = 0
                            ntmp -= 1
                        skiplist[j] = 1

    bestScore, bestComb = updateBestComb(bestScore, bestComb, np.nonzero(outblocks)[0], transBlocksData)
    return bestScore, bestComb


@cython.boundscheck(False)
@cython.wraparound(False)
cdef greedySubsetSelectorHeuristic(long[:] cluster, transBlocksData, long[:] seedBlocks, cpp_map[long, cpp_vec[long]] agroup, cpp_map[long, cpp_vec[long]] bgroup, int threshold):
    np.random.seed(1)
    cdef:
        Py_ssize_t                                      i, j, k, l
        cpp_bool                                        fnd, askip, bskip
        long                                            n = len(transBlocksData)
        long                                            ncls = len(cluster)
        long                                            length=0, ntmp=0            # number of temp cluster still need to be classified
        unsigned long[:]                                astart, aend, bstart, bend, aindex, bindex, clstrsize
        unsigned short int[:]                           meclass
        unsigned long[:]                                meto, mealist, meblist
        unsigned short int                              auni, buni, status
        unsigned short int[:]                           tempcluster, outblocks, skiplist
        unsigned short int[:]                           intrlist
        long[:]                                         outblockindex, skipindexmap

    bestScore = 0
    bestComb = []
    astart = np.array([transBlocksData[i].aStart for i in range(n)], dtype=np.uint)
    aend = np.array([transBlocksData[i].aEnd for i in range(n)], dtype=np.uint)
    bstart = np.array([transBlocksData[i].bStart for i in range(n)], dtype=np.uint)
    bend = np.array([transBlocksData[i].bEnd for i in range(n)], dtype=np.uint)
    aindex = np.array([transBlocksData[i].transGroupIndices[0] for i in range(n)], dtype=np.uint)
    bindex = np.array([transBlocksData[i].transGroupIndices[1] for i in range(n)], dtype=np.uint)

    garb = deque()
    for i in range(n):
        if not transBlocksData[i].aUni and not transBlocksData[i].bUni:
            garb.append(0)
        elif transBlocksData[i].status == 1:
            garb.append(0)
        elif not transBlocksData[i].aUni:
            garb.append(1)
        elif not transBlocksData[i].bUni:
            garb.append(2)
        else:
            garb.append(3)
    meclass = np.array(list(garb), np.uint16)
    ntmp = 0
    tempcluster = np.zeros(n, dtype=np.uint16)
    outblocks = np.zeros(n, dtype=np.uint16)
    skiplist = np.zeros(n, dtype=np.uint16)
    intrlist = np.zeros(n, dtype=np.uint16)
    for i in range(ncls):
        tempcluster[cluster[i]] = 1
        ntmp+=1
    for i in range(len(seedBlocks)):
        outblocks[seedBlocks[i]] = 1
        tempcluster[seedBlocks[i]] = 0
        ntmp-=1
    transBlocksScore = {}
    for i in range(n):
        if tempcluster[i] == 1:
            transBlocksScore[i] = aend[i] - astart[i] + bend[i] - bstart[i] + 2

    outblockindex = np.where(np.array(outblocks) == 1)[0]
    for i in range(n):
        if tempcluster[i] == 0:
            continue
        ## when the tempcluster element is overlapping with inplace blocks at one of the genomes
        if meclass[i] == 1:
            for j in outblockindex:
                if outblocks[j] == 1:
                    if bstart[j] - threshold < bstart[i]:
                        if bend[j] + threshold > bend[i]:
                            tempcluster[i] = 0
                            ntmp-=1
                            skiplist[i]=1
                            break
        elif meclass[i] == 2:
            for j in outblockindex:
                if outblocks[j] == 1:
                    if astart[j] - threshold < astart[i]:
                        if aend[j] + threshold > aend[i]:
                            tempcluster[i] = 0
                            ntmp-=1
                            skiplist[i]=1
                            break
        elif meclass[i] == 3:
            askip = False
            bskip = False
            for j in outblockindex:
                if outblocks[j] == 1:
                    if not askip:
                        if astart[j] - threshold < astart[i]:
                            if aend[j] + threshold > aend[i]:
                                askip = True
                    if not bskip:
                        if bstart[j] - threshold < bstart[i]:
                            if bend[j] + threshold > bend[i]:
                                bskip = True
                    if askip and bskip:
                        tempcluster[i] = 0
                        ntmp-=1
                        skiplist[i]=1
                        break

    skipindexmap = np.zeros(n, np.int)      ## smallest index for a tempcluster for which to check the skiplist
    ## For a given candidate, if all of its ME candidates have already been added to the skiplist, then that candidate will be selected as part of output
    for i in range(n):
        if tempcluster[i] ==0:
            continue
        fnd = True
        for j in range(skipindexmap[i], n):
            if skiplist[j] == 0:
                ## when the tempcluster element is overlapping with inplace blocks at one of the genomes
                if meclass[i] == 1:
                    if bstart[j] - threshold < bstart[i]:
                        if bend[j] + threshold > bend[i]:
                            fnd = False
                            skipindexmap[i] = j
                            break
                elif meclass[i] == 2:
                    if astart[j] - threshold < astart[i]:
                        if aend[j] + threshold > aend[i]:
                            fnd = False
                            skipindexmap[i] = j
                            break
                ## when the tempcluster element does not overlap with inplace blocks
                elif meclass[i]==3:
                    if astart[j] - threshold < astart[i]:
                        if aend[j] + threshold > aend[i]:
                            fnd = False
                            skipindexmap[i] = j
                            break
                    if bstart[j] - threshold < bstart[i]:
                        if bend[j] + threshold > bend[i]:
                            fnd = False
                            skipindexmap[i] = j
                            break
        if fnd:
            tempcluster[i]=0
            ntmp-=1
            outblocks[i]=1
            skipindexmap[i] = n

    while ntmp>0:
        ## select one of the twenty highest scoring candidate randomly to break the deadlock
        topblocks = deque()
        for i in range(n):
            if tempcluster[i] == 1:
                topblocks.append(i)
        topblocks = sorted(list(topblocks), key = lambda x: transBlocksScore[x], reverse = True)[:20]
        totalscore = sum(transBlocksScore[i] for i in topblocks)
        prob = [transBlocksScore[i]/totalscore for i in topblocks]
        i = int(np.random.choice(topblocks, size = 1, p = prob))

        ## update lists to include the selected candidate
        tempcluster[i] = 0
        ntmp-=1
        outblocks[i] = 1

        ## Remove ME candidates
        if meclass[i] == 1 or meclass[i] == 2:
            meto = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]], bgroup[bindex[i]], i, threshold, meclass[i])
            for j in meto:
                if tempcluster[j] == 1:
                    tempcluster[j] = 0
                    ntmp-=1
                skiplist[j] = 1
        elif meclass[i] == 3:
            mealist, meblist = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[i]], bgroup[bindex[i]], i, threshold, meclass[i])
            fnd = False
            for j in mealist:
                if outblocks[j] == 1:
                    fnd = True
                    for k in meblist:
                        if tempcluster[k] == 1:
                            tempcluster[k] = 0
                            ntmp -=1
                        skiplist[k] = 1
                    break
            if not fnd:
                for j in meblist:
                    if outblocks[j] == 1:
                        for k in mealist:
                            if tempcluster[k] == 1:
                                tempcluster[k] = 0
                                ntmp -=1
                            skiplist[k] = 1
                        break
            ## Also remove candidates which are ME to the selected candidate on both genomes
            for j in range(n):
                intrlist[j] = 0

            for j in mealist:
                intrlist[j] = 1

            for j in meblist:
                if intrlist[j] != 0:
                    if tempcluster[j] == 1:
                        tempcluster[j] = 0
                        ntmp -= 1
                    skiplist[j] = 1


        outblockindex = np.where(np.array(outblocks) == 1)[0]
        ## Check which of the candidates in the same group can't be added now
        for j in agroup[aindex[i]]:
            if tempcluster[j] == 0:
                continue
            ## when the tempcluster element is overlapping with inplace blocks at one of the genomes
            if meclass[j] == 1:
                for k in outblockindex:
                    if outblocks[k] == 1:
                        if bstart[k] - threshold < bstart[j]:
                            if bend[k] + threshold > bend[j]:
                                tempcluster[j] = 0
                                ntmp-=1
                                skiplist[j]=1
                                break
            elif meclass[j] == 2:
                for k in outblockindex:
                    if outblocks[k] == 1:
                        if astart[k] - threshold < astart[j]:
                            if aend[k] + threshold > aend[j]:
                                tempcluster[j] = 0
                                ntmp-=1
                                skiplist[j]=1
                                break
            elif meclass[j] == 3:
                askip = False
                bskip = False
                for k in outblockindex:
                    if outblocks[k] == 1:
                        if not askip:
                            if astart[k] - threshold < astart[j]:
                                if aend[k] + threshold > aend[j]:
                                    askip = True
                        if not bskip:
                            if bstart[k] - threshold < bstart[j]:
                                if bend[k] + threshold > bend[j]:
                                    bskip = True
                        if askip and bskip:
                            tempcluster[j] = 0
                            ntmp-=1
                            skiplist[j]=1
                            break

        for j in bgroup[bindex[i]]:
            if tempcluster[j] == 0:
                continue
            ## when the tempcluster element is overlapping with inplace blocks at one of the genomes
            if meclass[j] == 1:
                for k in outblockindex:
                    if outblocks[k] == 1:
                        if bstart[k] - threshold < bstart[j]:
                            if bend[k] + threshold > bend[j]:
                                tempcluster[j] = 0
                                ntmp-=1
                                skiplist[j]=1
                                break
            elif meclass[j] == 2:
                for k in outblockindex:
                    if outblocks[k] == 1:
                        if astart[k] - threshold < astart[j]:
                            if aend[k] + threshold > aend[j]:
                                tempcluster[j] = 0
                                ntmp-=1
                                skiplist[j]=1
                                break
            elif meclass[j] == 3:
                askip = False
                bskip = False
                for k in outblockindex:
                    if outblocks[k] == 1:
                        if not askip:
                            if astart[k] - threshold < astart[j]:
                                if aend[k] + threshold > aend[j]:
                                    askip = True
                        if not bskip:
                            if bstart[k] - threshold < bstart[j]:
                                if bend[k] + threshold > bend[j]:
                                    bskip = True
                        if askip and bskip:
                            tempcluster[j] = 0
                            ntmp-=1
                            skiplist[j]=1
                            break
    bestScore, bestComb = updateBestComb(bestScore, bestComb, np.nonzero(outblocks)[0], transBlocksData)
    return bestScore, bestComb



@cython.boundscheck(False)
@cython.wraparound(False)
cdef greedySubsetSelector(cluster, transBlocksData, seedblocks, iterCount = 100):
    cdef:
        Py_ssize_t                                      i, j, k, l
        cpp_bool                                        fnd, outchanged, changed
        long                                   n = len(transBlocksData)
        long                                   ncls = len(cluster)
        long                                            length=0, ntmp=0            # number of temp cluster still need to be classified
        unsigned long                                   newblock
        unsigned long[:]                                astart, aend, bstart, bend, aindex, bindex, clstrsize
        unsigned short int[:]                           meclass
        unsigned long[:]                                meto, mealist, meblist
        unsigned short int                              auni, buni, status
        unsigned short int[:]                           tempcluster, outblocks, skiplist
        unsigned short int[:]                           intrlist


    np.random.seed(1)
    bestScore = 0
    bestComb = []

    ## initiate transBlocksScore
    tempcluster = np.zeros(n, dtype=np.uint16)
    outblocks = np.zeros(n, dtype=np.uint16)
    skiplist = np.zeros(n, dtype=np.uint16)
    for i in range(ncls):
        tempcluster[cluster[i]] = 1
    for i in range(len(seedblocks)):
        outblocks[seedblocks[i]] = 1
        tempcluster[seedblocks[i]] = 0
    transBlocksScore = {}
    for i in np.nonzero(tempcluster)[0].astype(np.int):
        transBlocksScore[i] = (transBlocksData[i].aEnd - transBlocksData[i].aStart) + (transBlocksData[i].bEnd - transBlocksData[i].bStart)

    # Find best cluster subset
    for i in range(iterCount):

        # Initialise arrays
        ntmp = 0
        tempcluster = np.zeros(n, dtype=np.uint16)
        outblocks = np.zeros(n, dtype=np.uint16)
        skiplist = np.zeros(n, dtype=np.uint16)
        for i in range(ncls):
            tempcluster[cluster[i]] = 1
            ntmp+=1
        for i in range(len(seedblocks)):
            outblocks[seedblocks[i]] = 1
            tempcluster[seedblocks[i]] = 0
            ntmp-=1
        length = 0

        while ntmp > 0:
            # Run the loop as long as there are some changes happening
            while ntmp != length:
                length = ntmp
                # Remove blocks for which corresponding mutually exclusive elements have already been selected
                for j in range(n):
                    if tempcluster[j] == 0:
                        continue
                    if hasattr(transBlocksData[j],"meTo"):
                        for k in transBlocksData[j].meTo:
                            if outblocks[k] == 1:
                                tempcluster[j] = 0
                                ntmp-=1
                                skiplist[j]=1
                    elif hasattr(transBlocksData[j], "meAlist"):
                        for k in transBlocksData[j].meAlist:
                            if outblocks[k] == 1:
                                for l in transBlocksData[j].meBlist:
                                    if outblocks[l] == 1:
                                        tempcluster[j] = 0
                                        ntmp-=1
                                        skiplist[j] = 1
                                        break
                                break

                # Select blocks for which corresponding mutually exclusive elements have already been rejected
                for j in range(n):
                    if tempcluster[j]==0:
                        continue
                    if hasattr(transBlocksData[j],"meTo"):
                        fnd = False
                        for k in transBlocksData[j].meTo:
                            if skiplist[k] == 0:
                                fnd = True
                                break
                        if not fnd:
                            tempcluster[j] = 0
                            ntmp-=1
                            outblocks[j] = 1
                    elif hasattr(transBlocksData[j], "meAlist"):
                        fnd = False
                        for k in transBlocksData[j].meAlist:
                            if skiplist[k] == 0:
                                fnd = True
                                break
                        if not fnd:
                            for k in transBlocksData[j].meBlist:
                                if skiplist[k] == 0:
                                    fnd = True
                                    break
                        if not fnd:
                            tempcluster[j] = 0
                            ntmp-=1
                            outblocks[j] = 1

            # Select one of the top 20 blocks to break the deadlock
            if ntmp > 0:
                topBlocks = sorted(np.nonzero(tempcluster)[0], key = lambda x: transBlocksScore[x], reverse = True)[:20]
                totalScore = sum(transBlocksScore[i] for i in topBlocks)
                prob = [transBlocksScore[i]/totalScore for i in topBlocks]
                newblock = int(np.random.choice(topBlocks, size = 1, p = prob))
                outblocks[newblock] = 1
                tempcluster[newblock] = 0
                ntmp-=1

                # Remove blocks contradicting the selected block
                if hasattr(transBlocksData[newblock],"meTo"):
                    for j in transBlocksData[newblock].meTo:
                        if tempcluster[j] == 1:
                            tempcluster[j] = 0
                            ntmp-=1
                        skiplist[j] = 1

                elif hasattr(transBlocksData[newblock],"meAlist"):
                    fnd = False
                    for j in transBlocksData[newblock].meAlist:
                        if outblocks[j] == 1:
                            fnd = True
                            for k in transBlocksData[newblock].meBlist:
                                if tempcluster[k] == 1:
                                    tempcluster[k] = 0
                                    ntmp-=1
                                skiplist[k] = 1
                            break
                    if not fnd:
                        for j in transBlocksData[newblock].meBlist:
                            if outblocks[j] == 1:
                                for k in transBlocksData[newblock].meAlist:
                                    if tempcluster[k] == 1:
                                        tempcluster[k] = 0
                                        ntmp-=1
                                    skiplist[k] = 1
                                break
                    for j in transBlocksData[newblock].meAlist:
                        if j in transBlocksData[newblock].meBlist:
                            if tempcluster[j] == 1:
                                tempcluster[j] = 0
                                ntmp-=1
                            skiplist[j] = 1
        bestScore, bestComb = updateBestComb(bestScore, bestComb, np.nonzero(outblocks)[0], transBlocksData)
    return(bestScore, bestComb)


def getBestClusterSubset(cluster, transBlocksData, bRT, chromo='', aGroups=None, bGroups=None, threshold=None):
    logger = logging.getLogger('tdcluster'+chromo)
    if len(cluster) == 0:
        return
    seedBlocks = [i for i in cluster if transBlocksData[i].status == 1]
    if len(cluster) > 100000:
        logger.info('Massive (>100000 candidates) TD cluster (with '+ str(len(cluster)) +' candidate TDs) identified. Using low-memory high-runtime approach. Iterative sampling disabled. Using less stringent progressive elimination.')
        output = greedySubsetSelectorHeuristic(np.array(cluster, np.int), transBlocksData, np.array(seedBlocks, np.int), aGroups, bGroups, threshold)
    elif len(cluster) > 10000:
        logger.info('Large (>10000 candidates) TD cluster (with '+ str(len(cluster)) +' candidate TDs) identified. Using low-memory high-runtime approach. Iterative sampling disabled.')
        output = greedySubsetSelector2(np.array(cluster, np.int), transBlocksData, np.array(seedBlocks, np.int), aGroups, bGroups, threshold)
    elif len(cluster) < 50:
        output = bruteSubsetSelector(cluster, transBlocksData, seedBlocks, bRT)
        if output == "Failed":
            output = greedySubsetSelector(cluster, transBlocksData, seedBlocks)
    else:
        output = greedySubsetSelector(cluster, transBlocksData, seedBlocks)
    return output


def getTransClasses(clusterSolutionBlocks, transData, transagroups, transbgroups, astart, aend, bstart, bend, aindex, bindex, agroup, bgroup, threshold, meclass):
    logger = logging.getLogger('gettransclasses')
    def settl(j):
        if transData[j].dir == 1:
            transClasses["translocation"].append(j)
        elif transData[j].dir == -1:
            transClasses["invTranslocation"].append(j)
        else:
            logger.info("Wrong alignment direction" + j)

    def setdup(j):
        if transData[j].dir == 1:
             transClasses["duplication"].append(j)
        elif transData[j].dir == -1:
            transClasses["invDuplication"].append(j)
        else:
            logger.info("Wrong alignment direction" + j)

    transClasses = {"translocation":[],
                    "invTranslocation":[],
                    "duplication":[],
                    "invDuplication":[]}

    for i in clusterSolutionBlocks:
        for j in i:
            if not transData[j].aUni and not transData[j].bUni:
                logger.error("Redundant candidate selected as TD" + str(j))
            elif transData[j].status == 1:
                if not transData[j].aUni or not transData[j].bUni:
                    setdup(j)
                elif transData[j].aUni and transData[j].bUni:
                    if transData[j].genomeAUni and transData[j].genomeBUni:
                        settl(j)
                    elif not transData[j].genomeAUni:
                        istrans = 1
                        for k in transData[j].getoverlappingregions(groups = transagroups,genome="a"):
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
                        for k in transData[j].getoverlappingregions(groups = transbgroups,genome="b"):
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
                    if len(np.intersect1d(transData[j].meAlist, i)) > 0 or len(np.intersect1d(transData[j].meBlist, i)) > 0:
                        setdup(j)
                    else:
                        settl(j)
                else:
                    if meclass[j] == 1 or meclass[j] == 2:
                        meto = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[j]], bgroup[bindex[j]], j, threshold, meclass[j])
                        if len(np.intersect1d(meto, i)) > 0:
                            setdup(j)
                        else:
                            settl(j)
                    elif meclass[j] == 3:
                        mealist, meblist = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[j]],  bgroup[bindex[j]], j, threshold, meclass[j])
                        if len(np.intersect1d(mealist, i)) > 0 or len(np.intersect1d(meblist, i)) > 0:
                            setdup(j)
                        else:
                            settl(j)
                    else:
                        logger.info("Wrong candidate class" + j)
    return transClasses



def getDupGenome(dupData, allTransBlocksData, transClasses, astart, aend, bstart, bend, aindex, bindex, agroup, bgroup, threshold, meclass):
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
        if meclass[row.Index] != 3:
            print('Error in dup class identification ', row.Index)
            continue
        mealist, meblist = getmeblocks2(astart, aend, bstart, bend, agroup[aindex[row.Index]],  bgroup[bindex[row.Index]], row.Index, threshold, meclass=3)
        for i in mealist:
            if i in transClasses["translocation"] or i in transClasses["invTranslocation"]:
                found = True
                dupGenomes.append("B")
                break
        if not found:
            dupGenomes.append("A")
    dupData["dupGenomes"] = pd.Series(dupGenomes, index = dupData.index)
    return(dupData)

