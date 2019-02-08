# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:54:53 2017

@author: goel
"""

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

cimport numpy as np
cimport cython

np.random.seed(1)

from syri.pyxFiles.function cimport getmeblocks, getOverlapWithSynBlocks

def readCoords(coordsfin, chrmatch, cigar = False):
    logger = logging.getLogger('Reading Coords')
    try:
        coords = pd.read_table(coordsfin, header = None)
    except pd.errors.ParserError:
        coords = pd.read_table(coordsfin, header = None, engine = "python")

    if not cigar:
        if coords.shape[1] >= 12:
            coords = coords.iloc[:, 0:11]
        coords.columns = ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr"]
    else:
        if coords.shape[1] > 12:
            coords = coords.iloc[:, 0:12]
        coords.columns = ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']


    # Sanity check input file
    try:
        coords.aStart = coords.aStart.astype('int')
    except ValueError:
        logger.error('astart is not int')
        sys.exit()

    try:
        coords.aEnd = coords.aEnd.astype('int')
    except ValueError:
        logger.error('aend is not int')
        sys.exit()

    try:
        coords.bStart = coords.bStart.astype('int')
    except ValueError:
        logger.error('bstart is not int')
        sys.exit()

    try:
        coords.bEnd = coords.bEnd.astype('int')
    except ValueError:
        logger.error('abend is not int')
        sys.exit()

    try:
        coords.aLen = coords.aLen.astype('int')
    except ValueError:
        logger.error('alen is not int')
        sys.exit()

    try:
        coords.bLen = coords.bLen.astype('int')
    except ValueError:
        logger.error('blen is not int')
        sys.exit()

    try:
        coords.iden = coords.iden.astype('float')
    except ValueError:
        logger.error('iden is not float')
        sys.exit()

    try:
        coords.aDir = coords.aDir.astype('int')
    except ValueError:
        logger.error('aDir is not int')
        sys.exit()

    if any(coords.aDir != 1):
        logger.error('aDir can only have values 1')
        sys.exit()

    try:
        coords.bDir = coords.bDir.astype('int')
    except ValueError:
        logger.error('bDir is not int')
        sys.exit()

    for i in coords.bDir:
        if i not in [1,-1]:
            logger.error('bDir can only have values 1/-1')
            sys.exit()

    #check for bstart > bend when bdir is -1
    check = np.unique(coords.loc[coords.bDir == -1, 'bStart'] > coords.loc[coords.bDir == -1, 'bEnd'])
    if len(check) > 1:
        logger.error('Inconsistent start and end position for inverted alignment in genome B. For inverted alignments, either all bstart < bend or all bend > bstart')
        sys.exit()
    elif check[0] == True:
        pass
    else:
        logger.warning('For inverted alignments, bstart was less than bend. Swapping them.')
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] + coords.loc[coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bEnd'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[coords.bDir == -1, 'bEnd']

    coords.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)

    ## Ensure that chromosome IDs are same for the two genomes.
    ## Either find best query match for every reference genome.
    ## Or if --no-chrmatch is set then remove non-matching chromosomes.

    if np.unique(coords.aChr).tolist() != np.unique(coords.bChr).tolist():
        logger.warning('Chromosomes IDs do not match.')
        if not chrmatch:
            if len(np.unique(coords.aChr)) != len(np.unique(coords.bChr)):
                logger.error("Unequal number of chromosomes in the genomes. Exiting")
                sys.exit()
            else:
                logger.warning("Matching them automatically. For each reference genome, most similar query genome will be selected")
                chromMaps = defaultdict(dict)
                for i in np.unique(coords.bChr):
                    for j in np.unique(coords.aChr):
                        a = np.array(coords.loc[(coords.bChr == i) & (coords.aChr == j), ["aStart", "aEnd"]])
                        a = mergeRanges(a)
                        chromMaps[j][i] = len(a) + (a[:, 1] - a[:, 0]).sum()

                assigned = []
                for chrom in np.unique(coords.aChr):
                    maxid = max(chromMaps[chrom].items(), key=lambda x: x[1])[0]
                    if maxid in assigned:
                        logger.error("{} in genome B is best match for two chromosomes in genome A. Cannot assign chromosomes automatically.".format(maxid))
                        sys.exit()
                    assigned.append(maxid)
                    logger.info("setting {} as {}".format(maxid, chrom))
                    coords.loc[coords.bChr == maxid, "bChr"] = chrom
        else:
            logger.warning("--no-chrmatch is set. Not matching chromosomes automatically.")
            aChromo = set(coords["aChr"])
            bChromo = set(coords["bChr"])
            badChromo = list(aChromo - bChromo) + list(bChromo - aChromo)
            logger.warning(", ".join(badChromo) + " present in only one genome. Removing corresponding alignments")
            coords = coords.loc[~coords.aChr.isin(badChromo) & ~coords.bChr.isin(badChromo)]

    return coords

def startSyri(args):
    coordsfin = args.infile.name
    nCores = args.nCores
    bRT = args.bruteRunTime
    threshold = 50  ##args.threshold
    cwdPath = args.dir
    prefix = args.prefix
    tUC = args.TransUniCount
    tUP = args.TransUniPercent
    chrmatch = args.chrmatch

    # LOG_FORMAT = "%(asctime)s — %(name)s — %(levelname)s — %(funcName)s:%(lineno)d — %(message)s"
    # logging.basicConfig(filename=cwdPath+args.log_fin.name, level=args.log, format=LOG_FORMAT)
    logger = logging.getLogger("syri")
    logger.warning("starting")
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    coords = readCoords(coordsfin, chrmatch)

    uniChromo = list(np.unique(coords.aChr))
    logger.info('Analysing chromosomes: {}'.format(uniChromo))
    # Identify intra-chromosomal events (synteny, inversions, intra-trans, intra-dup) for each chromosome as a separate
    # process in parallel
    with Pool(processes = nCores) as pool:
        pool.map(partial(syri,threshold=threshold,coords=coords, cwdPath= cwdPath, bRT = bRT, prefix = prefix, tUC=tUC, tUP=tUP), uniChromo)

    # Merge output of all chromosomes
    mergeOutputFiles(uniChromo,cwdPath, prefix)

    #Identify cross-chromosomal events in all chromosomes simultaneously
    from syri.tdfunc import getCTX
    getCTX(coords, cwdPath, uniChromo, threshold, bRT, prefix, tUC, tUP, nCores)

    # Recalculate syntenic blocks by considering the blocks introduced by CX events
    outSyn(cwdPath, threshold, prefix)
    
def syri(chromo, threshold, coords, cwdPath, bRT, prefix, tUC, tUP):
    logger = logging.getLogger("syri."+chromo)

    coordsData = coords[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == 1)]

    logger.info(chromo+" " + str(coordsData.shape))
    logger.info("Identifying Synteny for chromosome " + chromo)

    df = pd.DataFrame(apply_TS(coordsData.aStart.values,coordsData.aEnd.values,coordsData.bStart.values,coordsData.bEnd.values, threshold), index = coordsData.index.values, columns = coordsData.index.values)
    nrow = df.shape[0]
    blocks = [alignmentBlock(i, np.where(df.iloc[i,] == True)[0], coordsData.iloc[i]) for i in range(nrow)]
    for block in blocks:
        i = 0
        while(i < len(block.children)):
            block.children = list(set(block.children) - set(blocks[block.children[i]].children))
            i+=1
        block.children.sort()

        for child in block.children:
            blocks[child].addParent(block.id)

        scores = [blocks[parent].score for parent in block.parents]
        if len(scores) > 0:
            block.bestParent(block.parents[scores.index(max(scores))], max(scores))
    synPath = getSynPath(blocks)
    synData = coordsData.iloc[synPath].copy()
    del(coordsData, blocks, df)
    collect()

    ##########################################################################
    #   Finding Inversions
    ##########################################################################
    logger.info("Identifying Inversions for chromosome " + chromo)

    from syri.inversions import getInversions
    invertedCoordsOri, profitable, bestInvPath, invData, synInInv, badSyn = getInversions(coords,chromo, threshold, synData, synPath)



    ##########################################################
    #### Identify Translocation and duplications
    ##########################################################
    logger.info("Identifying translocation and duplication for chromosome " + chromo)

    # Import functions
    from syri.tdfunc import getTransSynOrientation, findOrderedTranslocations, mergeTransBlocks, makeTransGroupList, getTransCluster, transBlock, getBestClusterSubset, getTransClasses, getDupGenome

    chromBlocks = coords[(coords.aChr == chromo) & (coords.bChr == chromo)]
    inPlaceIndices = sorted(list(synData.index.values) + list(invData.index.values))
    inPlaceBlocks = chromBlocks[chromBlocks.index.isin(sorted(list(synData.index.values)))].copy()
    
    
    for i in bestInvPath:
        invPos = profitable[i].invPos
        invBlockData = invertedCoordsOri.iloc[invPos]
        invCoord = [invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]]
        invCoord.append(invCoord[1] - invCoord[0])
        invCoord.append(invCoord[3] - invCoord[2])
        invCoord.append(sum((invBlockData.aLen+invBlockData.bLen)*invBlockData.iden)/(invCoord[-2] + invCoord[-1]))
        invCoord.extend([1,-1,chromo,chromo])
        for j in range(profitable[i].neighbours[0]+1,profitable[i].neighbours[1]):
            inPlaceBlocks = inPlaceBlocks[inPlaceBlocks.index != synData.iloc[j].name]
            try:
                inPlaceIndices.remove(synData.iloc[j].name)
            except:
                pass
        inPlaceBlocks = inPlaceBlocks.append(pd.Series(invCoord, index = inPlaceBlocks.columns, name = invPos[0]))
        
    inPlaceBlocks.sort_values(["aChr","aStart","aEnd","bChr","bStart","bEnd"], inplace = True)
    inPlaceBlocks.index = range(inPlaceBlocks.shape[0])
    outPlaceBlocks = chromBlocks[~chromBlocks.index.isin(inPlaceIndices)]
    
    logger.debug("Translocations : found blocks" + chromo)
    ## Should not filter redundant alignments as they "can" be part of bigger translocations
    ## filtering them may lead to removal of those translocations
    
    outPlaceBlocksFiltered = outPlaceBlocks.copy() 
        
    ## Create connectivity tree for directed and inverted blocks
    #### find all translocations which don't have large gaps between its alignments
    #### and are not overlappign with the syntenic blocks
    
    orderedBlocks = outPlaceBlocksFiltered[outPlaceBlocksFiltered.bDir == 1]
    invertedBlocks = outPlaceBlocksFiltered[outPlaceBlocksFiltered.bDir == -1]

    if len(orderedBlocks) > 0:
        transBlocksNeighbours = getTransSynOrientation(inPlaceBlocks, orderedBlocks, threshold)
        outOrderedBlocks = pd.DataFrame(makeBlocksTree(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, orderedBlocks.bDir.values, orderedBlocks.aChr.values, orderedBlocks.bChr.values, orderedBlocks.index.values, threshold, transBlocksNeighbours[0].values, transBlocksNeighbours[1].values))
        transBlocks = findOrderedTranslocations(outOrderedBlocks, orderedBlocks, inPlaceBlocks, threshold, tUC, tUP, ctx = False)
    else:
        transBlocks = []
    
    if len(invertedBlocks) > 0:
        invertedCoords = invertedBlocks.copy()
        invertedCoords.bStart = invertedCoords.bStart + invertedCoords.bEnd
        invertedCoords.bEnd = invertedCoords.bStart - invertedCoords.bEnd
        invertedCoords.bStart = invertedCoords.bStart - invertedCoords.bEnd
        invTransBlocksNeighbours = getTransSynOrientation(inPlaceBlocks, invertedBlocks, threshold)
        invertedCoords = invertedBlocks.copy()
        maxCoords = np.max(np.max(invertedCoords[["bStart","bEnd"]]))
        invertedCoords.bStart = maxCoords + 1 - invertedCoords.bStart 
        invertedCoords.bEnd = maxCoords + 1 - invertedCoords.bEnd
        outInvertedBlocks = pd.DataFrame(makeBlocksTree(invertedCoords.aStart.values, invertedCoords.aEnd.values, invertedCoords.bStart.values, invertedCoords.bEnd.values, invertedCoords.bDir.values, invertedCoords.aChr.values, invertedCoords.bChr.values, invertedCoords.index.values, threshold, invTransBlocksNeighbours[0].values, invTransBlocksNeighbours[1].values))
        invTransBlocks = findOrderedTranslocations(outInvertedBlocks, invertedCoords, inPlaceBlocks, threshold, tUC, tUP,ctx = False)
    else:
        invTransBlocks = []

    logger.debug("Translocations : found orderedBlocks " + chromo)
    logger.debug("Translocations : merging blocks " + chromo)

    allTransBlocks, allTransIndexOrder = mergeTransBlocks(transBlocks, orderedBlocks, invTransBlocks, invertedBlocks)
    allTransGenomeAGroups = makeTransGroupList(allTransBlocks, "aStart", "aEnd", threshold)
    allTransGenomeBGroups = makeTransGroupList(allTransBlocks, "bStart", "bEnd", threshold)
    
    allTransGroupIndices = {}
    for i in range(len(allTransGenomeAGroups)):
        for block in allTransGenomeAGroups[i].member:
            allTransGroupIndices[block] = [i]
    for i in range(len(allTransGenomeBGroups)):
        for block in allTransGenomeBGroups[i].member:
            allTransGroupIndices[block].append(i)
    
    logger.debug("Translocations : getting clusters " + chromo)
    allTransCluster = getTransCluster(allTransGroupIndices, allTransGenomeAGroups, allTransGenomeBGroups)
    
    allTransClusterIndices = dict()
    for i in range(len(allTransCluster)):
        allTransClusterIndices.update(dict.fromkeys(allTransCluster[i], i))

    logger.debug("Translocations : making blocks data " + chromo +" " + str(datetime.now()))
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    if len(allTransBlocks) > 0:
        auni = getOverlapWithSynBlocks(np.array(allTransBlocks.aStart),
                                       np.array(allTransBlocks.aEnd),
                                       np.array([chromo]*allTransBlocks.shape[0]),
                                       np.array(inPlaceBlocks.aStart),
                                       np.array(inPlaceBlocks.aEnd),
                                       np.array([chromo]*inPlaceBlocks.shape[0]),
                                       50,
                                       allTransBlocks.shape[0],
                                       tUC,
                                       tUP)
        sortedInPlace = inPlaceBlocks.sort_values(["bStart","bEnd"])
        buni = getOverlapWithSynBlocks(np.array(allTransBlocks.bStart), np.array(allTransBlocks.bEnd), np.array([chromo]*allTransBlocks.shape[0]), np.array(sortedInPlace.bStart), np.array(sortedInPlace.bEnd), np.array([chromo]*inPlaceBlocks.shape[0]), 50, allTransBlocks.shape[0], tUC, tUP)

    genomeGroupLengths = ([len(i.member) for i in allTransGenomeAGroups], [len(i.member) for i in allTransGenomeBGroups])

    allTransBlocksData = [transBlock(i) for i in range(allTransBlocks.shape[0])]
    count = 0
    for row in allTransBlocks.itertuples(index = False):
        allTransBlocksData[count].aStart = row.aStart
        allTransBlocksData[count].aEnd = row.aEnd
        allTransBlocksData[count].bStart = row.bStart
        allTransBlocksData[count].bEnd = row.bEnd
        allTransBlocksData[count].dir = row.dir
        allTransBlocksData[count].transClusterIndex = allTransClusterIndices[count]
        allTransBlocksData[count].transGroupIndices = allTransGroupIndices[count]
        allTransBlocksData[count].aUni = auni[count]
        allTransBlocksData[count].bUni = buni[count]
        if genomeGroupLengths[0][allTransBlocksData[count].transGroupIndices[0]] == 1:
            allTransBlocksData[count].genomeAUni = True
        if genomeGroupLengths[1][allTransBlocksData[count].transGroupIndices[1]] == 1:
            allTransBlocksData[count].genomeBUni = True
        if (allTransBlocksData[count].aUni and allTransBlocksData[count].genomeAUni) or (allTransBlocksData[count].bUni and allTransBlocksData[count].genomeBUni):
            allTransBlocksData[count].setStatus(1)
        count+=1

    logger.debug("Translocations : finished making blocks data on" + chromo)
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    aUni = np.array([allTransBlocksData[i].aUni for i in range(allTransBlocks.shape[0])], dtype="int")
    bUni = np.array([allTransBlocksData[i].bUni for i in range(allTransBlocks.shape[0])], dtype="int")
    status = np.array([allTransBlocksData[i].status for i in range(allTransBlocks.shape[0])], dtype="int")
    aIndex = np.array([allTransBlocksData[i].transGroupIndices[0] for i in range(allTransBlocks.shape[0])], dtype="int")
    bIndex = np.array([allTransBlocksData[i].transGroupIndices[1] for i in range(allTransBlocks.shape[0])], dtype="int")
    aGroups = {i:np.array(allTransGenomeAGroups[i].member, dtype="int") for i in range(len(allTransGenomeAGroups))}
    bGroups = {i: np.array(allTransGenomeBGroups[i].member, dtype="int") for i in range(len(allTransGenomeBGroups))}

    if len(allTransBlocks) > 0:
        out = getmeblocks(np.array(allTransBlocks.aStart), np.array(allTransBlocks.aEnd), np.array(allTransBlocks.bStart), np.array(allTransBlocks.bEnd), 50, allTransBlocks.shape[0], aUni, bUni, status, aIndex, bIndex, aGroups, bGroups)

        for i in range(len(out[0])):
            if out[0][i]:
                allTransCluster[allTransClusterIndices[i]].remove(i)

        for i in out[1].keys():
            if len(out[1][i]) > 0:
                allTransBlocksData[i].addMEBlock(list(out[1][i]))

        for i in out[2].keys():
            allTransBlocksData[i].setMEList(list(out[2][i][0]),list(out[2][i][1]))

        del(aUni, bUni, status, aIndex, bIndex, aGroups, bGroups, out)
        collect()

    logger.debug("Translocations : finding solutions "+ chromo + str(datetime.now()))
    clusterSolutions = []
    for i in range(len(allTransCluster)):
        if len(allTransCluster[i]) > 0:
            clusterSolutions.append(getBestClusterSubset(allTransCluster[i], allTransBlocksData, bRT))
    
    clusterSolutionBlocks = [i[1] for i in clusterSolutions]
    #clusterBlocks = unlist(clusterSolutionBlocks)
    
    logger.debug("Translocations : processing translocations " + chromo + str(datetime.now()))
    
    transClasses = getTransClasses(clusterSolutionBlocks, allTransBlocksData, allTransGenomeAGroups, allTransGenomeBGroups)
    
    dupData = allTransBlocks.iloc[transClasses["duplication"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    invDupData = allTransBlocks.iloc[transClasses["invDuplication"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    TLData = allTransBlocks.iloc[transClasses["translocation"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    invTLData = allTransBlocks.iloc[transClasses["invTranslocation"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])  
    
    dupData = getDupGenome(dupData, allTransBlocksData, transClasses)
    invDupData = getDupGenome(invDupData, allTransBlocksData, transClasses)
    
    
    fout = open(cwdPath+prefix+chromo+"_invOut.txt","w")
    tempInvBlocks = []
    for i in bestInvPath:
        invPos = profitable[i].invPos
        tempInvBlocks.append([invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]])
        fout.write("\t".join(map(str,["#",invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],"-",invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2],"\n"])))
        for j in invPos:
            fout.write("\t".join(map(str,invertedCoordsOri.iloc[j][:4])))
            fout.write("\n")
    fout.close()
    
    
    ## Grouping Syn blocks : Final synblock identification is done after ctx identification.
    allBlocks, outClusters = groupSyn(tempInvBlocks, dupData, invDupData, invTLData, TLData, threshold, synData, badSyn)
    
########################################################################################################################
    fout = open(cwdPath+prefix+chromo+"_synOut.txt","w")
    for i in outClusters:
        fout.write("\t".join(map(str,["#",allBlocks.at[i[0],"aStart"],allBlocks.at[i[-1],"aEnd"],"-",allBlocks.at[i[0],"bStart"],allBlocks.at[i[-1],"bEnd"],"\n"])))
        for j in i:
            fout.write("\t".join(map(str,allBlocks.loc[j][:-1])))
            if j in synInInv:
                fout.write("\tSyn_in_Inv\n")
            else:
                fout.write("\n")
    fout.close()
########################################################################################################################
    
    fout = open(cwdPath+prefix+chromo+"_dupOut.txt","w")
    for i in dupData.index.values:
        fout.write("\t".join(map(str,["#",dupData.at[i,"aStart"],dupData.at[i,"aEnd"],"-",dupData.at[i,"bStart"],dupData.at[i,"bEnd"],"-", dupData.at[i,"dupGenomes"],"\n"])))
        for j in transBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,orderedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################    
    
    fout = open(cwdPath+prefix+chromo+"_invDupOut.txt","w")
    for i in invDupData.index.values:
        fout.write("\t".join(map(str,["#",invDupData.at[i,"aStart"],invDupData.at[i,"aEnd"],"-",invDupData.at[i,"bStart"],invDupData.at[i,"bEnd"],"-", invDupData.at[i,"dupGenomes"],"\n"])))
        for j in invTransBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,invertedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################
    
    fout = open(cwdPath+prefix+chromo+"_TLOut.txt","w")
    for i in TLData.index.values:
        fout.write("\t".join(map(str,["#",TLData.at[i,"aStart"],TLData.at[i,"aEnd"],"-",TLData.at[i,"bStart"],TLData.at[i,"bEnd"],"\n"])))
        for j in transBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,orderedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################
    
    fout = open(cwdPath+prefix+chromo+"_invTLOut.txt","w")
    for i in invTLData.index.values:
        fout.write("\t".join(map(str,["#",invTLData.at[i,"aStart"],invTLData.at[i,"aEnd"],"-",invTLData.at[i,"bStart"],invTLData.at[i,"bEnd"],"\n"])))
        for j in invTransBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,invertedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################
    


cpdef apply_TS(np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.int threshold):
    assert(aStart.dtype == np.int and aEnd.dtype == np.int and bStart.dtype == np.int and bEnd.dtype == np.int)
    cdef Py_ssize_t i, j,  n = len(aStart)
    assert(n == len(aEnd) == len(bStart) == len(bEnd))
    cdef np.ndarray[object, ndim =2 ] df =  np.array([[np.nan]*n]*n, dtype=object)
    for i in range(n):
        for j in range(i+1,n):
            df[i][j] =  True if (aStart[j] - aStart[i]) > threshold and (aEnd[j] - aEnd[i]) > threshold and (bStart[j] - bStart[i]) > threshold and (bEnd[j] - bEnd[i]) > threshold else False
    return df


cpdef np.ndarray[object, ndim=2] makeBlocksTree(np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.ndarray bDir, np.ndarray aChr, np.ndarray bChr, np.ndarray index, np.int threshold, np.ndarray left, np.ndarray right):
    """Compute whether two alignments can be part of one translation block. For this:
        the alignments should not be separated by any inPlaceBlock on both ends and
        they should be collinear with respect to each other.
       
       Returns
       --------
       outOrderedBlocks: pandas DataFrame,
           Dataframe of type Object. Lower half is NA, upper half contains whether two
           alignments can be connected (True) or not (False).
    """

    assert(aStart.dtype==np.int and aEnd.dtype==np.int and bStart.dtype==np.int and bEnd.dtype==np.int and bDir.dtype==np.int and aChr.dtype==np.object and bChr.dtype==np.object and index.dtype==np.int and left.dtype==np.int and right.dtype==np.int)
    cdef Py_ssize_t i,j, n = len(aStart)
    assert(n == len(aEnd) == len(bStart) == len(bEnd) == len(index) == len(bDir) == len(aChr) == len(bChr) == len(left) == len(right))
    cdef np.ndarray[object, ndim=2] outOrderedBlocks =  np.array([[np.nan]*n]*n, dtype=object)
    cdef np.ndarray allRanges = np.array([range(left[i]+1, right[i]) for i in range(n)])

    for i in range(n):
        for j in range(i,n):
            #if len(np.intersect1d(range(left[i]+1,right[i]),range(left[j]+1,right[j]))) == 0:
            if bDir[i] != bDir[j]:
                sys.exit("ERROR: bDir not matching")
            elif not any([k in allRanges[i] for k in allRanges[j]]):
                    outOrderedBlocks[i][j] = False
            elif bDir[i] == bDir[j]:
                if (aStart[j] - aStart[i]) > threshold and (aEnd[j] - aEnd[i]) > threshold and (bStart[j] - bStart[i]) > threshold and (bEnd[j] - bEnd[i]) > threshold:
                    outOrderedBlocks[i][j] = True
                else:
                    outOrderedBlocks[i][j] = False
#            elif(aStart[j] - aStart[i]) > threshold and (aEnd[j] - aEnd[i]) > threshold and (bStart[j] - bEnd[i]) > threshold and (bEnd[j] - bStart[i]) > threshold:
#                outOrderedBlocks[i][j] = True
            else:
                outOrderedBlocks[i][j] = False
    return outOrderedBlocks



def getSynPath(blocks):
    cdef list synPath = []
    scores = [block.score for block in blocks]
    cdef int lastBlock = scores.index(max(scores))
    while blocks[lastBlock].bestParentID != -1:
        synPath.append(lastBlock)
        lastBlock = blocks[lastBlock].bestParentID        
    synPath.append(lastBlock)
    return(synPath[::-1])


def outSyn(cwdPath, threshold, prefix):
    ctxAnnoDict = {"duplication":"dupCtx",
                   "invDuplication":"invDupCtx",
                   "translocation":"TLCtx",
                   "invTranslocation":"invTLCtx"}
    reCoords =  pd.DataFrame()
        
    synData = []
    with open(cwdPath+prefix+"synOut.txt","r") as fin:
        for line in fin:
            line = line.strip().split("\t")
            if line[0] == "#":
                chromo = line[1]
                continue
            if len(line) == 4:
                synData.append(list(map(int,line[:4]))+[chromo,chromo])
            elif len(line) == 5:
                synData.append(list(map(int,line[:4]))+[chromo,chromo] + [line[4]])

    synData = pd.DataFrame(synData)
    if len(synData.columns) == 6:
        synData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
    else:
        synData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr","isinInv"]
    synData["class"] = "syn"
       
    for i in ["invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt","ctxOut.txt"]:    
        data = []
        with open(cwdPath+prefix+i,"r") as fin: 
            if i != "ctxOut.txt":
                for line in fin:
                    line = line.strip().split("\t")
                    if line[0] == "#":
                        data.append(list(map(int,getValues(line,[2,3,6,7]))) + [line[1],line[5]])
                data = pd.DataFrame(data, columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"], dtype=object)
                data["class"] = i.split("Out.txt")[0]
                if len(data)>0:
                    reCoords = reCoords.append(data)
            else:
                for line in fin:
                    line = line.strip().split("\t")
                    if line[0] == "#":
                        data.append(list(map(int,getValues(line,[2,3,6,7]))) + [line[1],line[5],ctxAnnoDict[line[8]]])
                data = pd.DataFrame(data, columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr","class"], dtype=object)
                if len(data)>0:
                    reCoords = reCoords.append(data)
                
    allBlocks = synData[["aStart","aEnd","bStart","bEnd","aChr","bChr","class"]].append(reCoords)
    allBlocks.index = range(allBlocks.shape[0])
    allBlocks.sort_values(["aChr","aStart","aEnd","bChr","bStart","bEnd"], inplace= True)
    synLocs = {np.where(allBlocks.index.values == i)[0][0]:i for i in range(synData.shape[0])}

    allBlocks.index = range(allBlocks.shape[0])

    aClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():    
        if len(currentCluster) == 0:
            if row["class"] != "syn":
                continue
            elif row["class"] == "syn":
                curChr = row["aChr"]
                currentCluster.append(index)
        elif row["class"] == "syn":
            if row["aChr"] == curChr:
                currentCluster.append(index)
            else:
                aClusters.append(currentCluster)
                currentCluster = [index]
                curChr = row["aChr"]
        
        elif row["class"] in ["TL", "inv","invTL","TLCtx","invTLCtx"]:
            aClusters.append(currentCluster)
            currentCluster = []
            curChr = ""
        else:
            if row["aEnd"] < allBlocks.loc[currentCluster[-1]]["aEnd"] + threshold:
                continue
            else:
                allClasses = allBlocks["class"][index:]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if max(row["aStart"],allBlocks.loc[currentCluster[-1]]["aEnd"]) > allBlocks.loc[nextSyn]["aStart"] - threshold:
                        continue
                    else:
                        aClusters.append(currentCluster)
                        currentCluster = []
                else:
                    aClusters.append(currentCluster)
                    currentCluster = []
    aClusters.append(currentCluster)
                    
    allBlocks.sort_values(["bChr","bStart","bEnd","bChr","aStart","aEnd"],inplace = True)
    bClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():
        if len(currentCluster) == 0:
            if row["class"] != "syn":
                continue
            elif row["class"] == "syn":
                curChr = row["bChr"]
                currentCluster.append(index)
        elif row["class"] == "syn":
            if row["aChr"] == curChr:
                currentCluster.append(index)
            else:
                bClusters.append(currentCluster)
                currentCluster = [index]
                curChr = row["bChr"]
        elif row["class"] in ["TL", "inv","invTL","TLCtx","invTLCtx"]:
            bClusters.append(currentCluster)
            currentCluster = []
            curChr = ""
        else:
            if row["bEnd"] < allBlocks.loc[currentCluster[-1]]["bEnd"] + threshold:
                continue
            else:
                allClasses = allBlocks["class"][index:]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if max(row["bStart"], allBlocks.loc[currentCluster[-1]]["bEnd"]) > allBlocks.loc[nextSyn]["bStart"] - threshold:
                        continue
                    else:
                        bClusters.append(currentCluster)
                        currentCluster = []
                else:
                    bClusters.append(currentCluster)
                    currentCluster = []
    bClusters.append(currentCluster)
    allBlocks.sort_values(["aChr","aStart","aEnd","bChr", "bStart","bEnd"],inplace = True)
   
    outClusters = []
    aIndex = 0 
    bIndex = 0
    currentCluster = []
    for i in unlist(aClusters):
        if i in aClusters[aIndex] and i in bClusters[bIndex]:
            currentCluster.append(i)
        else:
            if i not in aClusters[aIndex]:
                aIndex+=1
            if i not in bClusters[bIndex]:
                bIndex+=1
            outClusters.append(currentCluster)
            currentCluster = [i]
    outClusters.append(currentCluster)
    
    hasSynInInv = "isinInv" in synData.columns
    
    with open(cwdPath+prefix+"synOut.txt","w", encoding="utf-8") as fout:
        for i in outClusters:
            fout.write("\t".join(map(str,["#",allBlocks.at[i[0],"aChr"],allBlocks.at[i[0],"aStart"],allBlocks.at[i[-1],"aEnd"],"-",allBlocks.at[i[0],"aChr"],allBlocks.at[i[0],"bStart"],allBlocks.at[i[-1],"bEnd"]])) +"\n")
            for j in i:
                fout.write("\t".join(map(str,allBlocks.loc[j][0:4])))
                if hasSynInInv and synData.loc[synLocs[j]]["isinInv"] == "Syn_in_Inv":
                    fout.write("\tSyn_in_Inv\n")
                else:
                    fout.write("\n")   
    return None
    
        
def groupSyn(tempInvBlocks, dupData, invDupData, invTLData, TLData, threshold, synData, badSyn):
    
    synData = synData.drop(synData.index.values[badSyn])
    allBlocks = synData[["aStart","aEnd","bStart","bEnd"]].copy()
    allBlocks["class"] = "syn"
    
    tempInvBlocks = pd.DataFrame(tempInvBlocks,columns =["aStart","aEnd","bStart","bEnd"], dtype= object)
    tempInvBlocks["class"] = "inv"
    
    tempDupData = dupData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempDupData["class"] = "dup"
    
    tempInvDupData = invDupData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempInvDupData["class"] = "invDup"
    
    tempInvTLData = invTLData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempInvTLData["class"] = "invTL"
    
    tempTLData = TLData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempTLData["class"] = "TL"
    
    allBlocks = pd.concat([allBlocks,tempInvBlocks, tempInvDupData, tempInvTLData, tempTLData, tempDupData])
    allBlocks.index = range(allBlocks.shape[0])
    
    """
    Take data of all blocks and create groups of syntenic blocks from syntenic alignments
    """
    
    allBlocks.sort_values(["aStart","aEnd","bStart","bEnd"],inplace = True)
    
    aClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():        
        if len(currentCluster) == 0 and row["class"] != "syn":
            continue
        
        if row["class"] == "syn":
            currentCluster.append(index)
        elif row["class"] in ["TL", "inv","invTL"]:
            aClusters.append(currentCluster)
            currentCluster = []
        else:
            if row["aEnd"] < allBlocks.loc[currentCluster[-1]]["aEnd"] + threshold:
                continue
            else:
                allClasses = allBlocks["class"][index:]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if row["aStart"] > allBlocks.loc[nextSyn]["aStart"] - threshold:
                        continue
                    else:
                        aClusters.append(currentCluster)
                        currentCluster = []
                else:
                    aClusters.append(currentCluster)
                    currentCluster = []
    aClusters.append(currentCluster)
    
    allBlocks.sort_values(["bStart","bEnd","aStart","aEnd"],inplace = True)
    
    bClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():
        
        if len(currentCluster) == 0 and row["class"] != "syn":
            continue
        
        if row["class"] == "syn":
            currentCluster.append(index)
        elif row["class"] in ["TL", "inv","invTL"]:
            bClusters.append(currentCluster)
            currentCluster = []
        else:
            if row["bEnd"] < allBlocks.loc[currentCluster[-1]]["bEnd"] + threshold:
                continue
            else:
                allClasses = allBlocks["class"][index:]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if row["bStart"] > allBlocks.loc[nextSyn]["bStart"] - threshold:
                        continue
                    else:
                        bClusters.append(currentCluster)
                        currentCluster = []
                else:
                    bClusters.append(currentCluster)
                    currentCluster = []
    bClusters.append(currentCluster)
    allBlocks.sort_values(["aStart","aEnd","bStart","bEnd"],inplace = True)
    
    outClusters = []
    aIndex = 0 
    bIndex = 0
    currentCluster = []
    for i in range(synData.shape[0]):
        if i in aClusters[aIndex] and i in bClusters[bIndex]:
            currentCluster.append(i)
        else:
            if i not in aClusters[aIndex]:
                aIndex+=1
            if i not in bClusters[bIndex]:
                bIndex+=1
            outClusters.append(currentCluster)
            currentCluster = [i]
    outClusters.append(currentCluster)
    return (allBlocks, outClusters)

def mergeOutputFiles(uniChromo,path,prefix):
    def addData(fName,anno, chromo):
        fPath = open(path+prefix+chromo+"_"+anno+"Out.txt","r")
        for line in fPath.readlines():
            line = line.strip().split("\t")
            if line[0] == "#":
                fName.write("\t".join(unlist([line[0], chromo, line[1:4], chromo, line[4:]])) + "\n")
            else:
                fName.write("\t".join(line) + "\n")
        fPath.close()
        fileRemove(path+prefix+chromo+"_"+anno+"Out.txt")
                
    fSyn = open(path+prefix+"synOut.txt","w")
    fInv = open(path+prefix+"invOut.txt","w")
    fTL = open(path+prefix+"TLOut.txt","w")
    fInvTL = open(path+prefix+"invTLOut.txt","w")
    fDup = open(path+prefix+"dupOut.txt","w")
    fInvDup = open(path+prefix+"invDupOut.txt","w")
    
    files = [fSyn, fInv, fTL, fInvTL, fDup, fInvDup]
    classes = ["syn","inv","TL","invTL","dup","invDup"]
    
    for chromo in uniChromo:
        for i in range(len(classes)):
            addData(files[i], classes[i], chromo)
            
    for f in files:
        f.close()


class alignmentBlock:
    def __init__(self, id, children, data):
        self.id = id
        self.children = list(children)
        self.parents = []
        self.score = (data.aLen + data.bLen) * data.iden
        self.bestParentID = -1
    
    def addParent(self, parentID):
        self.parents.append(parentID)
        
    def bestParent(self,parentID,maxScore):
        self.bestParentID = parentID
        self.score = self.score + maxScore

