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

# cython: linetrace=True
np.random.seed(1)
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
    try:
        coords = pd.read_table(coordsfin, header = None)
    except pd.errors.ParserError:
        coords = pd.read_table(coordsfin, header = None, engine = "python")
    coords.columns = ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr"]

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


    uniChromo = list(np.unique(coords.aChr))
    logger.info('Analysing chromosomes: {}'.format(uniChromo))
    # Identify intra-chromosomal events (synteny, inversions, intra-trans, intra-dup) for each chromosome as a separate
    # process in parallel
    with Pool(processes = nCores) as pool:
        pool.map(partial(syri,threshold=threshold,coords=coords, cwdPath= cwdPath, bRT = bRT, prefix = prefix, tUC=tUC, tUP=tUP), uniChromo)

    # Merge output of all chromosomes
    mergeOutputFiles(uniChromo,cwdPath, prefix)

    #Identify cross-chromosomal events in all chromosomes simultaneously
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
    blocks = [alingmentBlock(i, np.where(df.iloc[i,] == True)[0], coordsData.iloc[i]) for i in range(nrow)]
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

    invertedCoordsOri, profitable, bestInvPath, invData, synInInv, badSyn = getInversions(coords,chromo, threshold, synData, synPath)
    

    ##########################################################
    #### Identify Translocation and duplications
    ##########################################################
    logger.info("Identifying translocation and duplication for chromosome " + chromo)

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

    # with open("members.txt","w") as fout:
    #     for i in allTransCluster:    #         fout.write(",".join(map(str, i.member)) + "\n")
    #     fout.write("\n\n")
    #     for i in allTransGenomeBGroups:
    #         fout.write(",".join(map(str, i.member)) + "\n")

    logger.debug("Translocations : making blocks data " + chromo +" " + str(datetime.now()))
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    # isSynOverlap = getOverlapWithSynBlocks(np.array(allTransBlocks.aStart), np.array(allTransBlocks.aEnd),
    #                               np.array(allTransBlocks.bStart), np.array(allTransBlocks.bEnd),
    #                               np.array(inPlaceBlocks.aStart), np.array(inPlaceBlocks.aEnd),
    #                               np.array(inPlaceBlocks.bStart), np.array(inPlaceBlocks.bEnd), 50,
    #                               allTransBlocks.shape[0], tUC, tUP)


    # allTransBlocks.to_csv(cwdPath+"allTransBlocks.txt", sep="\t", index = False)
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



#     allTransBlocksData = []
#     for i in range(allTransBlocks.shape[0]):
#         temp_trans_block = transBlock(allTransBlocks.iat[i, 0], \
#                                       allTransBlocks.iat[i, 1], \
#                                       allTransBlocks.iat[i, 2], \
#                                       allTransBlocks.iat[i, 3], \
#                                       allTransBlocks.iat[i, 4], \
#                                       allTransClusterIndices[i], \
#                                       i)
#         temp_trans_block.addTransGroupIndices(allTransGroupIndices[i])
#         temp_trans_block.checkOverlapWithSynBlocks(inPlaceBlocks, threshold)
#         temp_trans_block.checkoverlaps(allTransGenomeAGroups, allTransGenomeBGroups)
#         if (temp_trans_block.aUni and temp_trans_block.genomeAUni) or (
#                 temp_trans_block.bUni and temp_trans_block.genomeBUni):
#             temp_trans_block.setStatus(1)
#         allTransBlocksData.append(temp_trans_block)

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

    # Code for finding mutually exclusive blocks
    # for i in range(allTransBlocks.shape[0]):
    #     # if i%20000 == 0:
    #     #     logger.debug("transBlockAdded" + str(i))
    #     temp = allTransBlocksData[i]            ## create temporary trans blocks
    #     if not temp.aUni and not temp.bUni:
    #         allTransCluster[allTransClusterIndices[i]].remove(i)
    #     elif temp.status == 1:
    #         continue
    #     elif not temp.aUni:
    #         for j in temp.getoverlappingregions(groups = allTransGenomeBGroups,genome="b"):
    #             if allTransBlocksData[j].bStart - threshold < temp.bStart and allTransBlocksData[j].bEnd + threshold > temp.bEnd:
    #                 temp.addMEBlock(j)
    #     elif not temp.bUni:
    #         for j in temp.getoverlappingregions(groups = allTransGenomeAGroups,genome="a"):
    #             if allTransBlocksData[j].aStart - threshold < temp.aStart and allTransBlocksData[j].aEnd + threshold > temp.aEnd:
    #                 temp.addMEBlock(j)
    #     else:
    #         me_a = []       ## list of mutually exclusive regions on "a" genome
    #         for j in temp.getoverlappingregions(groups = allTransGenomeAGroups,genome="a"):
    #             if allTransBlocksData[j].aStart - threshold < temp.aStart and allTransBlocksData[j].aEnd + threshold > temp.aEnd:
    #                 me_a.append(j)
    #         me_b = []       ## list of mutually exclusive regions on "b" genome
    #         for j in temp.getoverlappingregions(groups = allTransGenomeBGroups,genome="b"):
    #             if allTransBlocksData[j].bStart - threshold < temp.bStart and allTransBlocksData[j].bEnd + threshold > temp.bEnd:
    #                 me_b.append(j)
    #         temp.setMEList(me_a, me_b)

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
    
def getBlocks(orderedBlocks, annoCoords, threshold, tUC, tUP):
    if len(orderedBlocks) == 0:
        return([])
    outOrderedBlocks = pd.DataFrame(makeBlocksTree_ctx(orderedBlocks.aStart.values, orderedBlocks.aEnd.values, orderedBlocks.bStart.values, orderedBlocks.bEnd.values, orderedBlocks.bDir.values, orderedBlocks.aChr.values, orderedBlocks.bChr.values, orderedBlocks.index.values, threshold))
    transBlocks = findOrderedTranslocations(outOrderedBlocks, orderedBlocks, annoCoords, threshold, tUC, tUP, ctx = True)
    return(transBlocks)
    
    
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

    # ctxBlocksData = deque()
    # for i in ctxTransBlocks.index.values:
    #     tempTransBlock = transBlock(ctxTransBlocks.at[i,"aStart"],\
    #                                 ctxTransBlocks.at[i,"aEnd"],\
    #                                 ctxTransBlocks.at[i,"bStart"],\
    #                                 ctxTransBlocks.at[i,"bEnd"],\
    #                                 ctxTransBlocks.at[i,"bDir"],\
    #                                 ctxClusterIndices[i],\
    #                                 i)
    #     # tempTransBlock.addTransGroupIndices(ctxGroupIndices[i])
    #     tempTransBlock.checkOverlapWithSynBlocks_A(annoCoords.loc[annoCoords.aChr == ctxTransBlocks.at[i,"aChr"]], threshold)
    #     tempTransBlock.checkOverlapWithSynBlocks_B(annoCoords.loc[annoCoords.bChr == ctxTransBlocks.at[i,"bChr"]], threshold)
    #     # tempTransBlock.addGenomeGroupMembers(ctxTransGenomeAGroups, ctxTransGenomeBGroups)
    #     # if (tempTransBlock.aUni and tempTransBlock.genomeAUni) or (tempTransBlock.bUni and tempTransBlock.genomeBUni):
    #     #     tempTransBlock.setStatus(1)
    #     ctxBlocksData.append(tempTransBlock)
    # ctxBlocksData = list(ctxBlocksData)

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

    # for i in range(len(ctxBlocksData)):
    #     tempTransBlock = ctxBlocksData[i]
    #     index = tempTransBlock.transBlocksID
    #     if not tempTransBlock.aUni and not tempTransBlock.bUni:
    #         ctxCluster[ctxClusterIndices[index]].remove(index)
    #     elif tempTransBlock.status == 1:
    #         continue
    #     elif not tempTransBlock.aUni:
    #         for j in tempTransBlock.getoverlappingregions(groups = ctxTransGenomeBGroups,genome="b"):
    #             if ctxBlocksData[j].bStart - threshold < tempTransBlock.bStart and ctxBlocksData[j].bEnd + threshold > tempTransBlock.bEnd:
    #                 tempTransBlock.addMEBlock(j)
    #     elif not tempTransBlock.bUni:
    #         for j in tempTransBlock.getoverlappingregions(groups = ctxTransGenomeAGroups,genome="a"):
    #             if ctxBlocksData[j].aStart - threshold < tempTransBlock.aStart and ctxBlocksData[j].aEnd + threshold > tempTransBlock.aEnd:
    #                 tempTransBlock.addMEBlock(j)
    #     else:
    #         ME_A = []
    #         for j in tempTransBlock.getoverlappingregions(groups = ctxTransGenomeAGroups,genome="a"):
    #             if ctxBlocksData[j].aStart - threshold < tempTransBlock.aStart and ctxBlocksData[j].aEnd + threshold > tempTransBlock.aEnd:
    #                 ME_A.append(j)
    #         ME_B = []
    #         for j in tempTransBlock.getoverlappingregions(groups = ctxTransGenomeBGroups,genome="b"):
    #             if ctxBlocksData[j].bStart - threshold < tempTransBlock.bStart and ctxBlocksData[j].bEnd + threshold > tempTransBlock.bEnd:
    #                 ME_B.append(j)
    #         tempTransBlock.setMEList(ME_A, ME_B)

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

cpdef apply_TS(np.ndarray aStart, np.ndarray aEnd, np.ndarray bStart, np.ndarray bEnd, np.int threshold):
    assert(aStart.dtype == np.int and aEnd.dtype == np.int and bStart.dtype == np.int and bEnd.dtype == np.int)
    cdef Py_ssize_t i, j,  n = len(aStart)
    assert(n == len(aEnd) == len(bStart) == len(bEnd))
    cdef np.ndarray[object, ndim =2 ] df =  np.array([[np.nan]*n]*n, dtype=object)
    for i in range(n):
        for j in range(i+1,n):
            df[i][j] =  True if (aStart[j] - aStart[i]) > threshold and (aEnd[j] - aEnd[i]) > threshold and (bStart[j] - bStart[i]) > threshold and (bEnd[j] - bEnd[i]) > threshold else False
    return df


def getSynPath(blocks):
    cdef list synPath = []
    scores = [block.score for block in blocks]
    cdef int lastBlock = scores.index(max(scores))
    while blocks[lastBlock].bestParentID != -1:
        synPath.append(lastBlock)
        lastBlock = blocks[lastBlock].bestParentID        
    synPath.append(lastBlock)
    return(synPath[::-1])

cpdef getInvBlocks(invTree, invertedCoordsOri):
    cdef int nrow, i, child
    nrow = invTree.shape[0]
    invBlocks = [alingmentBlock(i, np.where(invTree.iloc[i,] == True)[0], invertedCoordsOri.iloc[i]) for i in range(nrow)]
    
    for block in invBlocks:
        i = 0
        while(i < len(block.children)):
            block.children = list(set(block.children) - set(invBlocks[block.children[i]].children))
            i+=1
        block.children.sort()
        
        for child in block.children:
            invBlocks[child].addParent(block.id)
    return(invBlocks)

def getConnectivityGraph(blocksList):
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


cdef getAllLongestPaths(n,sNode, eNode, np.ndarray[np.int32_t, ndim =1] source, np.ndarray[np.int32_t, ndim =1] target, np.ndarray[np.float32_t, ndim=1] weight, by="weight"):
    """Uses Bellman-Ford Algorithm to find the shortest path from node "sNode" in the 
    directed acyclic graph "graph" to all nodes in the list "eNode". Edges weighed 
    are negative, so shortest path from sNode to eNode corresponds to the longest path.
       
        Parameters
        ----------
        graph: directeed igraph Graph(),
            Directed acyclic graph containing all the nodes and edges in the graph.
           
        sNode: int, 
            index of the start node in the graph.
       
        eNode: int list,
            list of all end nodes. longest path from start node to end nodes will be
            calculated
        
        by: igraph edge weight
        
        Returns
        -------
        list of len(eNodes) longest paths from sNodes to eNodes
        
    """
    pathList = deque()
    cdef:
        cdef Py_ssize_t i, j, current
        np.ndarray[np.int32_t, ndim =1] pred = np.array([-1]*n, dtype = np.int32)
        np.ndarray[np.float32_t, ndim=1] dist = np.array([np.float32('inf')]*n, dtype = np.float32)

    dist[sNode] = 0
    changes = True
    for i in range(n-1):
        if not changes:
            break
        changes = False
        for j in range(len(source)):
            if dist[source[j]] + weight[j] < dist[target[j]]:
                changes = True
                dist[target[j]] = dist[source[j]] + weight[j]
                pred[target[j]] = source[j]

    for j in range(len(source)):
        if dist[source[j]] + weight[j] < dist[target[j]]:
            sys.exit("Negative weight cycle identified")


    # Store paths for enodes which have already been found
    nodepath = {}
    for i in eNode:
        current = i
        if dist[current] != float("inf"):
            path = deque()
            while current!=sNode:
                if current in nodepath.keys():
                    path.extend(nodepath[current].copy())
                    break
                path.append(current)
                current = pred[current]
            nodepath[i] = path.copy()
            path.append(sNode)
            pathList.append(np.array(path, dtype="int32")[::-1])
    return(pathList)


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
    print(short)
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
               
def getRedundantIndex(inPlaceBlocks, outPlaceBlocks, threshold):
    nrow = outPlaceBlocks.shape[0]
    redundant = []
    for i in range(nrow):
        aRed = set(np.where(outPlaceBlocks.iat[i,0] >= (inPlaceBlocks.aStart-threshold))[0]).intersection(
            np.where(outPlaceBlocks.iat[i,1] <= (inPlaceBlocks.aEnd+threshold))[0])
        if outPlaceBlocks.iat[i,8] == 1: 
            bRed = set(np.where(outPlaceBlocks.iat[i,2] >= (inPlaceBlocks.bStart-threshold))[0]).intersection(
                np.where(outPlaceBlocks.iat[i,3] <= (inPlaceBlocks.bEnd + threshold))[0])
        else:
            bRed = set(np.where(outPlaceBlocks.iat[i,3] >= (inPlaceBlocks.bStart-threshold))[0]).intersection(
                np.where(outPlaceBlocks.iat[i,2] <= (inPlaceBlocks.bEnd+threshold))[0])
        if len(aRed) > 0 and len(bRed) > 0:
            redundant.append(i)
    return(redundant)


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
        

def findOverlappingSynBlocks(inPlaceBlocks, aStart, aEnd, bStart, bEnd):
#    aBlocks = list(np.intersect1d(np.where(inPlaceBlocks.aStart.values < aEnd)[0],\
#                                      np.where(inPlaceBlocks.aEnd.values > aStart)[0]))
#    bBlocks = list(np.intersect1d(np.where(inPlaceBlocks.bStart.values < bEnd)[0],\
#                                      np.where(inPlaceBlocks.bEnd.values > bStart)[0]))
    aBlocks = list(np.where((inPlaceBlocks.aStart.values < aEnd) & (inPlaceBlocks.aEnd.values > aStart) == True)[0])
    bBlocks = list(np.where((inPlaceBlocks.bStart.values < bEnd) & (inPlaceBlocks.bEnd.values > bStart) == True)[0])
    return(aBlocks, bBlocks)
    


#
#%%cython
#import numpy as np
#cimport numpy as np
#import sys
#
#
#


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


# noinspection PyUnreachableCode
def findOrderedTranslocations(outOrderedBlocks, orderedBlocks, inPlaceBlocks, threshold, tUC, tUP, ctx = False):
    logger = logging.getLogger("findOrderedTranslocations")
    if not isinstance(ctx, bool):
        print("CTX status must be a boolean")
        sys.exit()
    def makeBlocksList(blocksTree, blocksData):
        _blocksList = [alingmentBlock(_i, np.where(blocksTree.iloc[_i] == True)[0],blocksData.iloc[_i]) for _i in range(blocksTree.shape[0])]
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
        
    #%%            
def getTransOverlapGroups(transBlocks, orderedBlocks, threshold):
    transBlocksData = []
    for i in transBlocks:
        aStart = orderedBlocks.iat[i[0],0]
        aEnd = orderedBlocks.iat[i[-1],1]
        if orderedBlocks.iat[i[0],8] == 1:
            bStart = orderedBlocks.iat[i[0],2]
            bEnd = orderedBlocks.iat[i[-1],3]
        else:
            bStart = orderedBlocks.iat[i[0],3]
            bEnd = orderedBlocks.iat[i[-1],2]
        transBlocksData.append([aStart, aEnd, bStart, bEnd])
    transBlocksTable = pd.DataFrame(transBlocksData)
    transBlocksTable.columns =  ["aStart","aEnd","bStart","bEnd"]
    transBlocksTable.sort_values(["aStart","aEnd"], inplace = True)
    genomeAGroups = makeTransGroupList(transBlocksTable, "aStart","aEnd", threshold)
    transBlocksTable.sort_values(["bStart","bEnd"], inplace = True)
    genomeBGroups = makeTransGroupList(transBlocksTable, "bStart","bEnd", threshold)
    return genomeAGroups, genomeBGroups


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
#%%

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
             
#
#def getTransCluster(transGroupIndices, transGenomeAGroups, transGenomeBGroups):
#    assert(list(transGroupIndices.keys()) == list(range(len(transGroupIndices))))
#    nodeStack = np.zeros(len(transGroupIndices), dtype='uint8')
#    visitedTransBlock = np.zeros(len(transGroupIndices), dtype='uint8')
#    transCluster = []
#    
#    j = 0
#    for key,value in transGroupIndices.items():
#        if visitedTransBlock[key] == 0:
#            visitedTransBlock[key]=1
#            newGroup = [key]
#            node1 = value[0]
#            node2 = value[1]
#            nodeStack[transGenomeAGroups[node1].member] = 1
#            nodeStack[transGenomeBGroups[node2].member] = 1
#            nodeStack[np.nonzero(visitedTransBlock)[0]] = 0
#            while len(np.nonzero(nodeStack)[0]) > 0:
#                newKey = np.where(nodeStack == 1)[0][0]
#                if visitedTransBlock[newKey]== 0:
#                    visitedTransBlock[newKey] = 1
#                    newGroup.append(newKey)
#                    nodeStack[transGenomeAGroups[transGroupIndices[newKey][0]].member] = 1 
#                    nodeStack[transGenomeBGroups[transGroupIndices[newKey][1]].member] = 1
#                    nodeStack[np.nonzero(visitedTransBlock)[0]] = 0
#            newGroup.sort()
#            transCluster.append(list(newGroup))
#    return(transCluster)


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

#
# %%cython
# import numpy as np
# cimport numpy as np
cpdef getOverlapWithSynBlocks(np.ndarray[np.int_t, ndim=1] start, np.ndarray[np.int_t, ndim=1] end, np.ndarray chrom, np.ndarray[np.int_t, ndim=1] in_start, np.ndarray[np.int_t, ndim=1] in_end, np.ndarray in_chrom, np.int_t threshold, np.int_t count, np.int_t tUC, np.float_t tUP):

    assert(len(start) == len(end) == len(chrom) ==count)
    assert(len(in_start) == len(in_end) == len(in_chrom))

    cdef Py_ssize_t i, j, n = len(in_start)
    cdef np.int_t blockuni, s, e
    cdef np.ndarray[np.npy_bool, ndim = 1, cast=True] uni = np.zeros(count, dtype="bool")
    cdef np.ndarray[np.int_t, ndim=1] blocks

    for i in range(count):
        blocks = np.zeros(n, dtype="int")
        blockuni = 0
        s = start[i]
        e = end[i]
        c = chrom[i]
        for j in range(n):
            if (in_start[j] < end[i]) and (in_end[j] > start[i]) and (in_chrom[j] == c):
                blocks[j] = 1


        for j in range(len(blocks)):
            if blocks[j] == 1:
                if (in_start[j] - s < threshold) and (e - in_end[j] < threshold):
                    s = e
                    break
                elif in_start[j] < s and in_end[j] < e:
                    s = in_end[j]
                else:
                    blockuni+= in_start[j]-s
                    if in_end[j] < e:
                        s = in_end[j]
                    else:
                        s = e
                        break
        blockuni+= (e-s)

        if (blockuni > tUC) or (blockuni > tUP*(end[i]-start[i])):
            uni[i]=True
    return uni



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

                
def getscafDict(scaffolds, data, scafSize):
    """
    Select the reference chromosome to which the scaffold align best.
    Four conditions are checked (progressively) to select the best aligned region. 
    1) If the scaffold align to only one refChrom then that chromosome is selected
    2) If all refChrom have similar distribution for alignment length the query scaffold,
    then the refChrom which have total longest alignment is selected
    3) If refChrom_A_alignment_length - refGenome_B_alignment_length > a given threshold,
    then refGenome A is selected
    4) refGenome which has the average alignment length is selected
    """
    scafChrDict = {}
    scafCountDict = {}
    scafSizeDict = {}
    percentCutoff = 0.1
    for i in scaffolds:
        scafCountDict[i] = Counter(data.iloc[np.where(data[10] == i)[0], 9])
        scafData = data.loc[data[10] == i]
        uniChr = np.unique(scafData[9])
        bpSizeDict = {}
        alignSizes = {}
        for j in uniChr:
            indices = np.array(scafData.iloc[np.where(scafData[9] == j)[0], [0,1]])
            bpSizeDict[j] = count_uniq_elems(indices)    
            alignSizes[j] = list(scafData.iloc[np.where(scafData[9]==j)[0],4])
        scafSizeDict[i] = bpSizeDict   
        if len(uniChr) == 1:
            scafChrDict[i] = uniChr[0]
        else: 
            try:
                krusStat = kruskal(*alignSizes.values())[1]
            except ValueError as e:
                print(e, " for scaffold", i,". Ignoring this scaffold")
                krusStat = -1
            except e:
                print("Exception:",e,"\n incorrect values for kruskal. Ignoring scaffold",i)
                krusStat = -1
                
            if krusStat == -1:
                continue
            elif krusStat > 0.05:
                scafChrDict[i] = max(scafSizeDict[i].items(), key=lambda x: x[1])[0]
            else:
                size = scafSize[i]
                percentCutoffSize = percentCutoff*size
                bpSizeValues = sorted(bpSizeDict.values())
                identified = 0
#                for j in range(len(bpSizeValues)-1,0,-1):
                l = len(bpSizeValues)-1
                if (bpSizeValues[l] - bpSizeValues[l-1]) > percentCutoffSize:
                    scafChrDict[i] = list(bpSizeDict.keys())[list(bpSizeDict.values()).index(bpSizeValues[l])]
                    identified = 1
#                        break
                if identified == 0:
                    meanAlignSize = {}
                    for j in uniChr:
                        meanAlignSize[j] = scafSizeDict[i][j]/scafCountDict[i][j]
                    scafChrDict[i] = max(meanAlignSize.items(), key=lambda x:x[1])[0]
    return (scafCountDict, scafSizeDict, scafChrDict)

def orderFromMummerplot(filePath):
    gpData = open(filePath,"r").readlines()
    orderedScafID = []
    inverted = []
    started = 0
    for line in gpData:
        line = line.strip()
        if started and line ==  ')':
            break
        if started:
            ID = line.split(" ")[0].replace('"',"")
            if len(ID) > 0:
                if "*" == ID[0]:
                    ID = ID.replace("*","")
                    inverted.append(ID)
                orderedScafID.append(ID)
            continue
        if "set ytics" in line:
            started = 1
            continue
    return orderedScafID,inverted

def invertAlignmentDirection(tempData, scafSize):
    size = scafSize[np.unique(tempData[10])[0]]
    a = size  - tempData[2] + 1
    b = size - tempData[3] + 1
    c = -1*tempData[8]
    newTempData = tempData.copy()
    newTempData[2] = a
    newTempData[3] = b
    newTempData[8] = c
    return(newTempData)
    
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

def outSyn(cwdPath, threshold, prefix):
#    reCoords = pd.DataFrame(columns=["aStart","aEnd","bStart","bEnd","aChr","bChr"])
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
#    fin.close()
    
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
            
#%%

class alingmentBlock:
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
        
class transBlock:
#     def __init__(self,aStart, aEnd, bStart, bEnd, Dir, transClusterIndex, i):
#         self.aStart = aStart
#         self.aEnd = aEnd
#         self.bStart = bStart
#         self.bEnd = bEnd
#         self.dir = Dir
# #        self.orderedBlocksIndex = orderedBlocksIndex
#         self.transBlocksID = i
#         self.transClusterIndex = transClusterIndex
#         self.status = 0
#         self.overlappingInPlaceBlocks = []
#         self.aUni = False
#         self.bUni = False

    def __init__(self, i):
        self.aStart = None
        self.aEnd = None
        self.bStart = None
        self.bEnd = None
        self.dir = None
        #        self.orderedBlocksIndex = orderedBlocksIndex
        self.transBlocksID = i
        self.transClusterIndex = None
        self.status = 0
        self.overlappingInPlaceBlocks = []
        self.aUni = False
        self.bUni = False
        self.transGroupIndices = None
        self.genomeAUni = False
        self.genomeBUni = False


    # def addGenomeGroupMembers(self,ctxTransGenomeAGroups, ctxTransGenomeBGroups):
    #     aMem = ctxTransGenomeAGroups[self.transGroupIndices[0]].member.copy()
    #     aMem.remove(self.transBlocksID)
    #     bMem = ctxTransGenomeBGroups[self.transGroupIndices[1]].member.copy()
    #     bMem.remove(self.transBlocksID)
    #     self.genomeAMembers = np.array(aMem, dtype = 'int32')
    #     self.genomeBMembers = np.array(bMem, dtype = 'int32')
    #     self.genomeAUni = True if len(self.genomeAMembers) == 0 else False
    #     self.genomeBUni = True if len(self.genomeBMembers) == 0 else False
    #

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
#        bBlocks = list(np.intersect1d(np.where(inPlaceBlocks.bStart.values < self.bEnd)[0],\
#                                      np.where(inPlaceBlocks.bEnd.values > self.bStart)[0]))
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

    # def addMEBlock(self, blockID):
    #     """List of Blocks which prohibit the entry of current block in the
    #     optimal solution"""
    #
    #     try:
    #         self.meTo.extend(blockID) if type(blockID) == list else self.meTo.append(blockID)
    #     except AttributeError:
    #         self.meTo = blockID if type(blockID) == list else [blockID]


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

# cpdef getmeblocks(np.ndarray[np.int_t, ndim=1] aStart, np.ndarray[np.int_t, ndim=1] aEnd, np.ndarray[np.int_t, ndim=1] bStart, np.ndarray[np.int_t, ndim=1] bEnd, np.int_t threshold, np.int_t count, np.ndarray[np.int_t, ndim=1] aUni, np.ndarray[np.int_t, ndim=1] bUni, np.ndarray[np.int_t, ndim=1] status, np.ndarray[np.int_t, ndim=1] aIndex, np.ndarray[np.int_t, ndim=1] bIndex, aGroups, bGroups):
#     # Function take the coordinates and cluster information of all translocated blocks and identifies mutually exclusive
#     #  blocks by comparing the coordinates of each block to the coordinates of the member blocks in its cluster
#     cdef np.ndarray[np.int_t, ndim=1] members
#     cdef np.int_t i, j
#     rem = deque()
#     meblock = {}            ## for blocks which are overlapping with inplace blocks
#     melist = {}             ## for blocks which are not overlapping with inplace blocks
#
#     assert(len(aStart) == len(aEnd) == len(bStart)== len(bEnd))
#     assert(count == len(aUni)== len(bUni)== len(status)== len(aIndex)== len(bIndex))
#
#     for i in range(count):
#         if i%50000 == 0:
#             print("Number of mutually exclusive blocks identified", i)
#         if not aUni[i] and not bUni[i]:
#             rem.append(i)
#         elif status[i] == 1:
#             continue
#         elif not aUni[i]:
#             meb = deque()               ## vector of mutually exclusive block
#             members = bGroups[bIndex[i]]
#             for j in members:
#                 if j!=i:
#                     if bStart[j] - threshold < bStart[i] and bEnd[j] + threshold > bEnd[i]:
#                         meb.append(j)
#             meblock[i] = meb
#         elif not bUni[i]:
#             meb = deque()                ## vector of mutually exclusive block
#             members = aGroups[aIndex[i]]
#             for j in members:
#                 if j!=i:
#                     if aStart[j] - threshold < aStart[i] and aEnd[j]+threshold > aEnd[i]:
#                         meb.append(j)
#             meblock[i] = meb
#         else:
#             meb_a = deque()             ## vector of mutually exclusive block on A genome
#             meb_b = deque()             ## vector of mutually exclusive block on B genome
#             members = aGroups[aIndex[i]]
#             for j in members:
#                 if j!=i:
#                     if aStart[j] - threshold < aStart[i] and aEnd[j] + threshold > aEnd[i]:
#                         meb_a.append(j)
#
#             members = bGroups[bIndex[i]]
#             for j in members:
#                 if j != i:
#                     if bStart[j] - threshold < bStart[i] and bEnd[j] + threshold > bEnd[i]:
#                         meb_b.append(j)
#             melist[i] = (meb_a, meb_b)
#     return rem, meblock, melist

# %%cython
# import numpy as np
# cimport numpy as np
# from collections import deque
cpdef getmeblocks(np.ndarray[np.int_t, ndim=1] aStart, np.ndarray[np.int_t, ndim=1] aEnd, np.ndarray[np.int_t, ndim=1] bStart, np.ndarray[np.int_t, ndim=1] bEnd, np.int_t threshold, np.int_t count, np.ndarray[np.int_t, ndim=1] aUni, np.ndarray[np.int_t, ndim=1] bUni, np.ndarray[np.int_t, ndim=1] status, np.ndarray[np.int_t, ndim=1] aIndex, np.ndarray[np.int_t, ndim=1] bIndex, aGroups, bGroups):
    # Function take the coordinates and cluster information of all translocated blocks and identifies mutually exclusive
    #  blocks by comparing the coordinates of each block to the coordinates of the member blocks in its cluster
    logger = logging.getLogger("getmeblocks")
    cdef np.ndarray[np.int_t, ndim=1] members, temp
    cdef np.ndarray[np.npy_bool, ndim=1, cast=True] meb, meb_a, meb_b, rem = np.zeros(count, dtype="bool")

    cdef np.int_t i, j, index
    # rem = deque()
    meblock = {}            ## for blocks which are overlapping with inplace blocks
    melist = {}             ## for blocks which are not overlapping with inplace blocks

    assert(len(aStart) == len(aEnd) == len(bStart)== len(bEnd))
    assert(count == len(aUni)== len(bUni)== len(status)== len(aIndex)== len(bIndex))

    for i in range(count):
        if i%50000 == 0:
            logger.debug("Number of mutually exclusive blocks identified " + str(i))
        if not aUni[i] and not bUni[i]:
            rem[i] = True
        elif status[i] == 1:
            continue
        elif not aUni[i]:
            members = bGroups[bIndex[i]]
            meb = np.zeros(len(members), dtype="bool")              ## vector of mutually exclusive block
            for index in range(len(members)):
                j = members[index]
                if j!=i:
                    if bStart[j] - threshold < bStart[i] and bEnd[j] + threshold > bEnd[i]:
                        meb[index]=True
            meblock[i] = np.array(members[meb], dtype="uint32")
        elif not bUni[i]:
            members = aGroups[aIndex[i]]
            meb = np.zeros(len(members), dtype="bool")               ## vector of mutually exclusive block
            for index in range(len(members)):
                j = members[index]
                if j!=i:
                    if aStart[j] - threshold < aStart[i] and aEnd[j]+threshold > aEnd[i]:
                        meb[index] = True
            meblock[i] = np.array(members[meb], dtype="uint32")
        else:
            members = aGroups[aIndex[i]]
            meb_a = np.zeros(len(members), dtype="bool")             ## vector of mutually exclusive block on A genome
            for index in range(len(members)):
                j = members[index]
                if j!=i:
                    if aStart[j] - threshold < aStart[i] and aEnd[j] + threshold > aEnd[i]:
                        meb_a[index] = True
            temp = members[meb_a]

            members = bGroups[bIndex[i]]
            meb_b = np.zeros(len(members), dtype="bool")             ## vector of mutually exclusive block on B genome
            for index in range(len(members)):
                j = members[index]
                if j != i:
                    if bStart[j] - threshold < bStart[i] and bEnd[j] + threshold > bEnd[i]:
                        meb_b[index] = True
            melist[i] = (np.array(temp, dtype="uint32"), np.array(members[meb_b], dtype="uint32"))
    return rem, meblock, melist

#%%
#################################################################
### SV identification functions
#################################################################

def readSRData(cwdPath, prefix, dup = False):
    if not isinstance(dup, bool):
        sys.exit("need boolean")
    if dup:
        fin = ["synOut.txt","invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt", "ctxOut.txt"]
    else:
        fin = ["synOut.txt","invOut.txt", "TLOut.txt", "invTLOut.txt","ctxOut.txt"]
        
    annoCoords = pd.DataFrame()
    for fileType in fin[:-1]:
        try:
            fileData = pd.read_table(cwdPath+prefix+fileType, header=None, dtype = object)
        except pd.errors.ParserError as _e:
            fileData = pd.read_table(cwdPath+prefix+fileType, header=None, dtype = object, engine ="python")
        except pd.io.common.EmptyDataError:
            print(fileType, " is empty. Skipping analysing it.")
            continue
        except Exception as _e:
            print("ERROR: while trying to read ", fileType, _e)
            continue
            
        annoIndices = np.where(fileData[0] == "#")[0]
        annoIndices = np.append(annoIndices,len(fileData))
        repCount = annoIndices[1:] - annoIndices[:-1] - 1
        
        annoData = fileData.loc[fileData[0] == "#"].copy()
        coordsData = fileData.loc[fileData[0] !="#"].copy()
        coordsData = coordsData[[0,1,2,3]].astype(dtype = "int64")
        
        reps = []
        for i in annoData[1].unique():
            reps.extend(list(range(len(np.where(annoData[1] == i)[0]))))
        reps = np.repeat(reps, repCount)
        
        coordsData["group"] = reps
        coordsData["aChr"] = list(np.repeat(annoData[1], repCount))
        coordsData["bChr"] = list(np.repeat(annoData[5], repCount))
        coordsData["state"] = fileType.split("Out.txt")[0]
        annoCoords = annoCoords.append(coordsData.copy())
                                   
    try:
        pass
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header = None, dtype = object)
    except pd.errors.ParserError as e:
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header=None, dtype = object, engine ="python")
    except pd.io.common.EmptyDataError:
        print("ctxOut.txt is empty. Skipping analysing it.")
    except Exception as e:
        print("ERROR: while trying to read ", fileType, "Out.txt", e)
    
    annoIndices = np.where(fileData[0] =="#")[0]
    states = list(fileData[8].loc[annoIndices])
    coordsData = fileData.loc[fileData[0] =="#"].copy()
    coordsData1 = fileData.loc[fileData[0] !="#", [0,1,2,3]].copy().astype(dtype="int")
    annoIndices = np.append(annoIndices,len(fileData))
    repCount = annoIndices[1:] - annoIndices[:-1] - 1    
    reps = np.repeat(range(len(annoIndices)-1), repCount)
    stateReps = np.repeat(states, repCount)
    coordsData1["aChr"] = np.repeat(coordsData[1], repCount).tolist()
    coordsData1["bChr"] = np.repeat(coordsData[5], repCount).tolist()
    coordsData1["group"] = reps
    coordsData1["state"] = stateReps
    coordsData1 = coordsData1[[0,1,2,3,"group","aChr","bChr","state"]]
    coordsData1.loc[coordsData1.state == "translocation","state"] = "ctx"
    coordsData1.loc[coordsData1.state == "invTranslocation","state"] = "invCtx"
    coordsData1.loc[coordsData1.state == "duplication","state"] = "ctxDup"
    coordsData1.loc[coordsData1.state == "invDuplication","state"] = "ctxInvDup"
    if not dup:
        coordsData1 = coordsData1.loc[coordsData1["state"].isin(["ctx","invCtx"])]
    annoCoords = annoCoords.append(coordsData1)
    annoCoords.columns = ["aStart","aEnd","bStart","bEnd","group","aChr","bChr","state"]
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))
    return annoCoords

#
# def getSV(cwdPath, allAlignments, prefix, offset):
#     fout = open(cwdPath+prefix+"sv.txt","w")
#     allAlignments["id"] = allAlignments.group.astype("str") + allAlignments.aChr + allAlignments.bChr + allAlignments.state
#     allBlocks = pd.unique(allAlignments.id)
#
#     for i in allBlocks:
#         blocksAlign = allAlignments.loc[allAlignments.id == i].copy()
#         ordered = 1 if "inv" not in blocksAlign.state.iloc[0] else 0
#         for j in range(len(blocksAlign) - 1):
#             m = blocksAlign.iat[j+1,0] - blocksAlign.iat[j,1] - 1
#             if ordered:
#                 n = blocksAlign.iat[j+1,2] - blocksAlign.iat[j,3] - 1
#             else:
#                 n = blocksAlign.iat[j,3] - blocksAlign.iat[j+1,2] - 1
#
#
#             ## No overlap of reference genome
#             if m == 0:
#
#                 ## No overlap in query genome
#                 if n == 0:
#                     continue
#
#                 elif n > 0:
#                     if ordered:
#                         fout.write("\t".join(["InDel",
#                                               str(blocksAlign.iat[j,1]),
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,3] + 1),
#                                               str(blocksAlign.iat[j+1, 2] - 1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                     else:
#                         fout.write("\t".join(["InDel",
#                                               str(blocksAlign.iat[j,1]),
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,3] - 1),
#                                               str(blocksAlign.iat[j+1, 2] + 1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#
#                 elif n < 0:
#                     if ordered:
#                         j_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
#                         j1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])
#                         sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                         fout.write("\t".join(["CNV",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                     else:
#                         j_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
#                         j1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
#                         sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                         fout.write("\t".join(["CNV",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#
#             elif m == 1:
#
#                 if n == 0:
#                     if ordered:
#                         fout.write("\t".join(["InDel",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                     else:
#                         fout.write("\t".join(["InDel",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]),
#                                               str(blocksAlign.iat[j+1, 2]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#
#                 elif n==1:
#                     if ordered:
#                         fout.write("\t".join(["SNP",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]+1),
#                                               str(blocksAlign.iat[j+1,2]-1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                     else:
#                         fout.write("\t".join(["SNP",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]-1),
#                                               str(blocksAlign.iat[j+1,2]+1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                 elif n>1:
#                     if ordered:
#                         fout.write("\t".join(["InDel+SNP",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]+1),
#                                               str(blocksAlign.iat[j+1,2]-1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "SNP:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
#                     else:
#                         fout.write("\t".join(["InDel+SNP",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]-1),
#                                               str(blocksAlign.iat[j+1,2]+1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "SNP:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
#                 elif n<0:
#                     if ordered:
#                         j_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
#                         j1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])
#                         sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                         fout.write("\t".join(["CNV+InDel",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
#                     else:
#                         j_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
#                         j1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
#                         sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                         fout.write("\t".join(["CNV+InDel",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
#
#             elif m>1:
#
#                 if n==0:
#                     if ordered:
#                         fout.write("\t".join(["InDel",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                     else:
#                         fout.write("\t".join(["InDel",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                 elif n==1:
#                     if ordered:
#                         fout.write("\t".join(["InDel+SNP",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]+1),
#                                               str(blocksAlign.iat[j+1,2]-1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "SNP:Q:"+str(blocksAlign.iat[j,3]+1)+"-"+str(blocksAlign.iat[j+1,2]-1)]) + "\n")
#                     else:
#                         fout.write("\t".join(["InDel+SNP",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]-1),
#                                               str(blocksAlign.iat[j+1,2]+1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "SNP:Q:"+str(blocksAlign.iat[j,3]-1)+"-"+str(blocksAlign.iat[j+1,2]+1)]) + "\n")
#                 elif n>1:
#                     if ordered:
#                         fout.write("\t".join(["HDR",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]+1),
#                                               str(blocksAlign.iat[j+1,2]-1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) +"\n")
#                     else:
#                         fout.write("\t".join(["HDR",
#                                               str(blocksAlign.iat[j,1]+1),
#                                               str(blocksAlign.iat[j+1,0]-1),
#                                               str(blocksAlign.iat[j,3]-1),
#                                               str(blocksAlign.iat[j+1,2]+1),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) +"\n")
#                 elif n<0:
#                     if ordered:
#                         j_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
#                         j1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])
#                         sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                         fout.write("\t".join(["CNV+InDel",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
#                     else:
#                         j_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
#                         j1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
#                         sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                         fout.write("\t".join(["CNV+InDel",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
#
#             elif m<0:
#
#                 j_prop = abs(m) / (blocksAlign.iat[j,1] - blocksAlign.iat[j,0])
#                 j1_prop = abs(m) / (blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])
#
#                 if n==0:
#                     if ordered:
#                         sCoord = round(blocksAlign.iat[j,3] - j_prop*(blocksAlign.iat[j,3] - blocksAlign.iat[j,2])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,2] + j1_prop*(blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])).astype(int)
#                         fout.write("\t".join(["CNV",
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,1]),
#                                               str(sCoord),
#                                               str(eCoord),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#                     else:
#                         sCoord = round(blocksAlign.iat[j,3] + j_prop*(blocksAlign.iat[j,2] - blocksAlign.iat[j,3])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,2] - j1_prop*(blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])).astype(int)
#                         fout.write("\t".join(["CNV",
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,1]),
#                                               str(sCoord),
#                                               str(eCoord),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6]]) + "\n")
#
#                 if n>0:
#                     if ordered:
#                         sCoord = round(blocksAlign.iat[j,3] - j_prop*(blocksAlign.iat[j,3] - blocksAlign.iat[j,2])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,2] + j1_prop*(blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])).astype(int)
#                         fout.write("\t".join(["CNV+InDel",
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,1]),
#                                               str(sCoord),
#                                               str(eCoord),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "InDel:Q:"+str(blocksAlign.iat[j,3]+1)+"-"+str(blocksAlign.iat[j+1,2]-1)]) + "\n")
#                     else:
#                         sCoord = round(blocksAlign.iat[j,3] + j_prop*(blocksAlign.iat[j,2] - blocksAlign.iat[j,3])).astype(int)
#                         eCoord = round(blocksAlign.iat[j+1,2] - j1_prop*(blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])).astype(int)
#                         fout.write("\t".join(["CNV+InDel",
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,1]),
#                                               str(sCoord),
#                                               str(eCoord),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "InDel:Q:"+str(blocksAlign.iat[j,3]-1)+"-"+str(blocksAlign.iat[j+1,2]+1)]) + "\n")
#
#                 if n<0:
#                     maxOverlap = max(abs(m),abs(n))
#                     if abs(m-n) < 0.1*maxOverlap: ## no SV if the overlap on both genomes is of similar size
#                         continue
#
#                     if abs(m) > abs(n):
#                         if ordered:
#                             sCoord = round(blocksAlign.iat[j,3] - j_prop*(blocksAlign.iat[j,3] - blocksAlign.iat[j,2])).astype(int)
#                             eCoord = round(blocksAlign.iat[j+1,2] + j1_prop*(blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])).astype(int)
#                             fout.write("\t".join(["CNV+Tandem",
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,1]),
#                                               str(sCoord),
#                                               str(eCoord),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "Tandem:Q:"+str(blocksAlign.iat[j,3]+1)+"-"+str(eCoord)]) + "\n")
#
#                         else:
#                             sCoord = round(blocksAlign.iat[j,3] + j_prop*(blocksAlign.iat[j,2] -blocksAlign.iat[j,3])).astype(int)
#                             eCoord = round(blocksAlign.iat[j+1,2] - j1_prop*(blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,2])).astype(int)
#                             fout.write("\t".join(["CNV+Tandem",
#                                               str(blocksAlign.iat[j+1,0]),
#                                               str(blocksAlign.iat[j,1]),
#                                               str(sCoord),
#                                               str(eCoord),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "Tandem:Q:"+str(blocksAlign.iat[j,3]-1)+"-"+str(eCoord)]) + "\n")
#                     else:
#                         if ordered:
#                             k_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
#                             k1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j,2])
#                             sCoord = round(blocksAlign.iat[j,1] - k_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                             eCoord = round(blocksAlign.iat[j+1,0] + k1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                             fout.write("\t".join(["CNV+Tandem",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "Tandem:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(eCoord)]) + "\n")
#                         else:
#                             k_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
#                             k1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
#                             sCoord = round(blocksAlign.iat[j,1] - k_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
#                             eCoord = round(blocksAlign.iat[j+1,0] + k1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
#                             fout.write("\t".join(["CNV+Tandem",
#                                               str(sCoord),
#                                               str(eCoord),
#                                               str(blocksAlign.iat[j+1,2]),
#                                               str(blocksAlign.iat[j,3]),
#                                               blocksAlign.iat[0,5],
#                                               blocksAlign.iat[0,6],
#                                               "Tandem:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(eCoord)]) + "\n")
#
#     fout.close()
#     return None


def getSV(cwdPath, allAlignments, prefix, offset):
    """
    inverted regions are output in reverse as well
    :param cwdPath:
    :param allAlignments:
    :param prefix:
    :param offset:
    :return:
    """
    offset = -abs(offset)
    fout = open(cwdPath + prefix + "sv.txt", "w")
    allAlignments["id"] = allAlignments.group.astype(
        "str") + allAlignments.aChr + allAlignments.bChr + allAlignments.state
    allBlocks = pd.unique(allAlignments.id)

    for i in allBlocks:
        blocksAlign = allAlignments.loc[allAlignments.id == i].copy()
        if len(pd.unique(blocksAlign["aChr"])) > 1 or len(pd.unique(blocksAlign["aChr"])) > 1:
            sys.exit("More than one chromosome found for a SR")
        fout.write("\t".join(["#",
                              str(blocksAlign[["aStart","aEnd"]].min().min()),
                              str(blocksAlign[["aStart", "aEnd"]].max().max()),
                              str(blocksAlign[["bStart", "bEnd"]].min().min()),
                              str(blocksAlign[["bStart", "bEnd"]].max().max()),
                              pd.unique(blocksAlign["aChr"])[0],
                              pd.unique(blocksAlign["bChr"])[0]]) + "\n")
        ordered = 1 if "inv" not in blocksAlign.state.iloc[0] else 0
        for j in range(len(blocksAlign) - 1):
            m = blocksAlign.iat[j + 1, 0] - blocksAlign.iat[j, 1] - 1
            if ordered:
                n = blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j, 3] - 1
            else:
                n = blocksAlign.iat[j, 3] - blocksAlign.iat[j + 1, 2] - 1

            if offset <= m <= 0:  ## No significant overlap of reference genome

                if offset <= n <= 0:  ## No significant overlap in query genome
                    continue

                elif n > 0:
                    if ordered:
                        fout.write("\t".join(["INS",
                                              str(blocksAlign.iat[j, 1]),
                                              str(blocksAlign.iat[j + 1, 0]),
                                              str(blocksAlign.iat[j, 3] + 1),
                                              str(blocksAlign.iat[j + 1, 2] - 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        fout.write("\t".join(["INS",
                                              str(blocksAlign.iat[j, 1]),
                                              str(blocksAlign.iat[j + 1, 0]),
                                              str(blocksAlign.iat[j, 3] - 1),
                                              str(blocksAlign.iat[j + 1, 2] + 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

                elif n < offset:
                    if ordered:
                        j_prop = abs(n) / (blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2])
                        j1_prop = abs(n) / (blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2])
                        sCoord = round(
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        j_prop = abs(n) / (blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3])
                        j1_prop = abs(n) / (blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3])
                        sCoord = round(
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

            elif m == 1:

                if offset <= n <= 0:
                    if ordered:
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3]),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3]),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

                elif n == 1:
                    if ordered:
                        fout.write("\t".join(["SNP",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3] + 1),
                                              str(blocksAlign.iat[j + 1, 2] - 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        fout.write("\t".join(["SNP",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3] - 1),
                                              str(blocksAlign.iat[j + 1, 2] + 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                elif n > 1:
                    if ordered:
                        fout.write("\t".join(["HDR",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3] + 1),
                                              str(blocksAlign.iat[j + 1, 2] - 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        fout.write("\t".join(["HDR",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3] - 1),
                                              str(blocksAlign.iat[j + 1, 2] + 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

                elif n < offset:
                    if ordered:
                        j_prop = abs(n) / (blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2])
                        j1_prop = abs(n) / (blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2])
                        sCoord = round(
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        j_prop = abs(n) / (blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3])
                        j1_prop = abs(n) / (blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3])
                        sCoord = round(
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

            elif m > 1:
                if offset <= n <= 0:
                    if ordered:
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3]),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3]),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                elif n > 0:
                    if ordered:
                        fout.write("\t".join(["HDR",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3] + 1),
                                              str(blocksAlign.iat[j + 1, 2] - 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        fout.write("\t".join(["HDR",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              str(blocksAlign.iat[j, 3] - 1),
                                              str(blocksAlign.iat[j + 1, 2] + 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                elif n < offset:
                    if ordered:
                        j_prop = abs(n) / (blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2])
                        j1_prop = abs(n) / (blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2])
                        sCoord = round(
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        j_prop = abs(n) / (blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3])
                        j1_prop = abs(n) / (blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3])
                        sCoord = round(
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

            elif m < offset:

                j_prop = abs(m) / (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])
                j1_prop = abs(m) / (blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])

                if n >= offset:
                    if ordered:
                        sCoord = round(
                            blocksAlign.iat[j, 3] - j_prop * (blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 2] + j1_prop * (
                                    blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2])).astype(int)
                        fout.write("\t".join(["CPG",
                                              str(blocksAlign.iat[j + 1, 0]),
                                              str(blocksAlign.iat[j, 1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        sCoord = round(
                            blocksAlign.iat[j, 3] + j_prop * (blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3])).astype(
                            int)
                        eCoord = round(blocksAlign.iat[j + 1, 2] - j1_prop * (
                                    blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3])).astype(int)
                        fout.write("\t".join(["CPG",
                                              str(blocksAlign.iat[j + 1, 0]),
                                              str(blocksAlign.iat[j, 1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

                if n < offset:
                    maxOverlap = max(abs(m), abs(n))
                    if abs(m - n) < 0.1 * maxOverlap:  ## no SV if the overlap on both genomes is of similar size
                        continue

                    if abs(m) > abs(n):
                        if ordered:
                            sCoord = round(blocksAlign.iat[j, 3] - j_prop * (
                                        blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2])).astype(int)
                            eCoord = round(blocksAlign.iat[j + 1, 2] + j1_prop * (
                                        blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2])).astype(int)
                            fout.write("\t".join(["TDM",
                                                  str(blocksAlign.iat[j + 1, 0]),
                                                  str(blocksAlign.iat[j, 1]),
                                                  str(sCoord),
                                                  str(eCoord),
                                                  blocksAlign.iat[0, 5],
                                                  blocksAlign.iat[0, 6]]) + "\n")

                        else:
                            sCoord = round(blocksAlign.iat[j, 3] + j_prop * (
                                        blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3])).astype(int)
                            eCoord = round(blocksAlign.iat[j + 1, 2] - j1_prop * (
                                        blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3])).astype(int)
                            fout.write("\t".join(["TDM",
                                                  str(blocksAlign.iat[j + 1, 0]),
                                                  str(blocksAlign.iat[j, 1]),
                                                  str(sCoord),
                                                  str(eCoord),
                                                  blocksAlign.iat[0, 5],
                                                  blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        if ordered:
                            k_prop = abs(n) / (blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2])
                            k1_prop = abs(n) / (blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2])
                            sCoord = round(blocksAlign.iat[j, 1] - k_prop * (
                                        blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(int)
                            eCoord = round(blocksAlign.iat[j + 1, 0] + k1_prop * (
                                        blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                            fout.write("\t".join(["TDM",
                                                  str(sCoord),
                                                  str(eCoord),
                                                  str(blocksAlign.iat[j + 1, 2]),
                                                  str(blocksAlign.iat[j, 3]),
                                                  blocksAlign.iat[0, 5],
                                                  blocksAlign.iat[0, 6]]) + "\n")
                        else:
                            k_prop = abs(n) / (blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3])
                            k1_prop = abs(n) / (blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3])
                            sCoord = round(blocksAlign.iat[j, 1] - k_prop * (
                                        blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0])).astype(int)
                            eCoord = round(blocksAlign.iat[j + 1, 0] + k1_prop * (
                                        blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0])).astype(int)
                            fout.write("\t".join(["TDM",
                                                  str(sCoord),
                                                  str(eCoord),
                                                  str(blocksAlign.iat[j + 1, 2]),
                                                  str(blocksAlign.iat[j, 3]),
                                                  blocksAlign.iat[0, 5],
                                                  blocksAlign.iat[0, 6]]) + "\n")

    fout.close()
    return None


def runss(_id, _sspath, _delta, allAlignments):
    from subprocess import Popen, PIPE

    _block = allAlignments.loc[allAlignments.id == _id].copy()

    if 1 < len(pd.unique(_block["aChr"])) or 1 < len(pd.unique(_block["bChr"])):
            sys.exit("More than one chromosome found for a SR")
    _p = Popen([_sspath + " -HrTS " + _delta], stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    _out = _p.communicate(input=_block[["aStart", "aEnd", "bStart", "bEnd", "aChr", "bChr"]].to_string(index=False, header=False).encode())
    return "\t".join(["#",
                      str(_block[["aStart", "aEnd"]].min().min()),
                      str(_block[["aStart", "aEnd"]].max().max()),
                      str(_block[["bStart", "bEnd"]].min().min()),
                      str(_block[["bStart", "bEnd"]].max().max()),
                      pd.unique(_block["aChr"])[0],
                      pd.unique(_block["bChr"])[0]])+ "\n" + _out[0].decode("UTF-8")

def getshv(args):

    # fin = "mumSNPIn.txt" if args.align.name is None else args.align.name
    cwdpath = args.dir
    prefix = args.prefix
    nc = args.nCores
    buff = args.buff
    sspath = args.sspath
    delta = args.delta.name

    allAlignments = readSRData(cwdpath, prefix, args.all)
    allAlignments["id"] = allAlignments.group.astype("str") + allAlignments.aChr + allAlignments.bChr + allAlignments.state
    allBlocks = pd.unique(allAlignments.id)

    blocklists = [allBlocks[_i:(_i+nc)] for _i in range(0, len(allBlocks), nc)]
    # count = 0
    with open(cwdpath + prefix + "snps.txt", "w") as fout:
        for _id in blocklists:
            with Pool(processes=nc) as pool:
                out = pool.map(partial(runss, _sspath=sspath, _delta=delta, allAlignments=allAlignments), _id)
            for _snp in out:
                fout.write(_snp)

    if buff > 0:
        with open("snps.txt", "r") as fin:
            with open("snps_buff"+str(buff)+".txt", "w") as fout:
                for line in fin:
                    if line[0] == "#":
                        fout.write(line)
                    else:
                        _l = line.strip().split("\t")
                        if _l[1] != "." and _l[2] != "." and int(_l[4]) < buff:
                            continue
                        else:
                            fout.write(line)
    return None

def getNotAligned(cwdPath, prefix, ref, qry):
    refSize = {fasta.id: len(fasta.seq) for fasta in parse(ref,'fasta')}
    qrySize = {fasta.id: len(fasta.seq) for fasta in parse(qry,'fasta')}
    
    annoCoords = pd.DataFrame()
    for fileType in ["syn","inv", "TL", "invTL","dup", "invDup"]:
        try:
            fileData = pd.read_table(cwdPath+prefix+fileType+"Out.txt", header=None, dtype = object)  
        except pd.errors.ParserError as e:
            fileData = pd.read_table(cwdPath+prefix+fileType+"Out.txt", header=None, dtype = object, engine ="python")
        except pd.io.common.EmptyDataError:
            print(fileType, "Out.txt is empty. Skipping analysing it.")
            continue        
        except Exception as e:
            print("ERROR: while trying to read ", fileType, "Out.txt", e)
            continue
        coordsData = fileData.loc[fileData[0] == "#",[2,3,6,7,1,5]].copy()
        coordsData[[2,3,6,7]] = coordsData[[2,3,6,7]].astype(dtype="int64")
        coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
        annoCoords = annoCoords.append(coordsData.copy())

    try:
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header = None, dtype = object)
    except pd.errors.ParserError as e:
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header=None, dtype = object, engine ="python")
    except pd.io.common.EmptyDataError:
        print("ctxOut.txt is empty. Skipping analysing it.")
    except Exception as e:
        print("ERROR: while trying to read ", fileType, "Out.txt", e)
#    fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header = None, names = list(range(11)), dtype = object, sep ="\t")
    coordsData = fileData.loc[fileData[0] == "#"]
    coordsData = coordsData[[2,3,6,7,1,5]].copy()
    coordsData[[2,3,6,7]] = coordsData[[2,3,6,7]].astype(dtype="int64")
    coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
    
    annoCoords = annoCoords.append(coordsData.copy())
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))
  
    with open(cwdPath + prefix+"notAligned.txt","w") as fout:
        df = annoCoords[["aStart","aEnd","aChr"]].copy()
        df.sort_values(["aChr", "aStart", "aEnd"], inplace = True)
        for chrom in sorted(annoCoords.aChr.unique()):
            chromData = df.loc[df.aChr == chrom]
            maxEnd = chromData.iloc[0,1]
            if chromData.iloc[0,0] > 1:
                fout.write("\t".join(["R",str(1),
                                          str(chromData.iloc[0,0] - 1),
                                          chrom]) + "\n")
            for row in chromData.itertuples(index = False):
                if row.aStart > maxEnd+1:
                    fout.write("\t".join(["R",str(maxEnd+1),
                                          str(row.aStart - 1),
                                          chrom]) + "\n")
                if row.aEnd > maxEnd:
                    maxEnd = row.aEnd
            if maxEnd < refSize[chrom]:
                fout.write("\t".join(["R",str(maxEnd+1),
                                          str(refSize[chrom]),
                                          chrom]) + "\n")
            
        
        df = annoCoords[["bStart","bEnd","bChr"]].copy()
        df.sort_values(["bChr", "bStart", "bEnd"], inplace = True)
        for chrom in sorted(annoCoords.bChr.unique()):
            chromData = df.loc[df.bChr == chrom]
            maxEnd = chromData.iloc[0,1]
            if chromData.iloc[0,0] > 1:
                fout.write("\t".join(["Q",str(1),
                                          str(chromData.iloc[0,0] - 1),
                                          chrom]) + "\n")
            for row in chromData.itertuples(index = False):
                if row.bStart > maxEnd+1:
                    fout.write("\t".join(["Q",str(maxEnd+1),
                                          str(row.bStart - 1),
                                          chrom]) + "\n")
                if row.bEnd > maxEnd:
                    maxEnd = row.bEnd
            if maxEnd < qrySize[chrom]:
                fout.write("\t".join(["Q",str(maxEnd+1),
                                          str(qrySize[chrom]),
                                          chrom]) + "\n")
                
#    fout.close()
    return None

##################################################################
### Generate formatted output
##################################################################

def getsrtable(cwdpath):
    import re

    files = ["synOut.txt", "invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt", "ctxOut.txt"]
    entries = defaultdict()

    re_dup = re.compile("dup", re.I)
    align_num = 1
    sr_num = 0
    for file in files:
        sr_type = file.split("O")[0]
        if sr_type == "syn":
            sr_type = "SYN"
        elif sr_type == "inv":
            sr_type = "INV"
        elif sr_type == "TL":
            sr_type = "TRANS"
        elif sr_type == "invTL":
            sr_type = "INVTR"
        elif sr_type == "dup":
            sr_type = "DUP"
        elif sr_type == "invDup":
            sr_type = "INVDP"

        with open(cwdpath + file, "r") as fin:
            for line in fin:
                line = line.strip().split("\t")
                if line[0] == "#":
                    sr_num += 1
                    if file == 'ctxOut.txt':
                        if line[8] == "translocation":
                            sr_type = "TRANS"
                        elif line[8] == "invTranslocation":
                            sr_type = "INVTR"
                        elif line[8] == "duplication":
                            sr_type = "DUP"
                        elif line[8] == "invDuplication":
                            sr_type = "INVDP"
                    entries[sr_type + str(sr_num)] = {
                        'achr': line[1],
                        'astart': line[2],
                        'aend': line[3],
                        'bchr': line[5],
                        'bstart': line[6],
                        'bend': line[7],
                        'vartype': sr_type,
                        'parent': "-",
                        'dupclass': "-",
                        'aseq': "-",
                        "bseq": "-"
                    }
                    if re_dup.search(file):
                        entries[sr_type + str(sr_num)]["dupclass"] = "copygain" if line[9] == "B" else "copyloss"
                    if file == "ctxOut.txt":
                        entries[sr_type + str(sr_num)]["vartype"] = sr_type
                        if re_dup.search(line[8]):
                            entries[sr_type + str(sr_num)]["dupclass"] = "copygain" if line[9] == "B" else "copyloss"
                else:
                    entries[sr_type + "AL" + str(align_num)] = {
                        'achr': entries[sr_type + str(sr_num)]['achr'],
                        'astart': line[0],
                        'aend': line[1],
                        'bchr': entries[sr_type + str(sr_num)]['bchr'],
                        'bstart': line[2],
                        'bend': line[3],
                        'vartype': sr_type + "AL",
                        'parent': sr_type + str(sr_num),
                        'dupclass': "-",
                        'aseq': "-",
                        "bseq": "-"
                    }
                    align_num += 1

    anno = pd.DataFrame.from_dict(entries, orient="index")
    anno.loc[:, ['astart', 'aend', 'bstart', 'bend']] = anno.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
    anno['id'] = anno.index.values
    anno = anno.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
    anno.sort_values(['achr', 'astart', 'aend'], inplace=True)
    return anno


def extractseq(_gen, _pos):
    chrs = defaultdict(dict)
    for fasta in parse(_gen, 'fasta'):
        if fasta.id in _pos.keys():
            chrs[fasta.id] = {_i:fasta.seq[_i-1] for _i in _pos[fasta.id]}
    return chrs

def getTSV(cwdpath, ref):
    """
    :param cwdpath: Path containing all input files
    :return: A TSV file containing genomic annotation for the entire genome
    """

    import pandas as pd
    import sys
    from collections import defaultdict
    logger = logging.getLogger("getTSV")
    anno = getsrtable(cwdpath)

    _svdata = pd.read_table(cwdpath + "sv.txt", header=None)
    _svdata.columns = ["vartype", "astart", 'aend', 'bstart', 'bend', 'achr', 'bchr']

    entries = defaultdict()
    count = 1
    for row in _svdata.itertuples(index=False):
        if row.vartype == "#":
            _parent = anno.loc[(anno.achr == row.achr) &
                               (anno.astart == row.astart) &
                               (anno.aend == row.aend) &
                               (anno.bchr == row.bchr) &
                               (anno.bstart == row.bstart) &
                               (anno.bend == row.bend) &
                               (anno.parent == "-"), "id"]
            if len(_parent) != 1:
                logger.error("Error in finding parent for SV")
                logger.error(row.to_string() + "\t" + _parent.to_string())
                sys.exit()
            else:
                _parent = _parent.to_string(index=False, header=False)
            continue
        entries[row.vartype + str(count)] = {
            'achr': row.achr,
            'astart': row.astart,
            'aend': row.aend,
            'bchr': row.bchr,
            'bstart': row.bstart,
            'bend': row.bend,
            'vartype': row.vartype,
            'parent': _parent,
            'id': row.vartype + str(count),
            'dupclass': "-",
            'aseq': "-",
            "bseq": "-"
        }
        count += 1

    sv = pd.DataFrame.from_dict(entries, orient="index")
    sv.loc[:, ['astart', 'aend', 'bstart', 'bend']] = sv.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
    sv = sv.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
    sv.sort_values(['achr', 'astart', 'aend'], inplace=True)

    # out = pd.concat([anno, sv])
    # out.sort_values(['achr', 'astart', 'aend'], inplace=True)

    notal = pd.read_table(cwdpath + "notAligned.txt", header=None)
    notal.columns = ["gen", "start", "end", "chr"]
    notal[["start", "end"]] = notal[["start", "end"]].astype("int")
    entries = defaultdict()
    _c = 0
    for row in notal.itertuples(index=False):
        _c += 1
        if row.gen == "R":
            entries["NOTAL" + str(_c)] = {
                'achr': row.chr,
                'astart': row.start,
                'aend': row.end,
                'bchr': "-",
                'bstart': "-",
                'bend': "-",
                'vartype': "NOTAL",
                'parent': "-",
                'id': "NOTAL" + str(_c),
                'dupclass': "-",
                'aseq': "-",
                "bseq": "-"
            }
        elif row.gen == "Q":
            entries["NOTAL" + str(_c)] = {
                'achr': "-",
                'astart': "-",
                'aend': "-",
                'bchr': row.chr,
                'bstart': row.start,
                'bend': row.end,
                'vartype': "NOTAL",
                'parent': "-",
                'id': "NOTAL" + str(_c),
                'dupclass': "-",
                'aseq': "-",
                'bseq': "-"
            }
    notal = pd.DataFrame.from_dict(entries, orient="index")
    # notal.loc[:, ['astart', 'aend', 'bstart', 'bend']] = notal.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
    notal = notal.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
    notal.sort_values(['achr', 'astart', 'aend'], inplace=True)
    notal['selected'] = -1

    def p_indel():
        vtype = "INS" if indel == 1 else "DEL"
        entries[vtype + str(count)] = {
            'achr': _ac,
            'astart': _as,
            'aend': _ae,
            'bchr': _bc,
            'bstart': _bs,
            'bend': _be,
            'vartype': vtype,
            'parent': _p,
            'aseq': "-" if vtype == "INS" else _seq,
            'bseq': "-" if vtype == "DEL" else _seq,
            'id': vtype + str(count),
            'dupclass': "-"
        }

    entries = defaultdict()
    with open("snps.txt", "r") as fin:
        indel = 0                           # marker for indel status. 0 = no_indel, 1 = insertion, -1 = deletion
        _as = -1                            # astart
        _ae = -1
        _bs = -1
        _be = -1
        _ac = -1
        _bc = -1
        _p = -1
        _seq = ""

        for line in fin:
            line = line.strip().split("\t")
            try:
                if line[0] == "#" and len(line) == 7:
                    if indel != 0:
                        p_indel()
                        indel = 0
                        _seq = ""
                    _parent = anno.loc[(anno.achr == line[5]) &
                                       (anno.astart == int(line[1])) &
                                       (anno.aend == int(line[2])) &
                                       (anno.bchr == line[6]) &
                                       (anno.bstart == int(line[3])) &
                                       (anno.bend == int(line[4])) &
                                       (anno.parent == "-"), "id"]
                    if len(_parent) != 1:
                        logger.error("Error in finding parent for SNP")
                        logger.error("\t".join(line) + "\t" + _parent.to_string())
                        sys.exit()
                    else:
                        _parent = _parent.to_string(index=False, header=False)
                    continue
                elif line[1] != "." and line[2] != "." and len(line) == 12:
                    if indel != 0:
                        p_indel()
                        indel = 0
                        _seq = ""
                    count += 1
                    entries["SNP" + str(count)] = {
                        'achr': line[10],
                        'astart': int(line[0]),
                        'aend': int(line[0]),
                        'bchr': line[11],
                        'bstart': int(line[3]),
                        'bend': int(line[3]),
                        'vartype': "SNP",
                        'parent': _parent,
                        'aseq': line[1],
                        'bseq': line[2],
                        'id': "SNP" + str(count),
                        'dupclass': "-"
                    }
                elif indel == 0:
                    count += 1
                    _as = int(line[0])
                    _ae = int(line[0])
                    _bs = int(line[3])
                    _be = int(line[3])
                    _ac = line[10]
                    _bc = line[11]
                    _p = _parent
                    indel = 1 if line[1] == "." else -1
                    _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                elif indel == 1:
                    if int(line[0]) != _as or line[1] != "." or line[10] != _ac:
                        p_indel()
                        _seq = ""
                        count += 1
                        _as = int(line[0])
                        _ae = int(line[0])
                        _bs = int(line[3])
                        _be = int(line[3])
                        _ac = line[10]
                        _bc = line[11]
                        _p = _parent
                        indel = 1 if line[1] == "." else -1
                        _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                    else:
                        _be = int(line[3])
                        _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                elif indel == -1:
                    if int(line[3]) != _bs or line[2] != "." or line[11] != _bc:
                        p_indel()
                        _seq = ""
                        count += 1
                        _as = int(line[0])
                        _ae = int(line[0])
                        _bs = int(line[3])
                        _be = int(line[3])
                        _ac = line[10]
                        _bc = line[11]
                        _p = _parent
                        indel = 1 if line[1] == "." else -1
                        _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                    else:
                        _ae = int(line[0])
                        _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
            except IndexError as _e:
                logger.error(_e)
                logger.error("\t".join(line))
                sys.exit()

    snpdata = pd.DataFrame.from_dict(entries, orient="index")
    snpdata.loc[:, ['astart', 'aend', 'bstart', 'bend']] = snpdata.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
    snpdata['id'] = snpdata.index.values
    snpdata = snpdata.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
    snpdata.sort_values(['achr', 'astart', 'aend'], inplace=True)

    positions = defaultdict()
    for _chr in snpdata.achr.unique():
        positions[_chr] = snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "INS"), "astart"].tolist() + (snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "DEL"), "astart"] - 1).tolist()
        # positions[chrom] = (snpdata.loc[(snpdata.achr == chrom) & (snpdata.vartype == "DEL"), "astart"] - 1).tolist()

    seq = extractseq(ref, positions)

    for _chr in snpdata.achr.unique():
        ## Fix coordinates and sequence for insertions
        _indices = snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "INS")].index.values
        _seq = pd.Series([seq[_chr][_i] for _i in snpdata.loc[_indices, "astart"]], index = _indices)
        _dir = _indices[~snpdata.loc[_indices, "parent"].str.contains("INV")]
        _inv = _indices[snpdata.loc[_indices, "parent"].str.contains("INV")]

        # Add previous characted
        snpdata.loc[_indices, "aseq"] = _seq
        snpdata.loc[_dir, "bseq"] = _seq.loc[_dir] + snpdata.loc[_dir, "bseq"]
        snpdata.loc[_inv, "bseq"] = _seq.loc[_inv] + snpdata.loc[_inv, "bseq"].str[::-1]

        # Change coordinates
        snpdata.loc[_dir, "bstart"] = snpdata.loc[_dir, "bstart"] - 1
        snpdata.loc[_inv, "bend"] = snpdata.loc[_inv, "bstart"] + snpdata.loc[_inv, "bend"]
        snpdata.loc[_inv, "bstart"] = snpdata.loc[_inv, "bend"] - snpdata.loc[_inv, "bstart"] + 1
        snpdata.loc[_inv, "bend"] = snpdata.loc[_inv, "bend"] - snpdata.loc[_inv, "bstart"] + 1

        ## Fix coordinates and sequence for deletions
        _indices = snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "DEL")].index.values
        _seq = pd.Series([seq[_chr][_i-1] for _i in snpdata.loc[_indices, "astart"]], index=_indices)
        _dir = _indices[~snpdata.loc[_indices, "parent"].str.contains("INV")]
        _inv = _indices[snpdata.loc[_indices, "parent"].str.contains("INV")]

        # Add previous characted
        snpdata.loc[_indices, "aseq"] = _seq + snpdata.loc[_indices, "aseq"]
        snpdata.loc[_indices, "bseq"] = _seq

        # Change coordinates
        snpdata.loc[_indices, "astart"] = snpdata.loc[_indices, "astart"] - 1
        snpdata.loc[_inv, "bstart"] = snpdata.loc[_inv, "bstart"] + 1
        snpdata.loc[_inv, "bend"] = snpdata.loc[_inv, "bend"] + 1

    events = anno.loc[anno.parent == "-"]
    with open("syri.out", "w") as fout:
        notA = notal.loc[notal.achr != "-"].copy()
        notA.loc[:, ["astart", "aend"]] = notA.loc[:, ["astart", "aend"]].astype("int")
        notB = notal.loc[notal.bchr != "-"].copy()
        notB.loc[:, ["bstart", "bend"]] = notB.loc[:, ["bstart", "bend"]].astype("int")
        row_old = -1
        for row in events.itertuples(index=False):
            if row_old != -1 and row_old.achr != row.achr:
                _notA = notA.loc[(notA.achr == row_old.achr) & (notA.aend == row_old.aend+1) & (notA.selected != 1), notA.columns != 'selected']
                notA.loc[(notA.achr == row_old.achr) & (notA.aend == row_old.aend + 1), 'selected'] = 1
                if len(_notA) == 0:
                    pass
                elif len(_notA) == 1:
                    fout.write("\t".join(list(map(str, _notA.iloc[0]))) + "\n")
                else:
                    logger.error("too many notA regions")
                    sys.exit()

            _notA = notA.loc[(notA.achr == row.achr) & (notA.aend == row.astart-1) & (notA.selected != 1), notA.columns != 'selected']
            notA.loc[(notA.achr == row.achr) & (notA.aend == row.aend + 1), 'selected'] = 1
            # _notB = notB.loc[(notB.bchr == row.bchr) & (notB.bend == row.bstart-1)]
            if len(_notA) == 0:
                pass
            elif len(_notA) == 1:
                fout.write("\t".join(list(map(str, _notA.iloc[0]))) + "\n")
            else:
                logger.error("too many notA regions")
                sys.exit()

            # if len(_notB) == 0:
            #     pass
            # elif len(_notB) == 1:
            #     fout.write("\t".join(list(map(str, _notB.iloc[0]))) + "\n")
            # else:
            #     sys.exit("too many notB regions")
            fout.write("\t".join(list(map(str, row))) + "\n")
            row_old = row

            outdata = pd.concat([anno.loc[anno.parent == row.id], sv.loc[sv.parent == row.id], snpdata.loc[(snpdata.parent == row.id)]])
            outdata.sort_values(["astart", "aend"], inplace=True)
            fout.write(outdata.to_csv(sep="\t", index=False, header=False))
        fout.write(notB.loc[:, notB.columns != 'selected'].to_csv(sep="\t", index=False, header=False))
    return 0


def getVCF(finname, foutname):
    """
    does not output notal in qry genome
    :param finname:
    foutname:
    :return:
    """
    data = pd.read_table(finname, header=None)
    data.columns = ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']
    data = data.loc[data['achr'] != "-"].copy()
    data.sort_values(['achr', 'astart', 'aend'], inplace=True)
    data.loc[:, ['astart', 'aend', 'bstart', 'bend']] = data.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype(str)

    with open(foutname, 'w') as fout:
        fout.write('##fileformat=VCFv4.3\n')
        fout.write('##fileDate=' + str(date.today()).replace('-', '') + '\n')
        fout.write('##source=syri\n')
        fout.write('##ALT=<ID=SYN,Description="Syntenic region">' + '\n')
        fout.write('##ALT=<ID=INV,Description="Inversion">' + '\n')
        fout.write('##ALT=<ID=TRANS,Description="Translocation">' + '\n')
        fout.write('##ALT=<ID=INVTR,Description="Inverted Translocation">' + '\n')
        fout.write('##ALT=<ID=DUP,Description="Duplication">' + '\n')
        fout.write('##ALT=<ID=INVDP,Description="Inverted Duplication">' + '\n')
        fout.write('##ALT=<ID=SYNAL,Description="Syntenic alignment">' + '\n')
        fout.write('##ALT=<ID=INVAL,Description="Inversion alignment">' + '\n')
        fout.write('##ALT=<ID=TRANSAL,Description="Translocation alignment">' + '\n')
        fout.write('##ALT=<ID=INVTRAL,Description="Inverted Translocation alignment">' + '\n')
        fout.write('##ALT=<ID=DUPAL,Description="Duplication alignment">' + '\n')
        fout.write('##ALT=<ID=INVDPAL,Description="Inverted Duplication alignment">' + '\n')
        fout.write('##ALT=<ID=HDR,Description="Highly diverged regions">' + '\n')
        fout.write('##ALT=<ID=INS,Description="Insertion in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=DEL,Description="Deletion in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=CPG,Description="Copy gain in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=CPL,Description="Copy loss in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=SNP,Description="Single nucleotide polymorphism">' + '\n')
        fout.write('##ALT=<ID=TDM,Description="Tandem repeat">' + '\n')
        fout.write('##ALT=<ID=NOTAL,Description="Not Aligned region">' + '\n')
        fout.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position on reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=ChrB,Number=1,Type=String,Description="Chromoosme ID on the non-reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=StartB,Number=1,Type=Integer,Description="Start position on non-reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=EndB,Number=1,Type=Integer,Description="End position on non-reference genome">' + '\n')
        fout.write('##INFO=<ID=Parent,Number=1,Type=String,Description="ID of the parent SR">' + '\n')
        fout.write(
            '##INFO=<ID=VarType,Number=1,Type=String,Description="Start position on non-reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=DupType,Number=1,Type=String,Description="Copy gain or loss in the non-reference genome">' + '\n')
        # fout.write('##INFO=<ID=NotAlGen,Number=1,Type=String,Description="Genome containing the not aligned region">' + '\n')

        fout.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']) + '\n')
        for line in data.itertuples(index=False):
            pos = [line[0], line[1], ".", 'N', '<' + line[10] + '>', '.', 'PASS']

            if line[10] in ["SYN", "INV", "TRANS", "INVTR", "DUP", "INVDUP"]:
                _info = ';'.join(['END='+line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent=.', 'VarType='+'SR', 'DupType='+line[11]])
                pos.append(_info)
                fout.write('\t'.join(pos) + '\n')

            elif line[10] == "NOTAL":
                if line[0] != "-":
                    _info = ';'.join(
                        ['END=' + line[2], 'ChrB=.', 'StartB=.', 'EndB=.', 'Parent=.', 'VarType=.', 'DupType=.'])
                    pos.append(_info)
                    fout.write('\t'.join(pos) + '\n')

            elif line[10] in ["SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDUPAL"]:
                _info = ';'.join(
                    ['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=.', 'DupType=.'])
                pos.append(_info)
                fout.write('\t'.join(pos) + '\n')

            elif line[10] in ['CPG', 'CPL', 'TDM', 'HDR']:
                _info = ";".join(['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=ShV', 'DupType=.'])
                pos.append(_info)
                fout.write('\t'.join(pos) + '\n')

            elif line[10] in ['SNP', 'DEL', 'INS']:
                if (line[3] == '-') != (line[4] == '-'):
                    sys.exit("Inconsistency in annotation type. Either need seq for both or for none.")
                elif line[3] == '-' and line[4] == '-':
                    _info = ";".join(['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=ShV', 'DupType=.'])
                    pos.append(_info)
                    fout.write('\t'.join(pos) + '\n')
                elif line[3] != '-' and line[4] != '-':
                    pos = [line[0], line[1], ".", line[3], line[4], '.', 'PASS']
                    _info = ";".join(['END=' + line[2], 'ChrB=' + line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent=' + line[9], 'VarType=ShV', 'DupType=.'])
                    pos.append(_info)
                    fout.write('\t'.join(pos) + '\n')
    return 0



##################################################################
### Multi SV functions
##################################################################
        
def getBlocksData(filePaths, fileTypes):
    genomeData = pd.DataFrame()
    synData = pd.read_table(filePaths[0],header = None)
    synData = synData.loc[synData[0] != "#"]
    synData = synData[[0,1,2,3]]
    synData.columns = ["aStart", "aEnd","bStart","bEnd"]
    synData = synData.astype("int64")
    synData["state"] = fileTypes.pop(0)
    genomeData = genomeData.append(synData)
    
    for i in filePaths[1:]:
        fileData = pd.read_table(i, header = None)
        fileData = fileData.loc[fileData[0] == "#"]
        fileData = fileData[[1,2,4,5]]
        fileData.columns = ["aStart", "aEnd","bStart","bEnd"]
        fileData = fileData.astype("int64")
        fileData["state"] = fileTypes.pop(0)
        genomeData = genomeData.append(fileData)
    genomeData.sort_values(by = ["aStart", "aEnd","bStart","bEnd"], inplace = True)
    genomeData.index = range(len(genomeData))
    return(genomeData)
    

def getConservedRegions(dataList, isSyn = True):
    if not isinstance(isSyn, bool):
        raise TypeError("isSyn must be a bool")
    
    bGenomes = np.unique(dataList.bGenome)
    genomeCount = len(bGenomes)
    if isSyn:
        allCoords = dataList[["start","end","bGenome"]]
    else:
        allCoords = pd.DataFrame()
        for i in bGenomes:
            gData = dataList.loc[dataList.bGenome == i]
            start = list(gData.start)
            end = list(gData.end)
            s = start[0]
            e = end[0]
            endStack = [float('inf'),e]
            region = []
            for j in range(1,len(start)):
                s1 = start[j]
                e1 = end[j]
                if s1 < e:
                    region.append([s,s1])
                    s = s1
                    endStack.append(e1)
                    e = min(endStack)
                elif e <= s1:
                    region.append([s,e])
                    while True:
                        s = e
                        endStack.remove(e)
                        e = min(endStack)
                        if e <= s1:
                            region.append([s,e])
                        else:
                            break
                    if len(endStack) > 1:
                        region.append([s,s1])
                        s = s1
                        endStack.append(e1)
                        e = min(endStack)
                    else:
                        s = s1
                        endStack.append(e1)
                        e = min(endStack)
            while len(endStack) > 1:
                region.append([s,e])
                s = e
                endStack.remove(e)
                e = min(endStack)
            region = [a for a in region if a[0] < a[1]]
            region = pd.DataFrame(region)
            region.columns = ["start","end"]
            region["bGenome"] = i
            allCoords = allCoords.append(region)
        allCoords.sort_values(["start","end"],inplace = True)
        
    terminalData = pd.DataFrame(data = np.zeros([2,genomeCount], dtype = "int"), index = ["start","end"], columns = bGenomes)
    
    inGenomeCount = 0 
    regions = []
    count = 0
    for row in allCoords.itertuples(index=False):
        count+=1
        if row.end <= terminalData[row.bGenome].end:
            print("values must be sorted. Invalid Entry: ",row, terminalData[row.bGenome])
            sys.exit()
        if row.start <= terminalData[row.bGenome].end:
            terminalData[row.bGenome].start = terminalData[row.bGenome].end + 1
            terminalData[row.bGenome].end = row.end
        else:
            terminalData[row.bGenome].start = row.start
            terminalData[row.bGenome].end = row.end
        if max(terminalData.loc["start"]) < min(terminalData.loc["end"]):
            regions.append((max(terminalData.loc["start"]),min(terminalData.loc["end"])))
    regions = pd.DataFrame(regions, columns = ["start","end"])
    return regions


def getDataList(dataTables, genomeID, identity):
    start = []
    end = []
    bGenome = []
    bStart = []
    bEnd = []
    state = []
    
    if len(dataTables) != len(genomeID):
        print("need 1 identifier for each table")
        sys.exit()
    else:
        for i in range(len(genomeID)):
            if identity[i] == "a":
                start.extend(dataTables[i].aStart.tolist())
                end.extend(dataTables[i].aEnd.tolist())
                bStart.extend(dataTables[i].bStart.tolist())
                bEnd.extend(dataTables[i].bEnd.tolist())
            elif identity[i] == "b":
                start.extend(dataTables[i].bStart.tolist())
                end.extend(dataTables[i].bEnd.tolist())
                bStart.extend(dataTables[i].aStart.tolist())
                bEnd.extend(dataTables[i].aEnd.tolist())
            state.extend(dataTables[i].state.tolist())
            bGenome.extend([genomeID[i]]*len(dataTables[i]))
        outData = pd.DataFrame({"start":start,"end":end,"bStart":bStart,"bEnd":bEnd,"state":state,"bGenome":bGenome})
        outData.sort_values(["start","end","bStart","bEnd"], inplace = True)
        outData = outData[["start","end","bStart","bEnd","state","bGenome"]]
        outData.index = range(len(outData))
        return(outData)
        

def getCLQ(adjM,partSize):
    """
    Mirghorbani, M., & Krokhmal, P. (2013). On finding k-cliques in k-partite graphs. Optim Lett, 7, 1155–1165. https://doi.org/10.1007/s11590-012-0536-y
    """
    class startCLQ:
        def __init__(self, adjM, partSize):
            self.clq = []
            self.t = 0
            self.sub = []
            self.BsOut = list(range(len(partSize)))
            self.Bs = []
            self.Z = np.ones([len(partSize),sum(partSize)], dtype = "bool")
            self.Z0 = np.ones([len(partSize), sum(partSize)])
            self.adjM = adjM
            self.partSize = partSize
            self.partIndex = [sum(partSize[:i]) for i in range(len(partSize))]
            self.clqSize = len(partSize)
            self.S = []
            self.Q = []
            self.bitCLQ(self.t)
            
        def getBits(self, t, b):
            return self.Z[t, self.partIndex[b]:self.partIndex[b]+self.partSize[b]]

        def getBt(self, partSize, t):
            nodeCount = [sum(self.getBits(t,i)) for i in self.BsOut]
            return self.BsOut[nodeCount.index(min(nodeCount))]
        
        def bitCLQ(self, t):
            bt = self.getBt(self.partSize, t)
            sigBits = np.where(np.array(self.getBits(t,bt)) == True)[0]
            sigBitsLen = len(sigBits)
            count = 0
            for i in sigBits:
                count+=1
                if t == 0:
                    print(count, sigBitsLen, datetime.datetime.now())
                nt = self.partIndex[bt]+i
                self.Z[t,nt] = 0
                self.S.append(nt)
                if len(self.S) == self.clqSize:
                    self.Q.append(self.S.copy())
                    self.S.remove(nt)
                else:
                    self.Z[t+1] = self.Z[t] & self.adjM[nt]
                    self.Bs.append(bt)
                    self.BsOut.remove(bt)
                    P = sum([1 for i in self.BsOut if sum(self.getBits(t,i)) > 0])
                    if len(self.S) + P == self.clqSize:
                        self.bitCLQ(t+1)
                        self.S.remove(nt)
                        self.Bs.remove(bt)
                        self.BsOut.append(bt)
                    else:
                        self.S.remove(nt)
                        self.Bs.remove(bt)
                        self.BsOut.append(bt)
        
    def filterCLQ(clqList, partIndex):
        clqList = [sorted(i) for i in clqList]
        clqList = [[i[j] - partIndex[j] for j in range(len(i))] for i in clqList]
        return(clqList)
    CLQData = startCLQ(adjM, partSize)
    return(filterCLQ(CLQData.Q, CLQData.partIndex))
    

##%%
#def plotBlocks(blocksData):
#    blocksData = [orderedBlocks.iloc[[250,251,255]], orderedBlocks.iloc[[1370, 1371]]]
#    
#    blocksDataOri = [i.copy() for i in blocksData]
#
#    
#    blocksCoords = {}
#    for i in range(len(blocksData)):        
#        bData = blocksData[i].copy()
#        aMin = min(bData[['aStart','aEnd']].min())
#        aMax = max(bData[['aStart','aEnd']].max())
#        bMin = min(bData[['bStart','bEnd']].min())
#        bMax = max(bData[['bStart','bEnd']].max())
#        blocksCoords[i] = [aMin,aMax, bMin, bMax]        
#    
#    keyOrder = sorted(blocksCoords.keys(), key = lambda x : blocksCoords[x][0])
#    
#    gapsRemoved = []
##    startPosition = blocksCoords[keyOrder[0]][0]
#    maxEnd = blocksCoords[keyOrder[0]][1]
#    for i in range(1,len(keyOrder)):
#        if blocksCoords[keyOrder[i]][0] > maxEnd:
#            gapsRemoved.append(blocksCoords[keyOrder[i]][0] - blocksCoords[keyOrder[i-1]][1])
#        else:
#            gapsRemoved.append(0)
#            
#    for i in range(1, len(keyOrder)):
#        leftShift = sum(gapsRemoved[:i])
#        rightShift = max(0,log(leftShift,1.015))
#        blocksData[keyOrder[i]]['aStart'] = blocksData[keyOrder[i]]['aStart'] -  leftShift + rightShift
#        blocksData[keyOrder[i]]['aEnd'] = blocksData[keyOrder[i]]['aEnd'] -  leftShift + rightShift
#    
#    dataLimits = [[min(i[['aStart','aEnd']].min()), max(i[['aStart','aEnd']].max()),min(i[['bStart','bEnd']].min()), max(i[['bStart','bEnd']].max())] for i in blocksData]
#    aMin = min(unlist([x[0:2] for x in dataLimits]))
#    aMax = max(unlist([x[0:2] for x in dataLimits]))
#    bMin = min(unlist([x[2:4] for x in dataLimits]))
#    bMax = max(unlist([x[2:4] for x in dataLimits]))
#
#    
##    
##    for i in range(1,len(blocksCoords)):
##        if blocksCoords[i][0] > blocksCoords[i-1][1] 
##    
#    colList = ['aStart', 'aEnd', 'bStart', 'bEnd']    
#    for bData in blocksData:
#        for i in colList[:2]:
#            bData[i] = (bData[i] - aMin)/(aMax - aMin)
#            
#        for i in colList[2:]:
#            bData[i] = (bData[i] - bMin)/(bMax - bMin)
#    
#    colors = getColors(plt.cm.Dark2,len(blocksData))
#    
#    bLen = len(blocksData)
#    
#    for count in range(bLen):
#        bData = blocksData[count]
#        for i in range(bData.shape[0]):
#            row = bData.iloc[i]
#            plt.plot([row[0], row[1]], [0 - (0.1*count),0 - (0.1*count)], linewidth = 5, color = colors[count], path_effects = [pe.Stroke(linewidth = 10, foreground="k"),pe.Normal()])
#            plt.plot([row[2], row[3]], [1 + 0.1*(bLen-1-count),1 + 0.1*(bLen-1-count)], linewidth = 5, color = colors[count],path_effects = [pe.Stroke(linewidth = 10, foreground="k"),pe.Normal()])
#    plt.show()
#    plt.gca().add_patch(patches.Rectangle((10,10), 30,5))
#    plt.show()
#    plt.plot()
