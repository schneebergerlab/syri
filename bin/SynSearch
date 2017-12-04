#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 16:10:48 2017

@author: goel
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:05:51 2017

@author: goel
"""
import sys
import pandas as pd
import numpy as np
from igraph import *
from syri.methods.myUsefulFunctions import *
from syri.methods.synsearchFunctions import *
from collections import Counter
from matplotlib import pyplot as plt
from multiprocessing import Pool
from functools import partial
from time import time
import  multiprocessing
import os



def SynSearch(chromo, threshold, coords, cwdPath):
    coordsData = coords[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == 1)]
    print(coordsData.shape)
    nrow = coordsData.shape[0]
    
#    print("Accessing possible syn blocks")
    
    df = coordsData.apply( lambda X : coordsData.apply(TS, axis = 1, args = (X, threshold)), axis = 1)
    nrow = df.shape[0]
    synTree = {}
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
    synData = coordsData.iloc[synPath]
    
    
    
    
    
    ######################################################################
#    print("Finding Inversions")
##########################################################################
    
    invertedCoordsOri, profitable, bestInvPath, invData, synInInv = getInversions(coords,chromo, threshold, synData, synPath)
    

    ##########################################################
    #### Identify Translocation and duplications
    ##########################################################
    print("Finding translocation and duplication")
    
    chromBlocks = coords[(coords.aChr == chromo) & (coords.bChr == chromo)]
    
    inPlaceIndices = sorted(list(synData.index.values) + list(invData.index.values))
    
    inPlaceBlocks = chromBlocks[chromBlocks.index.isin(sorted(list(synData.index.values)))].copy()
    for i in bestInvPath:
        invPos = profitable[i-1].invPos
        invBlockData = invertedCoordsOri.iloc[invPos]
        invCoord = [invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]]
        invCoord.append(invCoord[1] - invCoord[0])
        invCoord.append(invCoord[3] - invCoord[2])
        invCoord.append(sum((invBlockData.aLen+invBlockData.bLen)*invBlockData.iden)/(invCoord[-2] + invCoord[-1]))
        invCoord.extend([1,-1,chromo,chromo])
        for j in range(profitable[i-1].neighbours[0]+1,profitable[i-1].neighbours[1]):
            inPlaceBlocks = inPlaceBlocks[inPlaceBlocks.index != synData.iloc[j].name]
            inPlaceIndices.remove(synData.iloc[j].name)
        inPlaceBlocks = inPlaceBlocks.append(pd.Series(invCoord, index = inPlaceBlocks.columns), ignore_index = True)
        
    inPlaceBlocks.sort_values(["aChr","aStart","aEnd","bChr","bStart","bEnd"], inplace = True)
    inPlaceBlocks.index = range(inPlaceBlocks.shape[0])
        

    
#    inPlaceBlocks = chromBlocks[chromBlocks.index.isin(sorted(list(synData.index.values) + list(invData.index.values)))]
    outPlaceBlocks = chromBlocks[~chromBlocks.index.isin(inPlaceIndices)]
    
    ## Should not filter redundant alignments as they "can" be part of bigger translocations
    ## filtering them may lead to removal of those translocations
    #redundant = getRedundantIndex(inPlaceBlocks, outPlaceBlocks, threshold)
    
    outPlaceBlocksFiltered = outPlaceBlocks.copy() #drop(outPlaceBlocks.index[redundant])
    orderedBlocks = outPlaceBlocksFiltered[outPlaceBlocksFiltered.bDir == 1]
    invertedBlocks = outPlaceBlocksFiltered[outPlaceBlocksFiltered.bDir == -1]
    
    ## Create connectivity tree for directed blocks
    transBlocksNeighbours = getTransSynOrientation(inPlaceBlocks, orderedBlocks, threshold)
    outOrderedBlocks = makeBlocksTree(inPlaceBlocks, orderedBlocks, threshold, transBlocksNeighbours)

    ## Create connectivity tree for inverted blocks
    
    invertedCoords = invertedBlocks.copy()
    invertedCoords.bStart = invertedCoords.bStart + invertedCoords.bEnd
    invertedCoords.bEnd = invertedCoords.bStart - invertedCoords.bEnd
    invertedCoords.bStart = invertedCoords.bStart - invertedCoords.bEnd
    invTransBlocksNeighbours = getTransSynOrientation(inPlaceBlocks, invertedBlocks, threshold)
    
    invertedCoords = invertedBlocks.copy()
#    minCoords = np.min(np.min(invertedCoords[["bStart","bEnd"]]))
    maxCoords = np.max(np.max(invertedCoords[["bStart","bEnd"]]))
    invertedCoords.bStart = maxCoords + 1 - invertedCoords.bStart 
    invertedCoords.bEnd = maxCoords + 1 - invertedCoords.bEnd
    outInvertedBlocks = makeBlocksTree(inPlaceBlocks, invertedCoords, threshold, invTransBlocksNeighbours)
    
    ## find all translocations which don't have large gaps between its alignments
    ## and are not overlappign with the syntenic blocks
    orderedBlocksTree = outOrderedBlocks 
    orderedBlocksList = makeBlocksList(orderedBlocksTree, orderedBlocks)
    outOG = getConnectivityGraph(orderedBlocksList)
    
    shortestOutOG = []
    for i in range(len(orderedBlocksList)):
        eNode = [i]
        eNode.extend(list(np.where(orderedBlocksTree.iloc[i] == True)[0]))
        shortestOutOG.append(getAllLongestPaths(outOG,i,eNode,"weight"))
    
    
    
    
#    ## find all shortest paths
#    shortestOutOG = []
#    if len(outOG.es) > 0:
#        for i in range(len(orderedBlocksList)):
#            eNode = [i]
#            eNode.extend(list(np.where(orderedBlocksTree.iloc[i] == True)[0]))
#            shortestOutOG.append(getAllLongestPaths(outOG,i,eNode,"weight"))
#            shortestOutOG = shortestOutOG + [outOG.get_all_shortest_paths(i,[i].extend(list(np.where(orderedBlocksTree.iloc[i] == True)[0])), "weight")]
        
    transScores = getTranslocationScore(shortestOutOG, orderedBlocks)
    transBlocks = getTransBlocks(transScores, shortestOutOG, orderedBlocks, inPlaceBlocks, threshold)
    
    
    
    invertedBlocksTree = outInvertedBlocks
    invertedBlocksList = makeBlocksList(invertedBlocksTree, invertedBlocks)
    outIG = getConnectivityGraph(invertedBlocksList)
    
    shortestOutIG = []
    for i in range(len(invertedBlocksList)):
        eNode = [i]
        eNode.extend(list(np.where(invertedBlocksTree.iloc[i] == True)[0]))
        shortestOutIG.append(getAllLongestPaths(outIG,i,eNode,"weight"))
    ## find all shortest paths, i.e. all possible inversions
#    shortestOutIG = []
#    if len(outIG.es) > 0:
#        for i in range(len(invertedBlocksList)):
#            eNode = [i]
#            eNode.extend(list(np.where(invertedBlocksTree.iloc[i] == True)[0]))
#            shortestOutIG.append(getAllLongestPaths(outIG,i,eNode,"weight"))
            
            
#            shortestOutIG = shortestOutIG + [outIG.get_all_shortest_paths(i,[i].extend(list(np.where(invertedBlocksTree.iloc[i] == True)[0])), "weight")]
        
    invTransScores = getTranslocationScore(shortestOutIG, invertedCoords)
    invTransBlocks = getTransBlocks(invTransScores, shortestOutIG, invertedBlocks, inPlaceBlocks, threshold)
    
    
    
    
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
    
    allTransCluster = getTransCluster(allTransGroupIndices, allTransGenomeAGroups, allTransGenomeBGroups)
    
    
    allTransClusterIndices = dict()
    for i in range(len(allTransCluster)):
        allTransClusterIndices.update(dict.fromkeys(allTransCluster[i], i))
    
    allTransBlocksData = []
    for i in range(allTransBlocks.shape[0]):
        tempTransBlock = transBlock(allTransBlocks.iat[i,0],\
                                    allTransBlocks.iat[i,1],\
                                    allTransBlocks.iat[i,2],\
                                    allTransBlocks.iat[i,3],\
                                    allTransBlocks.iat[i,4],\
    #						invTransBlocks[i],\
                                    allTransClusterIndices[i],\
                                    i)
    #	tempTransBlock.addOrderedData(invertedBlocks.iloc[tempTransBlock.orderedBlocksIndex])
        tempTransBlock.addTransGroupIndices(allTransGroupIndices[i])
        tempTransBlock.checkOverlapWithSynBlocks(inPlaceBlocks, threshold)
        tempTransBlock.addGenomeGroupMembers(allTransGenomeAGroups, allTransGenomeBGroups)
        if (tempTransBlock.aUni and tempTransBlock.genomeAUni)	or (tempTransBlock.bUni and tempTransBlock.genomeBUni):
            tempTransBlock.setStatus(1)
        allTransBlocksData.append(tempTransBlock)

    
    for i in range(allTransBlocks.shape[0]):
        tempTransBlock = allTransBlocksData[i]
        if not tempTransBlock.aUni and not tempTransBlock.bUni:
            allTransCluster[allTransClusterIndices[i]].remove(i)
        elif tempTransBlock.status == 1:
            continue
        elif not tempTransBlock.aUni:
            for j in tempTransBlock.genomeBMembers:
                if allTransBlocksData[j].bStart - 50 < tempTransBlock.bStart and allTransBlocksData[j].bEnd + 50 > tempTransBlock.bEnd:
                    tempTransBlock.addMEBlock(j)
        elif not tempTransBlock.bUni:
            for j in tempTransBlock.genomeAMembers:
                if allTransBlocksData[j].aStart - 50 < tempTransBlock.aStart and allTransBlocksData[j].aEnd + 50 > tempTransBlock.aEnd:
                    tempTransBlock.addMEBlock(j)
        else:
            ME_A = []
            for j in tempTransBlock.genomeAMembers:
                if allTransBlocksData[j].aStart - 50 < tempTransBlock.aStart and allTransBlocksData[j].aEnd + 50 > tempTransBlock.aEnd:
                    ME_A.append(j)
            ME_B = []
            for j in tempTransBlock.genomeBMembers:
                if allTransBlocksData[j].bStart - 50 < tempTransBlock.bStart and allTransBlocksData[j].bEnd + 50 > tempTransBlock.bEnd:
                    ME_B.append(j)
            tempTransBlock.setMEList(ME_A, ME_B)

    
    clusterSolutions = []
    for i in range(len(allTransCluster)):
        tempCluster = allTransCluster[i].copy()
        clusterSolutions.append(getBestClusterSubset(tempCluster, allTransBlocksData))
    
    clusterSolutionBlocks = [i[1] for i in clusterSolutions]
    clusterBlocks = unlist(clusterSolutionBlocks)
    
    transClasses = getTransClasses(clusterSolutionBlocks, allTransBlocksData)
    
    dupData = allTransBlocks.iloc[transClasses["duplication"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    invDupData = allTransBlocks.iloc[transClasses["invDuplication"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    TLData = allTransBlocks.iloc[transClasses["translocation"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    invTLData = allTransBlocks.iloc[transClasses["invTranslocation"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])  
    
    dupData = getDupGenome(dupData, allTransBlocksData, transClasses)
    invDupData = getDupGenome(invDupData, allTransBlocksData, transClasses)
    
    
    fout = open(cwdPath+chromo+"_invOut.txt","w")
    tempInvBlocks = []
    for i in bestInvPath:
        invPos = profitable[i-1].invPos
        tempInvBlocks.append([invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]])
        fout.write("\t".join(map(str,["#",invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],"-",invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2],"\n"])))
        for j in invPos:
            fout.write("\t".join(map(str,invertedCoordsOri.iloc[j][:4])))
            fout.write("\n")
    fout.close()
    
    ## Grouping Syn blocks
    allBlocks, outClusters = groupSyn(tempInvBlocks, dupData, invDupData, invTLData, TLData, threshold, synData)
    
########################################################################################################################
    fout = open(cwdPath+chromo+"_synOut.txt","w")
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
    
    fout = open(cwdPath+chromo+"_dupOut.txt","w")
    for i in dupData.index.values:
        fout.write("\t".join(map(str,["#",dupData.at[i,"aStart"],dupData.at[i,"aEnd"],"-",dupData.at[i,"bStart"],dupData.at[i,"bEnd"],"-", dupData.at[i,"dupGenomes"],"\n"])))
        for j in transBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,orderedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################    
    
    fout = open(cwdPath+chromo+"_invDupOut.txt","w")
    for i in invDupData.index.values:
        fout.write("\t".join(map(str,["#",invDupData.at[i,"aStart"],invDupData.at[i,"aEnd"],"-",invDupData.at[i,"bStart"],invDupData.at[i,"bEnd"],"-", invDupData.at[i,"dupGenomes"],"\n"])))
        for j in invTransBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,invertedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################
    
    fout = open(cwdPath+chromo+"_TLOut.txt","w")
    for i in TLData.index.values:
        fout.write("\t".join(map(str,["#",TLData.at[i,"aStart"],TLData.at[i,"aEnd"],"-",TLData.at[i,"bStart"],TLData.at[i,"bEnd"],"\n"])))
        for j in transBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,orderedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################
    
    fout = open(cwdPath+chromo+"_invTLOut.txt","w")
    for i in invTLData.index.values:
        fout.write("\t".join(map(str,["#",invTLData.at[i,"aStart"],invTLData.at[i,"aEnd"],"-",invTLData.at[i,"bStart"],invTLData.at[i,"bEnd"],"\n"])))
        for j in invTransBlocks[allTransIndexOrder[i]]:
            fout.write("\t".join(map(str,invertedBlocks.iloc[j][:4])))
            fout.write("\n")
    fout.close()

########################################################################################################################
def getCTX(coords, cwdPath, uniChromo):
    annoCoords = readAnnoCoords(cwdPath, uniChromo)
    ctxData = coords.loc[coords['aChr'] != coords['bChr']].copy()
#    ctxData = ctxBlocks.copy()
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
#    ctxBlocks = ctxData.copy()
    
    orderedBlocks = ctxData[ctxData.bDir == 1]
    invertedBlocks = ctxData[ctxData.bDir == -1]
    
    ## Create connectivity tree for directed blocks
#    transBlocksNeighbours = getTransSynOrientation(annoCoords, orderedBlocks, threshold, ctx = True)
    outOrderedBlocks = makeBlocksTree(annoCoords, orderedBlocks, threshold, ctx = True)

#    invTransBlocksNeighbours = getTransSynOrientation(annoCoords, invertedBlocks, threshold, ctx = True)
       
    outInvertedBlocks = makeBlocksTree(annoCoords, invertedBlocks, threshold, ctx = True)
 
    
"""
DEBUG FROM HERE
"""
    
    
    
    ## find all translocations which don't have large gaps between its alignments
    ## and are not overlappign with the syntenic blocks
    orderedBlocksTree = outOrderedBlocks   
    orderedBlocksList = makeBlocksList(orderedBlocksTree, orderedBlocks)
    outOG = getConnectivityGraph(orderedBlocksList)
    
    ## find all shortest paths, i.e. all possible inversions
    shortestOutOG = []
    if len(outOG.es) > 0:
        for i in range(len(orderedBlocksList)):
            eNode = [i]
            eNode.extend(list(np.where(orderedBlocksTree.iloc[i] == True)[0]))
            shortestOutOG.append(getAllLongestPaths(outOG,i,eNode,"weight"))
            
#            shortestOutOG = shortestOutOG + [outOG.get_all_shortest_paths(i,[i].extend(list(np.where(orderedBlocksTree.iloc[i] == True)[0])), "weight")]
        
    transScores = getTranslocationScore(shortestOutOG, orderedBlocks, ctx = True)
    transBlocks = getTransBlocks(transScores, shortestOutOG, orderedBlocks, annoCoords, threshold,ctx = True)
    

    invertedBlocksTree = outInvertedBlocks
    invertedBlocksList = makeBlocksList(invertedBlocksTree, invertedBlocks)
    outIG = getConnectivityGraph(invertedBlocksList)
    
    ## find all shortest paths, i.e. all possible inversions
    shortestOutIG = []
    if len(outIG.es) > 0:
        for i in range(len(invertedBlocksList)):
            eNode = [i]
            eNode.extend(list(np.where(invertedBlocksTree.iloc[i] == True)[0]))
            shortestOutIG.append(getAllLongestPaths(outIG,i,eNode,"weight"))

#            shortestOutIG = shortestOutIG + [outIG.get_all_shortest_paths(i,[i].extend(list(np.where(invertedBlocksTree.iloc[i] == True)[0])), "weight")]
        
    invTransScores = getTranslocationScore(shortestOutIG, invertedBlocks, ctx = True)
    invTransBlocks = getTransBlocks(invTransScores, shortestOutIG, invertedBlocks, annoCoords, threshold, ctx = True)

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
    
    ctxBlocksData = []
    for i in ctxTransBlocks.index.values:
        tempTransBlock = transBlock(ctxTransBlocks.at[i,"aStart"],\
                                    ctxTransBlocks.at[i,"aEnd"],\
                                    ctxTransBlocks.at[i,"bStart"],\
                                    ctxTransBlocks.at[i,"bEnd"],\
                                    ctxTransBlocks.at[i,"bDir"],\
    #						invTransBlocks[i],\
                                    ctxClusterIndices[i],\
                                    i)
    #	tempTransBlock.addOrderedData(invertedBlocks.iloc[tempTransBlock.orderedBlocksIndex])
        tempTransBlock.addTransGroupIndices(ctxGroupIndices[i])
        tempTransBlock.checkOverlapWithSynBlocks_A(annoCoords.loc[annoCoords.aChr == ctxTransBlocks.at[i,"aChr"]], threshold)
        tempTransBlock.checkOverlapWithSynBlocks_B(annoCoords.loc[annoCoords.bChr == ctxTransBlocks.at[i,"bChr"]], threshold)
        tempTransBlock.addGenomeGroupMembers(ctxTransGenomeAGroups, ctxTransGenomeBGroups)
        if (tempTransBlock.aUni and tempTransBlock.genomeAUni) or (tempTransBlock.bUni and tempTransBlock.genomeBUni):
            tempTransBlock.setStatus(1)
        ctxBlocksData.append(tempTransBlock)

    
    for i in range(len(ctxBlocksData)):
        tempTransBlock = ctxBlocksData[i]
        index = tempTransBlock.transBlocksID
        if not tempTransBlock.aUni and not tempTransBlock.bUni:
            ctxCluster[ctxClusterIndices[index]].remove(index)
        elif tempTransBlock.status == 1:
            continue
        elif not tempTransBlock.aUni:
            for j in tempTransBlock.genomeBMembers:
                j_pos = np.where(ctxTransBlocks.index.values == j)[0][0]
                if ctxBlocksData[j_pos].bStart - 50 < tempTransBlock.bStart and ctxBlocksData[j_pos].bEnd + 50 > tempTransBlock.bEnd:
                    tempTransBlock.addMEBlock(j)
        elif not tempTransBlock.bUni:
            for j in tempTransBlock.genomeAMembers:
                j_pos = np.where(ctxTransBlocks.index.values == j)[0][0]
                if ctxBlocksData[j_pos].aStart - 50 < tempTransBlock.aStart and ctxBlocksData[j_pos].aEnd + 50 > tempTransBlock.aEnd:
                    tempTransBlock.addMEBlock(j)
        else:
            ME_A = []
            for j in tempTransBlock.genomeAMembers:
                j_pos = np.where(ctxTransBlocks.index.values == j)[0][0]
                if ctxBlocksData[j_pos].aStart - 50 < tempTransBlock.aStart and ctxBlocksData[j_pos].aEnd + 50 > tempTransBlock.aEnd:
                    ME_A.append(j)
            ME_B = []
            for j in tempTransBlock.genomeBMembers:
                j_pos = np.where(ctxTransBlocks.index.values == j)[0][0]
                if ctxBlocksData[j_pos].bStart - 50 < tempTransBlock.bStart and ctxBlocksData[j_pos].bEnd + 50 > tempTransBlock.bEnd:
                    ME_B.append(j)
            tempTransBlock.setMEList(ME_A, ME_B)

    
    clusterSolutions = []
    for i in range(len(ctxCluster)):
#    for i in range(10):
        print(i)
        tempCluster = ctxCluster[i].copy()
        if len(tempCluster) == 0:
            continue
        else:
            clusterSolutions.append(getBestClusterSubset(tempCluster, ctxBlocksData))
        
    clusterSolutionBlocks = [i[1] for i in clusterSolutions]
    clusterBlocks = unlist(clusterSolutionBlocks)
    
    transClasses = getTransClasses(clusterSolutionBlocks, ctxBlocksData)
    
    indices = sorted(unlist(list(transClasses.values())))
    keys = [key for index in indices for key in list(transClasses.keys()) if index in transClasses[key]]
    blocksClasses = dict(zip(indices,keys))
    
    fout = open(cwdPath+"ctxOut.txt","w")
    
    for index in indices:
        if ctxBlocksData[index].dir == 1:
            alignIndices = transBlocks[ctxTransIndexOrder[index]]
            fout.write("#\t" + "\t".join(map(str,ctxTransBlocks.iloc[index])) + "\t" + blocksClasses[index]+ "\n")
            for i in alignIndices:
                fout.write("\t".join(map(str,orderedBlocks.iloc[i,0:-2]))+"\n")            
        elif ctxBlocksData[index].dir == -1:
            alignIndices = invTransBlocks[ctxTransIndexOrder[index]]
            fout.write("#\t" + "\t".join(map(str,ctxTransBlocks.iloc[index]))+ "\t" + blocksClasses[index] +"\n")
            for i in alignIndices:
                fout.write("\t".join(map(str,invertedBlocks.iloc[i,0:-2])) + "\n")
    fout.close()
    return 0
    
   
if __name__ == "__main__":
    if len(sys.argv) not in [2,3]:
        sys.exit("Usage: SynSearch <path_to_coords_file> <number_of_CPU_cores>")
    if len(sys.argv) == 2:
        fileLocation = sys.argv[1]   
        nCores = 1
    elif len(sys.argv) == 3:
        fileLocation = sys.argv[1]
        nCores = int(sys.argv[2])
        
#    chromo = sys.argv[2]
#    cwdPath = sys.argv[3]
    cwdPath = os.getcwd()+ os.sep
    coords = pd.read_table(fileLocation, header = None) 
    coords.columns  =   ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr"]
    aChromo = set(coords["aChr"])
    bChromo = set(coords["bChr"])
    uniChromo = list(aChromo if len(aChromo) < len(bChromo) else bChromo)
    uniChromo.sort()
    # Set this as an argument
    ## This is the number of non-overlapping BPs required to consider two overlapping
    ## alignments as different
    threshold = 50
    procs = []
    print(uniChromo)
    with Pool(processes = nCores) as pool:
        pool.map(partial(SynSearch,threshold=threshold,coords=coords, cwdPath= cwdPath), uniChromo) 
    mergeOutputFiles(uniChromo,cwdPath)
#    ctxBlocks = getCTX(coords, cwdPath, uniChromo)
    
    
    
    
#    ctxBlocks.to_csv(cwdPath+"ctxOut.txt", sep="\t",index=False) 
        
#    for chromo in uniChromo:
#        procs.append(multiprocessing.Process(target=SynSearch, args=(chromo,threshold,coords,cwdPath,)))
##    SynSearch(chromo, threshold, coords)
#    for p in procs:
#        p.start()
#    for p in procs:
#        p.join()
    print("End of program")
