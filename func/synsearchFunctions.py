# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:54:53 2017

@author: goel
"""
import numpy as np
from func.myUsefulFunctions import *
import sys
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from math import log
import matplotlib.patheffects as pe
from itertools import chain
import time
from operator import itemgetter
from igraph import *
from collections import Counter
from scipy.stats import *
import datetime


np.random.seed(1)


def SynAndOverlap(jData, iData, threshold):
    """Compute whether two alignments are syntenic to each other.
       
       Parameters
       ----------
       jData: pandas Series, 
           Series containing the information of the second (right) alignment
           Required elements: `aStart, aEnd, bStart, bEnd, aDir, bDir`
       iData: pandas Series,
           Series containing the information of the first (left) alignment
           Required elements: `aStart, aEnd, bStart, bEnd, aDir, bDir`
       threshold: int, 
           cut-off value.
      
       Returns
       --------
       Binary value
           True if two alignments can be syntenic to each other, False otherwise
    """
    if jData.bDir == iData.bDir:
        if (jData.aStart - iData.aStart) > threshold and (jData.aEnd - iData.aEnd) > threshold and (jData.bStart - iData.bStart) > threshold and (jData.bEnd - iData.bEnd) > threshold:
            return True
    else:
        if (jData.aStart - iData.aStart) > threshold and (jData.aEnd - iData.aEnd) > threshold and (jData.bStart - iData.bEnd) > threshold and (jData.bEnd - iData.bStart) > threshold:
            return True
    return False

def TS(jData, iData, threshold):
    if (jData.aStart - iData.aStart) > threshold and (jData.aEnd - iData.aEnd) > threshold and (jData.bStart - iData.bStart) > threshold and (jData.bEnd - iData.bEnd) > threshold:
            return True
    return False

def getSynPath(blocks):
    synPath = []
    scores = [block.score for block in blocks]
    #bestScore = max(scores)
    lastBlock = scores.index(max(scores))
    while blocks[lastBlock].bestParentID != -1:
        synPath.append(lastBlock)
        lastBlock = blocks[lastBlock].bestParentID        
    synPath.append(lastBlock)
    return(synPath[::-1])
    
def getInversions(coords,chromo, threshold, synData, synPath):
    
    class inversion:
        def __init__(self, cost, revenue, neighbours, invPos):
            self.cost = cost
            self.revenue = revenue
            self.profit = revenue - cost
            self.neighbours = list(neighbours)
            self.invPos = invPos
    
    def bestInv(profG):
        vcount = profG.vcount()
        for i in range(vcount):
            parents = profG.neighbors(i,IN)
            if len(parents) > 0:
                scores = [profG.vs[j]["score"] for j in parents]
                profG.vs[i]["parent"] = parents[scores.index(max(scores))]
                profG.vs[i]["score"] = profG.vs[i]["score"] + max(scores)
        nodeScores = profG.vs["score"]
    #    print(nodeScores)
        node = profG.vs[nodeScores.index(max(nodeScores))]
        bestInvPath = []
        while(node["parent"]):
            bestInvPath.append(node.index)
            node = profG.vs[node["parent"]]
    
        bestInvPath.append(node.index)
        return(bestInvPath[::-1])

    def getNeighbours(neighbourSyn, j):
        return(min(neighbourSyn[j[0]]+neighbourSyn[j[-1]]), max(neighbourSyn[j[0]]+neighbourSyn[j[-1]]))
    
    invertedCoordsOri = coords.loc[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == -1)]
    invertedCoords = invertedCoordsOri.copy()    
#    minCoords = np.min(np.min(invertedCoords[["bStart","bEnd"]]))
    maxCoords = np.max(np.max(invertedCoords[["bStart","bEnd"]]))
    
    invertedCoords.bStart = maxCoords + 1 - invertedCoords.bStart 
    invertedCoords.bEnd = maxCoords + 1 - invertedCoords.bEnd

    nrow = pd.Series(range(invertedCoords.shape[0]))
    invTree = nrow.apply( lambda x : invertedCoords.iloc[x:,].apply(TS, axis = 1, args = (invertedCoords.iloc[x], threshold)))
    del(invertedCoords)
    
    #######################################################################
    ###### Create list of inverted alignments
    #######################################################################

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
    
    #########################################################################
    ###### Finding profitable inversions (group of inverted blocks)
    #########################################################################
#    print("finding possible inversions")
    
    shortest = []
    invG = getConnectivityGraph(invBlocks)        
    if len(invG.es) > 0:        
        for i in range(len(invBlocks)):
            shortest.append(getAllLongestPaths(invG,i,range(len(invBlocks))))
            
    ## Get revenue of shortest paths, score of adding the inversion
    
    ##### NEED TO CHANGE THIS TO GIVE LOWER SCORE TO OVERLAPPING INVERSIONS
    ## THIS MAYBE CONTAIN SOME BUGS!!!
    revenue = []
    for i in shortest:
        values = []
        for j in i:
            if len(j) == 1:
                values.append(invBlocks[j[0]].score)
            else:              
                score = 0
                startA = [invertedCoordsOri.iat[j[0],0]]
                endA = [invertedCoordsOri.iat[j[0],1]]
                startB = [invertedCoordsOri.iat[j[0],3]]
                endB = [invertedCoordsOri.iat[j[0],2]]
                iden = [invertedCoordsOri.iat[j[0],6]]
                for k in j[1:]:
                    isMore = True if invertedCoordsOri.iat[k,6] > iden[-1] else False
                    if invertedCoordsOri.iat[k,0] < endA[-1]:
                        if isMore:
                            endA[-1] = invertedCoordsOri.iat[k,0]
                            startA.append(invertedCoordsOri.iat[k,0])
                            endA.append(invertedCoordsOri.iat[k,1])
                        else:
                            startA.append(endA[-1])
                            endA.append(invertedCoordsOri.iat[k,1])
                    else:
                        startA.append(invertedCoordsOri.iat[k,0])
                        endA.append(invertedCoordsOri.iat[k,1])
                    
                    if invertedCoordsOri.iat[k,2] > startB[-1]:
                        if isMore:
                            startB[-1] = invertedCoordsOri.iat[k,2]
                            startB.append(invertedCoordsOri.iat[k,3])
                            endB.append(invertedCoordsOri.iat[k,2])
                        else:
                            endB.append(startB[-1])
                            startB.append(invertedCoordsOri.iat[k,3])
                    else:
                        startB.append(invertedCoordsOri.iat[k,3])
                        endB.append(invertedCoordsOri.iat[k,2])
                    iden.append(invertedCoordsOri.iat[k,6])
                if len(startA) == len(endA) == len(startB) == len(endB) == len(iden):
                    for i in range(len(iden)):
                        score += iden[i]*((endA[i] - startA[i]) + (endB[i] - startB[i]))
                values.append(score)
        revenue = revenue + [values]
    
    ## Get syntenic neighbouring blocks of inversions
    neighbourSyn = dict()
    for i in range(invertedCoordsOri.shape[0]):
        index = invertedCoordsOri.index.values[i]
        upSyn = np.where(synData.index.values < index)[0]
        downSyn = np.where(synData.index.values > index)[0]
        
        upBlock  = -1
        downBlock = len(synData)    
        for j in upSyn[::-1]:
            if SynAndOverlap(invertedCoordsOri.loc[index], synData.iloc[j], threshold):
                upBlock = j
                break
        
        for j in downSyn:
            if SynAndOverlap(synData.iloc[j], invertedCoordsOri.loc[index], threshold):
                downBlock = j
                break
        neighbourSyn[i] = [upBlock, downBlock]
    
    ## Calculate score of individual synblock
    synBlockScore = [(i.aLen + i.bLen)*i.iden for index,i in synData.iterrows()]
    
    ## Calculate cost adding an inversion, i.e sum of all synblocks which need to be removed to accomodate teh synblocks
    
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
                if overlapLength < ((rightEnd - leftEnd)/2):
                    values.append(synCost + sys.maxsize)
                else:
                    values.append(synCost)
        cost = cost + [values]
    ## Calculate profit (or loss) associated with the addition of an inversion
    profit = []
    for i in range(len(revenue)):
        profit = profit + [[revenue[i][j] - cost[i][j]for j in range(len(revenue[i]))]]
    
    ## Create list of all profitable inversions
    
    ##invPos are 0-indexed positions of inverted alignments in the invertedCoordsOri object
    profitable = [inversion(cost[i][j], revenue[i][j],
                             getNeighbours(neighbourSyn, shortest[i][j]),shortest[i][j])
                             for i in range(len(profit)) for j in range(len(profit[i]))\
                                 if profit[i][j] > (0.1*cost[i][j])]     ##Select only those inversions for which the profit is more than  10% of the cost
    #####################################################################
    #### Find optimal set of inversions from all profitable inversions
    #####################################################################
    profitInvs = [p.profit for p in profitable]
    if len(profitInvs) > 0:  
        vcount = len(profitable)+2
        profG =  Graph().as_directed()
        profG.add_vertices(vcount)
        iAStart = []
        iAEnd = []
        iBStart = []
        iBEnd = []
        
        for i in profitable:
            iAStart.append(invertedCoordsOri.iat[i.invPos[0], 0])
            iAEnd.append(invertedCoordsOri.iat[i.invPos[-1], 1])
            iBStart.append(invertedCoordsOri.iat[i.invPos[-1], 3])
            iBEnd.append(invertedCoordsOri.iat[i.invPos[0], 2])
        
        nonOverLapA = [np.where(iAStart > i - threshold)[0] for i in iAEnd] 
        nonOverLapB = [np.where(iBStart > i - threshold)[0] for i in iBEnd]
        for i in range(len(profitable)):
            childNodes = np.intersect1d(nonOverLapA[i], nonOverLapB[i]) + 1               ## two inversions can co-exist only if the overlap between them is less than threshold on both genomes 
            profG.add_edges(zip([i+1]*len(childNodes), childNodes))
            profG.es[-len(childNodes):]["weight"] = profitable[i].profit
            profG.vs[i+1]["child"] = list(childNodes)
            profG.vs[i+1]["score"] = profitable[i].profit
    
        profG.vs[[0,vcount-1]]["score"] = 0
        
        noInEdge = np.where(np.array(profG.vs[1:-1].indegree()) == 0)[0] + 1
        noOutEdge = np.where(np.array(profG.vs[1:-1].outdegree()) == 0)[0] + 1
        
        for i in noInEdge:
            profG.add_edge(0,i)
            profG.es[-1:]["weight"] = 0
        profG.vs[0]["child"] = noInEdge
        for i in noOutEdge:
            profG.add_edge(i,vcount - 1)
            profG.es[-1:]["weight"] = profitable[i-1].profit         
            profG.vs[i]["child"] = vcount-1
        
                
        bestInvPath = bestInv(profG)
    else:
        bestInvPath = []

    invBlocksIndex = unlist([profitable[i-1].invPos for i in bestInvPath])
    invData = invertedCoordsOri.iloc[invBlocksIndex]
    invNeighbours  = [profitable[i-1].neighbours for i in bestInvPath]
    synInInv = unlist([list(range(i[0]+1, i[1])) for i in invNeighbours])
#    synInInvIndices = getValues(synData.index, synInInv)
    return(invertedCoordsOri, profitable, bestInvPath,invData, synInInv)
        
           
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
        return transPositions

    if ctx:
        
        inPlaceBlocks.sort_values(["aIndex"], inplace = True)
        transRowIndex = transData.index.values
        transPositions = dict()
        ## identify neighbours on reference genome
        for i in transRowIndex:
            upSyn = getValues(inPlaceBlocks.index.values,
                              np.intersect1d(
                                      np.where(inPlaceBlocks.aStart <= transData.at[i,"aEnd"])[0],
                                      np.where(inPlaceBlocks.aChr == transData.at[i,"aChr"])[0]))
            downSyn = getValues(inPlaceBlocks.index.values,
                                np.intersect1d(
                                        np.where(inPlaceBlocks.aEnd >= transData.at[i,"aStart"])[0],
                                        np.where(inPlaceBlocks.aChr == transData.at[i,"aChr"])[0]))
            
            upBlock = -1
            downBlock = len(inPlaceBlocks)
            
            for j in upSyn[::-1]:
                if (transData.at[i,"aStart"] - inPlaceBlocks.at[j,"aStart"]) > threshold and (transData.at[i,"aEnd"] - inPlaceBlocks.at[j,"aEnd"]) > threshold:
                    upBlock = j
                    break
            
            for j in downSyn:
                if (inPlaceBlocks.at[j,"aStart"] - transData.at[i,"aStart"]) > threshold and (inPlaceBlocks.at[j,"aEnd"] - transData.at[i,"aEnd"]) > threshold:
                    downBlock = j 
                    break
            
            transPositions[i] = [inPlaceBlocks.at[upBlock,"aIndex"]] if upBlock != -1 else [-1]
            transPositions[i].append(inPlaceBlocks.at[downBlock,"aIndex"]) if downBlock != len(inPlaceBlocks) else transPositions[i].append(len(inPlaceBlocks))
            
        
        ## identify neighbours on query genome
        inPlaceBlocks.sort_values("bIndex", inplace = True)
        for i in transRowIndex:
            bUpSyn = getValues(inPlaceBlocks.index.values,
                               np.intersect1d(
                                       np.where(inPlaceBlocks.bStart <= transData.at[i,"bEnd"])[0],
                                       np.where(inPlaceBlocks.bChr == transData.at[i,"bChr"])[0]))
            bDownSyn = getValues(inPlaceBlocks.index.values,
                                 np.intersect1d(
                                         np.where(inPlaceBlocks.bEnd >= transData.at[i,"bStart"])[0],
                                         np.where(inPlaceBlocks.bChr == transData.at[i, "bChr"])[0]))

            bUpBlock = -1
            bDownBlock = len(inPlaceBlocks)
            
            for j in bUpSyn[::-1]:
                if (transData.at[i,"bStart"] - inPlaceBlocks.at[j,"bStart"]) > threshold and (transData.at[i,"bEnd"] - inPlaceBlocks.at[j,"bEnd"]) > threshold:
                    bUpBlock = j
                    break
            
            for j in bDownSyn:
                if (inPlaceBlocks.at[j,"bStart"] - transData.at[i,"bStart"]) > threshold and (inPlaceBlocks.at[j,"bEnd"] - transData.at[i,"bEnd"]) > threshold:
                    bDownBlock = j 
                    break 
            transPositions[i].append(inPlaceBlocks.at[bUpBlock,"bIndex"]) if bUpBlock != -1 else transPositions[i].append(-1)
            transPositions[i].append(inPlaceBlocks.at[bDownBlock,"bIndex"]) if bDownBlock != len(inPlaceBlocks) else transPositions[i].append(len(inPlaceBlocks))
        return transPositions
    
        

#%%

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
        return(transBlocksData, orderedIndex)
        

def findOverlappingSynBlocks(inPlaceBlocks, aStart, aEnd, bStart, bEnd):
    aBlocks = list(np.intersect1d(np.where(inPlaceBlocks.aStart.values < aEnd)[0],\
                                      np.where(inPlaceBlocks.aEnd.values > aStart)[0]))
    bBlocks = list(np.intersect1d(np.where(inPlaceBlocks.bStart.values < bEnd)[0],\
                                      np.where(inPlaceBlocks.bEnd.values > bStart)[0]))
    return(aBlocks, bBlocks)
    

def getConnectivityGraph(blocksList):
    outOG = Graph().as_directed()
    outOG.add_vertices(len(blocksList))
    if len(blocksList) == 0:
        return outOG
    ## Add edges and edge weight
    for i in blocksList:
        if len(i.children) > 0:
            outOG.add_edges(zip([i.id]*len(i.children), i.children))
            outOG.es[-len(i.children):]["weight"] = [-i.score]*len(i.children)
    return outOG

def getAllLongestPaths(graph,sNode, eNode,by="weight"):
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
    pathList = []
    dist = {}
    pred = {}
    
    allNodes = graph.vs.indices
    for i in allNodes:
        dist[i] = float("inf")
        pred[i] = None
        
    dist[sNode] = 0
    changes = 1
    
    for i in range(len(allNodes)-1):
        if changes == 0:
            break
        changes = 0
        for e in graph.es:
            if dist[e.source] + e[by] < dist[e.target]:
                changes = 1
                dist[e.target] = dist[e.source] + e[by]
                pred[e.target] = e.source
    
    for e in graph.es:
        if dist[e.source] + e[by] < dist[e.target]:
            sys.exit("Negative weight cycle identified")
    
    for key in eNode:
        if dist[key] != float("inf"):
            path = []
            while key!=sNode:
                path.append(key)
                key = pred[key]
            path.append(sNode)
            pathList.append(path[::-1])
    return(pathList)



def findOrderedTranslocations(outOrderedBlocks, orderedBlocks, inPlaceBlocks, threshold, ctx = False):
    if not isinstance(ctx, bool):
        print("CTX status must be a boolean")
        sys.exit()
    def makeBlocksList(blocksTree, blocksData):
        nrow = blocksTree.shape[0]
        blocksList = [alingmentBlock(i, np.where(blocksTree.iloc[i] == True)[0],blocksData.iloc[i]) for i in range(nrow)]
        for block in blocksList:
            i = 0
            while(i < len(block.children)):
                block.children = list(set(block.children) - set(blocksList[block.children[i]].children))
                i+=1
            block.children.sort()
            
            for child in block.children:
                blocksList[child].addParent(block.id)
        return blocksList
        
    def getTranslocationScore(translocations, transData, ctx):
        """Function to score the proposed translocation block based on the number of
            basepairs it explains and the gaps between alignments of the block
        """
        if not isinstance(ctx, bool):
            print("CTX status must be a boolean")
            sys.exit()
        if not ctx:
            transScores = []
            for blocks in translocations:
                blocksScores = []
                for block in blocks:
                    aScore = transData.iat[block[0],4]
                    bScore = transData.iat[block[0],5]
                    aGap = 0
                    bGap = 0
                    if len(block) > 1:
                        for i in range(1, len(block)):
                            aScore += transData.iat[block[i],4]
                            bScore += transData.iat[block[i],5]
                            aGap += max(0,transData.iat[block[i],0] - transData.iat[block[i-1],1])
                            bGap += max(0,transData.iat[block[i],2] - transData.iat[block[i-1],3])
                        blockScore = min(((aScore - aGap)/aScore),((bScore - bGap)/bScore))
                        blocksScores.append(blockScore)
                    else:
                        blocksScores.append(1)
                transScores.append(blocksScores)
            return transScores
        if ctx:
            transScores = []
            for blocks in translocations:
        #        print(len(blocks))
                blocksScores = []
                for block in blocks:
                    indices = getValues(transData.index.values,block)
                    aScore = transData.at[indices[0],"aLen"]
                    bScore = transData.at[indices[0],"bLen"]
                    aGap = 0
                    bGap = 0
                    if len(block) > 1:
                        for i in range(1, len(block)):
                            aScore += transData.at[indices[i],"aLen"]
                            bScore += transData.at[indices[i],"bLen"]
                            aGap += max(0,transData.at[indices[i],"aStart"] - transData.at[indices[i-1],"aEnd"])
                            bGap += max(0,transData.at[indices[i],"bStart"] - transData.at[indices[i-1],"bEnd"]) if transData.at[indices[i],"bDir"] == 1 else max(0, transData.at[indices[i-1],"bStart"] - transData.at[indices[i],"bEnd"]) 
                        blockScore = min(((aScore - aGap)/aScore),((bScore - bGap)/bScore))
                        blocksScores.append(blockScore)
                    else:
                        blocksScores.append(1)
                transScores.append(blocksScores)
            return transScores
    
    def getTransBlocks(transScores, shortTrans, transData, inPlaceBlocks, threshold, ctx):
        """This method filters possible translocation blocks to select those which have a posivitive gap based score
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
        transBlocks = [getValues(shortTrans[i], positiveTransScores[i]) for i in range(len(shortTrans))]
        transBlocks = [i for j in transBlocks for i in j]
        outBlocks = []
        if not isinstance(ctx,bool):
            print("CTX status must be a boolean")
            sys.exit()
        if not ctx:
            for block in transBlocks:
                blockAlength = 0
                blockBlength = 0
                blockAUni = 0
                blockBUni = 0
                
                for almnt in block:
                    aStart = transData.iat[almnt,0]
                    aEnd = transData.iat[almnt,1]
                    if transData.iat[almnt,8] == 1:
                        bStart = transData.iat[almnt,2]
                        bEnd = transData.iat[almnt,3]
                    else:
                        bStart = transData.iat[almnt,3]
                        bEnd = transData.iat[almnt,2]      
                    blockAlength += transData.iat[almnt,4]
                    blockBlength += transData.iat[almnt,5]
                    
                    aBlocks, bBlocks = findOverlappingSynBlocks(inPlaceBlocks, aStart, aEnd, bStart, bEnd)
        
                    for aBlock in aBlocks:
                        if inPlaceBlocks.iat[aBlock,0] - aStart < threshold and aEnd - inPlaceBlocks.iat[aBlock,1] < threshold:
                            aStart = aEnd
                            break
                        elif inPlaceBlocks.iat[aBlock,0] < aStart and inPlaceBlocks.iat[aBlock,1] < aEnd:
                            aStart = inPlaceBlocks.iat[aBlock,1]
                        else:
                            blockAUni += inPlaceBlocks.iat[aBlock,0] - aStart
                            if inPlaceBlocks.iat[aBlock,1] < aEnd:
                                aStart = inPlaceBlocks.iat[aBlock,1]
                            else:
                                aStart = aEnd
                                break
                    blockAUni += aEnd - aStart
                    
                    
                    for bBlock in bBlocks:
                        bBlockEnd = inPlaceBlocks.iat[bBlock,3]
                        bBlockStart = inPlaceBlocks.iat[bBlock,2]
                        if bBlockStart - bStart < threshold and bEnd - bBlockEnd < threshold:
                            bStart = bEnd
                            break
                        elif bBlockStart < bStart and bBlockEnd < bEnd:
                            bStart = bBlockEnd
                        else:
                            blockBUni += bBlockStart - bStart
                            if bBlockEnd < bEnd:
                                bStart = bBlockEnd
                            else:
                                bStart = bEnd
                                break
                    blockBUni += bEnd - bStart
        #Trans block is selected IFF either the unique region on any genome is larger than 1kb
        # or length of unique region on a genome is larger than 0.5 times the length of
        # the overlapping region on that genome
                if blockAUni > 1000 or blockBUni > 1000 or blockAUni > 0.5*blockAlength or blockBUni > 0.5*blockBlength:
                    outBlocks.append(block)
            return(outBlocks)
        ##########
        ## With CTX
        ##########
        if ctx:
            for block in transBlocks:
                blockAlength = 0
                blockBlength = 0
                blockAUni = 0
                blockBUni = 0
                
                for almnt in block:
                    index = transData.index.values[almnt]
                    aStart = transData.at[index,"aStart"]
                    aEnd = transData.at[index,"aEnd"]
                    bStart = transData.at[index,"bStart"]
                    bEnd = transData.at[index,"bEnd"]
                    
                    if bEnd < bStart:
                        print("CTX Input: bStart must be less than bEnd")
                        sys.exit()
     
                    blockAlength += transData.at[index,"aLen"]
                    blockBlength += transData.at[index,"bLen"]
                    
    #                aBlocks, bBlocks = findOverlappingSynBlocks(inPlaceBlocks, aStart, aEnd, bStart, bEnd)
                    aBlocks = list(np.intersect1d(np.where(inPlaceBlocks.aStart.values < aEnd)[0], np.where(inPlaceBlocks.aEnd.values > aStart)[0]))
                    aBlocks = list(np.intersect1d(aBlocks, np.where(inPlaceBlocks.aChr == transData.at[index,"aChr"])[0]))
                    aBlocks = getValues(inPlaceBlocks.index.values,aBlocks)
        
                    for aBlock in aBlocks:
                        if inPlaceBlocks.at[aBlock,"aStart"] - aStart < threshold and aEnd - inPlaceBlocks.at[aBlock,"aEnd"] < threshold:
                            aStart = aEnd
                            break
                        elif inPlaceBlocks.at[aBlock,"aStart"] < aStart and inPlaceBlocks.at[aBlock,"aEnd"] < aEnd:
                            aStart = inPlaceBlocks.at[aBlock,"aEnd"]
                        else:
                            blockAUni += inPlaceBlocks.at[aBlock,"aStart"] - aStart
                            if inPlaceBlocks.at[aBlock,"aEnd"] < aEnd:
                                aStart = inPlaceBlocks.at[aBlock,"aEnd"]
                            else:
                                aStart = aEnd
                                break
                    blockAUni += aEnd - aStart
                    
                    bBlocks = list(np.intersect1d(np.where(inPlaceBlocks.bStart.values < bEnd)[0], np.where(inPlaceBlocks.bEnd.values > bStart)[0]))
                    bBlocks = list(np.intersect1d(bBlocks, np.where(inPlaceBlocks.bChr == transData.at[index,"bChr"])[0]))
                    bBlocks = getValues(inPlaceBlocks.index.values, bBlocks)
                    
                    for bBlock in bBlocks:
                        bBlockStart = inPlaceBlocks.at[bBlock,"bStart"]
                        bBlockEnd = inPlaceBlocks.at[bBlock,"bEnd"]
                        if bStart - bBlockStart < threshold and bEnd - bBlockEnd < threshold:
                            bStart = bEnd
                            break
                        elif bBlockStart < bStart and bBlockEnd < bEnd:
                            bStart = bBlockEnd
                        else:
                            blockBUni += bBlockStart - bStart
                            if bBlockEnd < bEnd:
                                bStart = bBlockEnd
                            else:
                                bStart = bEnd
                                break
                    blockBUni += bEnd - bStart
        #Trans block is selected IFF either the unique region on any genome is larger than 1kb
        # or length of unique region on a genome is larger than 0.5 times the length of
        # the overlapping region on that genome
                if blockAUni > 1000 or blockBUni > 1000 or blockAUni > 0.5*blockAlength or blockBUni > 0.5*blockBlength:
                    outBlocks.append(block)
            return(outBlocks)
    
    orderedBlocksList = makeBlocksList(outOrderedBlocks, orderedBlocks)
    outOG = getConnectivityGraph(orderedBlocksList)
    shortestOutOG = []
    for i in range(len(orderedBlocksList)):
        eNode = [i]
        eNode.extend(list(np.where(outOrderedBlocks.iloc[i] == True)[0]))
        shortestOutOG.append(getAllLongestPaths(outOG,i,eNode,"weight"))      
    transScores = getTranslocationScore(shortestOutOG, orderedBlocks, ctx)
    transBlocks = getTransBlocks(transScores, shortestOutOG, orderedBlocks, inPlaceBlocks, threshold, ctx)
    return(transBlocks)
        
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
    return(genomeAGroups, genomeBGroups)


def makeTransGroupList(transBlocksData, start, end, threshold):
    transBlocksTable = transBlocksData.sort_values([start,end])
    indices = transBlocksTable.index.values
    if len(transBlocksData) > 0:
        genomeGroups = [transGroups(transBlocksTable.at[indices[0],start],\
                                    transBlocksTable.at[indices[0],end], indices[0], threshold)]
        for i in indices[1:]:
            if transBlocksTable.at[i, start] > genomeGroups[-1].rightEnd:
                genomeGroups.append(transGroups(transBlocksTable.at[i,start],\
                                                transBlocksTable.at[i,end], i, threshold))
            elif genomeGroups[-1].checkOverlap(transBlocksTable.at[i,start],\
                             transBlocksTable.at[i,end]):
                genomeGroups[-1].addMember(transBlocksTable.at[i,start],\
                            transBlocksTable.at[i,end], i)
            else:
                genomeGroups.append(transGroups(transBlocksTable.at[i,start],\
                                                transBlocksTable.at[i,end], i, threshold))
        return genomeGroups
    else:
        return []
#%%

def makeBlocksTree(inPlaceBlocks, blocksData, threshold, transBlocksNeighbours = None, ctx = False):
    """Compute whether two alignments can be part of one translation block. For this:
        the alignments should not be separated by any inPlaceBlock on both ends and
        they should be syntenic with respect to each other.
       
       Parameters
       ----------
       inPlaceBlocks: pandas DataFrame, 
           dataframe containing coordinates and other information for all inPlaceBlocks
       blocksData: pandas DataFrame,
           dataframe containing coordinates of alignments which are needed to be connected.
           Requires that `alignment_start` < `alingment_end` 
       threshold: int, 
           cut-off value.
       transBlocksNeighbours: list,
           list containing indices of inPlaceBlocks neighbours of each alignment. Output
           of getTransSynOrientation.
    
       Returns
       --------
       outOrderedBlocks: pandas DataFrame,
           Dataframe of type Object. Lower half is NA, upper half contains whether two
           alignments can be connected (True) or not (False).
    """
    if not isinstance(ctx, bool):
        print("CTX status must be a boolean")
        sys.exit()
    if not ctx:
        if transBlocksNeighbours == None:
            sys.exit("ERROR: MISSING transBlocksNeighbours")
        orderedBlocksLen = len(blocksData)
        outOrderedBlocks = np.zeros((orderedBlocksLen,orderedBlocksLen), dtype = object)
        outOrderedBlocks[np.tril_indices(outOrderedBlocks.shape[0],-1)] = np.nan
        outOrderedBlocks[np.triu_indices(outOrderedBlocks.shape[0],0)] = np.False_
        for i in range(orderedBlocksLen):
            for j in range(i+1, orderedBlocksLen):
                if len(np.intersect1d(range(transBlocksNeighbours[i][0]+1,transBlocksNeighbours[i][1]),\
                    range(transBlocksNeighbours[j][0]+1,transBlocksNeighbours[j][1]))) == 0:
                    continue
                outOrderedBlocks[i][j] = SynAndOverlap(blocksData.iloc[j],blocksData.iloc[i], threshold)
        outOrderedBlocks = pd.DataFrame(outOrderedBlocks)
        return outOrderedBlocks
    
    if ctx:
        indices = blocksData.index.values
        orderedBlocksLen = len(indices)
        outOrderedBlocks = np.zeros((orderedBlocksLen, orderedBlocksLen), dtype = object)
        outOrderedBlocks[np.tril_indices(outOrderedBlocks.shape[0],-1)] = np.nan
        outOrderedBlocks[np.triu_indices(outOrderedBlocks.shape[0],0)] = np.False_
        for i in range(len(indices)):
            index = indices[i]
            iData = blocksData.loc[index]
            for j in range(i+1, len(indices)):
                index_2 = indices[j]
                if blocksData.at[index_2,"aChr"] != iData.aChr or blocksData.at[index_2, "bChr"] != iData.bChr:
                    continue
                if blocksData.at[index,"bDir"] != blocksData.at[index_2,"bDir"]:
                    sys.exit("ERROR: bDir not matching")
                if blocksData.at[index, "bDir"] == 1:   ##When input contains ordered blocks
                    outOrderedBlocks[i][j] = SynAndOverlap(blocksData.loc[index_2], blocksData.loc[index], threshold)
                elif blocksData.at[index, "bDir"] == -1:    ##When input contains inverted blocks
                    iData = blocksData.loc[index]
                    jData = blocksData.loc[index_2]
                    if (jData.aStart - iData.aStart) > threshold and (jData.aEnd - iData.aEnd) > threshold and (iData.bStart - jData.bStart) > threshold and (iData.bEnd - jData.bEnd) > threshold:
                        outOrderedBlocks[i][j] = True
                else:
                    sys.exit("ERROR: ILLEGAL BDIR VALUE")
        outOrderedBlocks = pd.DataFrame(outOrderedBlocks)
        outOrderedBlocks.index = blocksData.index
        outOrderedBlocks.columns = blocksData.index.values
        return outOrderedBlocks
                
                
   #%%     
			

    
    



            
def getTransCluster(transGroupIndices, transGenomeAGroups, transGenomeBGroups):
	nodeStack = []
	visitedTransBlock = []
	transCluster = []
	
	for key,value in transGroupIndices.items():
		if key not in visitedTransBlock:
			newGroup = [key]
			visitedTransBlock.append(key)
			node1 = value[0]
			node2 = value[1]
			nodeStack.extend(transGenomeAGroups[node1].member)
			nodeStack.extend(transGenomeBGroups[node2].member)
			
			while len(nodeStack) != 0:
				newKey = nodeStack.pop()
				if newKey not in visitedTransBlock:
					visitedTransBlock.append(newKey)
					newGroup.append(newKey)
					nodeStack.extend(transGenomeAGroups[transGroupIndices[newKey][0]].member)
					nodeStack.extend(transGenomeBGroups[transGroupIndices[newKey][1]].member)
			newGroup.sort()
			transCluster.append(newGroup)
	return(transCluster)


def getBestClusterSubset(cluster, transBlocksData):
    seedBlocks = [i for i in cluster if transBlocksData[i].status == 1]
    if len(cluster) < 50:
        output = bruteSubsetSelector(cluster, transBlocksData, seedBlocks)
        if output == "Failed":
            output = greedySubsetSelector(cluster, transBlocksData, seedBlocks)
    else:
        output = greedySubsetSelector(cluster, transBlocksData, seedBlocks)            
    return output

def bruteSubsetSelector(cluster, transBlocksData, seedBlocks):
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
        if (timeTaken > 10 and timeTaken*(1.5**remainingIterations) > 600):
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
    bestScore = 0
    bestComb = []
    for i in range(iterCount):
        tempCluster = cluster.copy()
        length = len(tempCluster)
        outBlocks =  seedBlocks.copy()
        [tempCluster.remove(i) for i in outBlocks]
        skipList = []
        transBlocksScore = {}
        for i in tempCluster:
            transBlocksScore[i] = (transBlocksData[i].aEnd - transBlocksData[i].aStart) + (transBlocksData[i].bEnd - transBlocksData[i].bStart)
        while len(tempCluster) > 0:
            while len(tempCluster) != length:
                length = len(tempCluster)
                for i in tempCluster:
                    if hasattr(transBlocksData[i],"meTo"):
                        if any(j in outBlocks for j in transBlocksData[i].meTo):
                            tempCluster.remove(i)
                            skipList.append(i)
                    elif hasattr(transBlocksData[i], "meAlist"):
                        if any(j in outBlocks for j in transBlocksData[i].meAlist) and any(j in outBlocks for j in transBlocksData[i].meBlist):
                            tempCluster.remove(i)
                            skipList.append(i)
                            
                
                for i in tempCluster:
                    if hasattr(transBlocksData[i],"meTo"):
                        if all(j in skipList for j in transBlocksData[i].meTo):
                            tempCluster.remove(i)
                            outBlocks.append(i)
                    elif hasattr(transBlocksData[i], "meAlist"):
                        if all(j in skipList for j in transBlocksData[i].meAlist) and all(j in outBlocks for j in transBlocksData[i].meBlist):
                            tempCluster.remove(i)
                            outBlocks.append(i)
 
            if len(tempCluster) > 0:
                topBlocks = sorted(tempCluster, key = lambda x: transBlocksScore[x], reverse = True)[:20]
                totalScore = sum(transBlocksScore[i] for i in topBlocks)
                prob = [transBlocksScore[i]/totalScore for i in topBlocks]
                newBlock = int(np.random.choice(topBlocks, size = 1, p = prob))
                outBlocks.append(newBlock)
                tempCluster.remove(newBlock)
                if hasattr(transBlocksData[newBlock],"meTo"):
                    for i in transBlocksData[newBlock].meTo:
                        if i in tempCluster:
                            tempCluster.remove(i)
                        skipList.append(i)
                elif hasattr(transBlocksData[newBlock],"meAlist"):
                    if any(j in outBlocks for j in transBlocksData[newBlock].meAlist):
                        for k in transBlocksData[newBlock].meBlist:
                            if k in tempCluster:
                                tempCluster.remove(k)
                        skipList.extend(transBlocksData[newBlock].meBlist)
                    elif any(j in outBlocks for j in transBlocksData[newBlock].meBlist):
                        for k in transBlocksData[newBlock].meAlist:
                            if k in tempCluster:
                                tempCluster.remove(k)
                        skipList.extend(transBlocksData[newBlock].meAlist)
                    for meElement in transBlocksData[newBlock].meAlist:
                        if meElement in transBlocksData[newBlock].meBlist:
                            if meElement in tempCluster:
                                tempCluster.remove(meElement)
                            skipList.append(meElement)
#        print(outBlocks)
        bestScore, bestComb = updateBestComb(bestScore, bestComb, outBlocks, transBlocksData)
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
    subs = a[1:,0] - a[:-1,1]    
    overf = (a[:-1,1] - a[1:,1])
    return (a[:,1] - a[:,0]).sum() + subs[subs < 0].sum() + overf[overf > 0].sum()

                
def getscafDict(scaffolds, data, scafSize):
    scafChrDict = {}
    scafCountDict = {}
    scafSizeDict = {}
    percentCutoff = 0.1
    for i in scaffolds:
        scafCountDict[i] = Counter(data.iloc[np.where(data[10] == i)[0], 9])
        #scafData = data.iloc[np.where(data[10] == i)[0],]
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
#            print(alignSizes.values())
            if kruskal(*alignSizes.values())[1] > 0.05:
                scafChrDict[i] = max(scafSizeDict[i].items(), key=lambda x: x[1])[0]
            else:
                size = scafSize[i]
                percentCutoffSize = percentCutoff*size
                bpSizeValues = sorted(bpSizeDict.values())
                identified = 0
                for j in range(len(bpSizeValues)-1,0,-1):
                    if (bpSizeValues[j] - bpSizeValues[j-1]) > percentCutoffSize:
                        scafChrDict[i] = list(bpSizeDict.keys())[list(bpSizeDict.values()).index(bpSizeValues[j])]
                        identified = 1
                        break
                if identified == 0:
                    meanAlignSize = {}
                    for j in uniChr:
                        meanAlignSize[j] = scafSizeDict[i]/scafCountDict[i][j]
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
    
def getTransClasses(clusterSolutionBlocks, transData):
    def setTL(j):
        if transData[j].dir == 1:                   
            transClasses["translocation"].append(j)
        elif transData[j].dir == -1:
            transClasses["invTranslocation"].append(j)
        else:
            print("ERROR ERROR ERROR", j)
            
    def setDup(j):
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
                    setDup(j)
                elif transData[j].aUni and transData[j].bUni:
                    if transData[j].genomeAUni and transData[j].genomeBUni:
                        setTL(j)
                    elif not transData[j].genomeAUni:
                        isTrans = 1
                        for k in transData[j].genomeAMembers:
                            if k in i:
                                if getScore([k], transData) > getScore([j],transData):
                                    isTrans = 0
                                    break
                        if isTrans:
                            setTL(j)
                        else:
                            setDup(j)
                    elif not transData[j].genomeBUni:
                        isTrans = 1
                        for k in transData[j].genomeBMembers:
                            if k in i:
                                if getScore([k], transData) > getScore([j],transData):
                                    isTrans = 0
                                    break
                        if isTrans:
                            setTL(j)
                        else:
                            setDup(j)
            elif not transData[j].aUni or not transData[j].bUni:
                setDup(j)
            elif transData[j].aUni and transData[j].bUni:
                if hasattr(transData[j],"meTo"):
                    if len(np.intersect1d(transData[j].meTo, i)) > 0:
                        setDup(j)
                    else:
                        setTL(j)
                elif hasattr(transData[j],"meAlist"):
                    if len(np.intersect1d(transData[j].meAlist, i)) > 0 or\
                    len(np.intersect1d(transData[j].meBlist, i)) > 0:
                        setDup(j)
                    else:
                        setTL(j)
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

def outSyn(cwdPath, threshold):
    reCoords = pd.DataFrame(columns=["aStart","aEnd","bStart","bEnd","aChr","bChr"])
    ctxAnnoDict = {"duplication":"dupCtx",
                   "invDuplication":"invDupCtx",
                   "translocation":"TLCtx",
                   "invTranslocation":"invTLCtx"}
    reCoords =  pd.DataFrame()
        
    synData = []
    with open(cwdPath+"synOut.txt","r") as fin:
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
    if max([len(i) for i in synData]) == 7:
        synData = pd.DataFrame(synData,columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr","isinInv"])
    else:
        synData = pd.DataFrame(synData,columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"])
        synData["isinInv"] = ""
        
    synData["class"] = "syn"
       
    for i in ["invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt","ctxOut.txt"]:    
        data = []
        with open(cwdPath+i,"r") as fin: 
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
    
    with open(cwdPath+"synOut.txt","w") as fout:
        for i in outClusters:
            fout.write("\t".join(map(str,["#",allBlocks.at[i[0],"aChr"],allBlocks.at[i[0],"aStart"],allBlocks.at[i[-1],"aEnd"],"-",allBlocks.at[i[0],"aChr"],allBlocks.at[i[0],"bStart"],allBlocks.at[i[-1],"bEnd"],"\n"])))
            for j in i:
                fout.write("\t".join(map(str,allBlocks.loc[j][0:4])))
                if synData.loc[synLocs[j]]["isinInv"] == "Syn_in_Inv":
                    fout.write("\tSyn_in_Inv\n")
                else:
                    fout.write("\n")   
    return None
    
        
def groupSyn(tempInvBlocks, dupData, invDupData, invTLData, TLData, threshold, synData):
    
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
    
    for i in [tempInvBlocks, tempInvDupData, tempInvTLData, tempTLData, tempDupData]:
        if len(i) > 0:
            allBlocks = allBlocks.append(i)
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


def mergeOutputFiles(uniChromo,path):
    def addData(fName,anno, chromo):
        fPath = open(path+chromo+"_"+anno+"Out.txt","r")
        for line in fPath.readlines():
            line = line.strip().split("\t")
            if line[0] == "#":
                fName.write("\t".join(unlist([line[0], chromo, line[1:4], chromo, line[4:]])) + "\n")
            else:
                fName.write("\t".join(line) + "\n")
        fPath.close()
        fileRemove(path+chromo+"_"+anno+"Out.txt")
                
    fSyn = open(path+"synOut.txt","w")
    fInv = open(path+"invOut.txt","w")
    fTL = open(path+"TLOut.txt","w")
    fInvTL = open(path+"invTLOut.txt","w")
    fDup = open(path+"dupOut.txt","w")
    fInvDup = open(path+"invDupOut.txt","w")
    
    files = [fSyn, fInv, fTL, fInvTL, fDup, fInvDup]
    classes = ["syn","inv","TL","invTL","dup","invDup"]
    
    for chromo in uniChromo:
        for i in range(len(classes)):
            addData(files[i], classes[i], chromo)
            
    for f in files:
        f.close()

def readAnnoCoords(cwdPath, uniChromo):
    annoCoords = pd.DataFrame(columns=["aStart","aEnd","bStart","bEnd","aChr","bChr"])
    synData = []
    fin = open(cwdPath+"synOut.txt","r")
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
        fin = open(cwdPath+i,"r")
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
    def __init__(self,aStart, aEnd, bStart, bEnd, Dir, transClusterIndex, i):
        self.aStart = aStart
        self.aEnd = aEnd
        self.bStart = bStart
        self.bEnd = bEnd
        self.dir = Dir
#        self.orderedBlocksIndex = orderedBlocksIndex
        self.transBlocksID = i
        self.transClusterIndex = transClusterIndex
        self.status = 0
        self.overlappingInPlaceBlocks = []
    
    def addGenomeGroupMembers(self,transGenomeAGroups, transGenomeBGroups):
        self.genomeAMembers = list(set(transGenomeAGroups[self.transGroupIndices[0]].member)\
                                       - set([self.transBlocksID]))
        self.genomeBMembers = list(set(transGenomeBGroups[self.transGroupIndices[1]].member)\
                                       - set([self.transBlocksID]))
        self.genomeAUni = True if len(self.genomeAMembers) == 0 else False
        self.genomeBUni = True if len(self.genomeBMembers) == 0 else False
        
    def addOrderedData(self, orderedData):
        self.orderedData = orderedData
        
    def addTransGroupIndices(self, indices):
        self.transGroupIndices = indices
        
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
        aBlocks = list(np.intersect1d(np.where(inPlaceBlocks.aStart.values < self.aEnd)[0],\
                                      np.where(inPlaceBlocks.aEnd.values > self.aStart)[0]))

        blockAUni = 0
        
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
        
        self.overlappingInPlaceBlocks.append(aBlocks)
        self.aUni = True if blockAUni > 1000 or blockAUni > 0.5*(self.aEnd-self.aStart) else False       
        
    def checkOverlapWithSynBlocks_B(self,inPlaceBlocks, threshold):
        bBlocks = list(np.intersect1d(np.where(inPlaceBlocks.bStart.values < self.bEnd)[0],\
                                      np.where(inPlaceBlocks.bEnd.values > self.bStart)[0]))    
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
        
        try:
            self.meTo.extend(blockID) if type(blockID) == list else self.meTo.append(blockID)
        except AttributeError:
            self.meTo = blockID if type(blockID) == list else [blockID]
            
    def setMEList(self, meAlist, meBlist):
        """Lists of a-overlap and b-overlap blocks. If at least 1 block has
        been selected from both lists then this block would become redundant"""
        self.meAlist = meAlist
        self.meBlist = meBlist
        
    def setStatus(self,stat):
        """stat = 1 ==> transBlock is important/necessary/unique"""
        self.status = stat
#%%
#################################################################
### SV identification functions
#################################################################
        
def readSVData(cwdPath):
    annoCoords = pd.DataFrame()
    for fileType in ["syn","inv","TL","invTL"]:
        try:
            fileData = pd.read_table(cwdPath+fileType+"Out.txt", header = None, dtype = object)
        except pd.io.common.EmptyDataError:
            print(fileType+"Out.txt is empty. Skipping analysing it.")
            continue
        except Exception as e:
            print("ERROR: while trying to read ", fileType, "Out.txt", e)
            continue
            
        annoIndices = np.where(fileData[0] =="#")[0]
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
        coordsData["aChr"] = list(np.repeat(annoData[1],repCount))
        coordsData["bChr"] = list(np.repeat(annoData[5],repCount))
        coordsData["state"] = fileType
        annoCoords = annoCoords.append(coordsData.copy())
                                
    try:
        fileData = pd.read_table(cwdPath+"ctxOut.txt", header = None, dtype = object)
        annoIndices = np.where(fileData[0] =="#")[0]
        states = np.array(fileData[8].loc[annoIndices], dtype = "str")
        aChr = np.array(fileData[1].loc[annoIndices], dtype = "str")
        bChr = np.array(fileData[5].loc[annoIndices], dtype = "str")
        coordsData = fileData.loc[fileData[0] !="#"].copy()
        annoIndices = np.append(annoIndices,len(fileData))
        repCount = annoIndices[1:] - annoIndices[:-1] - 1
            
        reps = np.repeat(range(len(annoIndices)-1), repCount)
        stateReps = np.repeat(states, repCount)
        aChrReps = np.repeat(aChr, repCount)
        bChrReps = np.repeat(bChr, repCount)
        
        coordsData1 = coordsData[[0,1,2,3]].astype(dtype = "int64")
        coordsData1["aChr"] = aChrReps
        coordsData1["bChr"] = bChrReps
        coordsData1["group"] = reps
        coordsData1["state"] = stateReps
        coordsData1 = coordsData1[[0,1,2,3,"group","aChr","bChr","state"]]
        coordsData1 = coordsData1.loc[coordsData1["state"].isin(["translocation","invTranslocation"])]
        coordsData1.loc[coordsData1.state == "translocation","state"] = "ctx"
        coordsData1.loc[coordsData1.state == "invTranslocation","state"] = "invCtx"
        annoCoords = annoCoords.append(coordsData1)
    except pd.io.common.EmptyDataError:
        print("ctxOut.txt is empty. Skipping analysing it.")
    except Exception as e:
        print("ERROR: while trying to read ctxOut.txt", e)
        
    annoCoords.columns = ["aStart","aEnd","bStart","bEnd","group","aChr","bChr","state"]
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))
    return annoCoords



def getSV(cwdPath, allAlignments):
    fout = open(cwdPath+"sv.txt","w")
    allAlignments["id"] = allAlignments.group.astype("str") + allAlignments.aChr + allAlignments.bChr + allAlignments.state
    allBlocks = pd.unique(allAlignments.id)

    for i in allBlocks:
        blocksAlign = allAlignments.loc[allAlignments.id == i].copy()
        ordered = 1 if "inv" not in blocksAlign.state.iloc[0] else 0
        for j in range(len(blocksAlign) - 1):
            m = blocksAlign.iat[j+1,0] - blocksAlign.iat[j,1] - 1
            if ordered:
                n = blocksAlign.iat[j+1,2] - blocksAlign.iat[j,3] - 1
            else:
                n = blocksAlign.iat[j,3] - blocksAlign.iat[j+1,2] - 1
                
            if m == 0:
                
                if n == 0:
                    continue
                
                elif n > 0:
                    if ordered:
                        fout.write("\t".join(["InDel",
                                              str(blocksAlign.iat[j,1]),
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,3] + 1),
                                              str(blocksAlign.iat[j+1, 2] - 1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                    else:
                        fout.write("\t".join(["InDel",
                                              str(blocksAlign.iat[j,1]),
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,3] - 1),
                                              str(blocksAlign.iat[j+1, 2] + 1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                
                elif n < 0:
                    if ordered:
                        j_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
                        j1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])
                        sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                        fout.write("\t".join(["CNV",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                    else:
                        j_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
                        j1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
                        sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                        fout.write("\t".join(["CNV",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                    
            elif m == 1:
                
                if n == 0:
                    if ordered:
                        fout.write("\t".join(["InDel",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]),
                                              str(blocksAlign.iat[j+1,2]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                    else:
                        fout.write("\t".join(["InDel",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]),
                                              str(blocksAlign.iat[j+1, 2]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                        
                elif n==1:
                    if ordered:
                        fout.write("\t".join(["SNP",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]+1),
                                              str(blocksAlign.iat[j+1,2]-1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                    else:
                        fout.write("\t".join(["SNP",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]-1),
                                              str(blocksAlign.iat[j+1,2]+1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                elif n>1:
                    if ordered:
                        fout.write("\t".join(["InDel+SNP",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]+1),
                                              str(blocksAlign.iat[j+1,2]-1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "SNP:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
                    else:
                        fout.write("\t".join(["InDel+SNP",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]-1),
                                              str(blocksAlign.iat[j+1,2]+1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "SNP:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
                elif n<0:
                    if ordered:
                        j_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
                        j1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])
                        sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                        fout.write("\t".join(["CNV+InDel",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
                    else:
                        j_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
                        j1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
                        sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                        fout.write("\t".join(["CNV+InDel",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
        
            elif m>1:
                
                if n==0:
                    if ordered:
                        fout.write("\t".join(["InDel",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]),
                                              str(blocksAlign.iat[j+1,2]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                    else:
                        fout.write("\t".join(["InDel",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]),
                                              str(blocksAlign.iat[j+1,2]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                elif n==1:
                    if ordered:
                        fout.write("\t".join(["InDel+SNP",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]+1),
                                              str(blocksAlign.iat[j+1,2]-1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "SNP:Q:"+str(blocksAlign.iat[j,3]+1)+"-"+str(blocksAlign.iat[j+1,2]-1)]) + "\n")
                    else:
                        fout.write("\t".join(["InDel+SNP",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]-1),
                                              str(blocksAlign.iat[j+1,2]+1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "SNP:Q:"+str(blocksAlign.iat[j,3]-1)+"-"+str(blocksAlign.iat[j+1,2]+1)]) + "\n")
                elif n>1:
                    if ordered:
                        fout.write("\t".join(["HDR",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]+1),
                                              str(blocksAlign.iat[j+1,2]-1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) +"\n")
                    else:
                        fout.write("\t".join(["HDR",
                                              str(blocksAlign.iat[j,1]+1),
                                              str(blocksAlign.iat[j+1,0]-1),
                                              str(blocksAlign.iat[j,3]-1),
                                              str(blocksAlign.iat[j+1,2]+1),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) +"\n")
                elif n<0:
                    if ordered:
                        j_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
                        j1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])
                        sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                        fout.write("\t".join(["CNV+InDel",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
                    else:
                        j_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
                        j1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
                        sCoord = round(blocksAlign.iat[j,1] - j_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,0] + j1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                        fout.write("\t".join(["CNV+InDel",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "InDel:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(blocksAlign.iat[j+1,0]-1)]) + "\n")
            
            elif m<0:
                
                j_prop = abs(m) / (blocksAlign.iat[j,1] - blocksAlign.iat[j,0])
                j1_prop = abs(m) / (blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])
                
                if n==0:
                    if ordered:
                        sCoord = round(blocksAlign.iat[j,3] - j_prop*(blocksAlign.iat[j,3] - blocksAlign.iat[j,2])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,2] + j1_prop*(blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])).astype(int)
                        fout.write("\t".join(["CNV",
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                    else:
                        sCoord = round(blocksAlign.iat[j,3] + j_prop*(blocksAlign.iat[j,2] - blocksAlign.iat[j,3])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,2] - j1_prop*(blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])).astype(int)
                        fout.write("\t".join(["CNV",
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6]]) + "\n")
                
                if n>0:
                    if ordered:
                        sCoord = round(blocksAlign.iat[j,3] - j_prop*(blocksAlign.iat[j,3] - blocksAlign.iat[j,2])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,2] + j1_prop*(blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])).astype(int)
                        fout.write("\t".join(["CNV+InDel",
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "InDel:Q:"+str(blocksAlign.iat[j,3]+1)+"-"+str(blocksAlign.iat[j+1,2]-1)]) + "\n")
                    else:
                        sCoord = round(blocksAlign.iat[j,3] + j_prop*(blocksAlign.iat[j,2] - blocksAlign.iat[j,3])).astype(int)
                        eCoord = round(blocksAlign.iat[j+1,2] - j1_prop*(blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])).astype(int)
                        fout.write("\t".join(["CNV+InDel",
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "InDel:Q:"+str(blocksAlign.iat[j,3]-1)+"-"+str(blocksAlign.iat[j+1,2]+1)]) + "\n")
    
                if n<0:
                    maxOverlap = max(abs(m),abs(n))
                    if abs(m-n) < 0.1*maxOverlap: ## no SV if the overlap on both genomes is of similar size
                        continue
                    
                    if abs(m) > abs(n):
                        if ordered:                   
                            sCoord = round(blocksAlign.iat[j,3] - j_prop*(blocksAlign.iat[j,3] - blocksAlign.iat[j,2])).astype(int)
                            eCoord = round(blocksAlign.iat[j+1,2] + j1_prop*(blocksAlign.iat[j+1,3] - blocksAlign.iat[j+1,2])).astype(int)
                            fout.write("\t".join(["CNV+Tandem",
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "Tandem:Q:"+str(blocksAlign.iat[j,3]+1)+"-"+str(eCoord)]) + "\n")
                        
                        else:
                            sCoord = round(blocksAlign.iat[j,3] + j_prop*(blocksAlign.iat[j,2] -blocksAlign.iat[j,3])).astype(int)
                            eCoord = round(blocksAlign.iat[j+1,2] - j1_prop*(blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,2])).astype(int)
                            fout.write("\t".join(["CNV+Tandem",
                                              str(blocksAlign.iat[j+1,0]),
                                              str(blocksAlign.iat[j,1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "Tandem:Q:"+str(blocksAlign.iat[j,3]-1)+"-"+str(eCoord)]) + "\n")
                    else:
                        if ordered:
                            k_prop = abs(n) / (blocksAlign.iat[j,3] - blocksAlign.iat[j,2])
                            k1_prop = abs(n) / (blocksAlign.iat[j+1,3] - blocksAlign.iat[j,2])
                            sCoord = round(blocksAlign.iat[j,1] - k_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                            eCoord = round(blocksAlign.iat[j+1,0] + k1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                            fout.write("\t".join(["CNV+Tandem",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "Tandem:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(eCoord)]) + "\n")
                        else:
                            k_prop = abs(n) / (blocksAlign.iat[j,2] - blocksAlign.iat[j,3])
                            k1_prop = abs(n) / (blocksAlign.iat[j+1,2] - blocksAlign.iat[j+1,3])
                            sCoord = round(blocksAlign.iat[j,1] - k_prop*(blocksAlign.iat[j,1] - blocksAlign.iat[j,0])).astype(int)
                            eCoord = round(blocksAlign.iat[j+1,0] + k1_prop*(blocksAlign.iat[j+1,1] - blocksAlign.iat[j+1,0])).astype(int)
                            fout.write("\t".join(["CNV+Tandem",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j+1,2]),
                                              str(blocksAlign.iat[j,3]),
                                              blocksAlign.iat[0,5],
                                              blocksAlign.iat[0,6],
                                              "Tandem:R:"+str(blocksAlign.iat[j,1]+1)+"-"+str(eCoord)]) + "\n")
                            
    fout.close()
    return None                



def getNotAligned(cwdPath):    
    annoCoords = pd.DataFrame()
    for fileType in ["syn","inv", "TL", "invTL","dup", "invDup"]:
        try:
            fileData = pd.read_table(cwdPath+fileType+"Out.txt", header = None, dtype = object)  
        except pd.io.common.EmptyDataError:
            print(fileType, "Out.txt is empty. Skipping analysing it.")
            continue
        except Exception as e:
            print("ERROR: while trying to read ", fileType, "Out.txt", e)
            continue
        
        
        
#        fileData = pd.read_table(cwdPath+fileType+"Out.txt", header = None, dtype = object)
        coordsData = fileData.loc[fileData[0] == "#",[2,3,6,7,1,5]].copy()
        coordsData[[2,3,6,7]] = coordsData[[2,3,6,7]].astype(dtype="int64")
        coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
        annoCoords = annoCoords.append(coordsData.copy())

    fileData = pd.read_table(cwdPath+"ctxOut.txt", header = None, names = list(range(11)), dtype = object, sep ="\t")
    coordsData = fileData.loc[fileData[0] == "#"]
    coordsData = coordsData[[2,3,6,7,1,5]].copy()
    coordsData[[2,3,6,7]] = coordsData[[2,3,6,7]].astype(dtype="int64")
    coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
    
    annoCoords = annoCoords.append(coordsData.copy())
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))
  
    fout = open(cwdPath + "notAligned.txt","w")
    
#    fout.write("#Reference Genome\n")
    df = annoCoords[["aStart","aEnd","aChr"]].copy()
    df.sort_values(["aChr", "aStart", "aEnd"], inplace = True)
    for chrom in sorted(annoCoords.aChr.unique()):
        chromData = df.loc[df.aChr == chrom]
        maxEnd = chromData.iloc[0,1]
        for row in chromData.itertuples(index = False):
            if row.aStart > maxEnd+1:
                fout.write("\t".join(["R",str(maxEnd+1),
                                      str(row.aStart - 1),
                                      chrom]) + "\n")
            if row.aEnd > maxEnd:
                maxEnd = row.aEnd
    
#    fout.write("\n")
    
#    fout.write("#Query Genome\n")
    df = annoCoords[["bStart","bEnd","bChr"]].copy()
    df.sort_values(["bChr", "bStart", "bEnd"], inplace = True)
    for chrom in sorted(annoCoords.bChr.unique()):
        chromData = df.loc[df.bChr == chrom]
        maxEnd = chromData.iloc[0,1]
        for row in chromData.itertuples(index = False):
            if row.bStart > maxEnd+1:
                fout.write("\t".join(["Q",str(maxEnd+1),
                                      str(row.bStart - 1),
                                      chrom]) + "\n")
            if row.bEnd > maxEnd:
                maxEnd = row.bEnd
                
    fout.close()
    return None


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
#        print(fileData.head())
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
#        print(row.end, terminalData[row.bGenome].end)
        if row.end <= terminalData[row.bGenome].end:
            print("values must be sorted. Invalid Entry: ",row, terminalData[row.bGenome])
            sys.exit()
        if row.start <= terminalData[row.bGenome].end:
            terminalData[row.bGenome].start = terminalData[row.bGenome].end + 1
            terminalData[row.bGenome].end = row.end
        else:
            terminalData[row.bGenome].start = row.start
            terminalData[row.bGenome].end = row.end
#        print(terminalData)
        if max(terminalData.loc["start"]) < min(terminalData.loc["end"]):
            regions.append((max(terminalData.loc["start"]),min(terminalData.loc["end"])))
#        if count == 20:
##            break
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
#            nodeCount = [sum(self.Z[t,:partSize[0]])]
            nodeCount = [sum(self.getBits(t,i)) for i in self.BsOut]
            return self.BsOut[nodeCount.index(min(nodeCount))]
        
        def bitCLQ(self, t):
            bt = self.getBt(self.partSize, t)
#            print(bt)
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
#                    print("found clique:", self.S)
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
    

#%%
def plotBlocks(blocksData):
	blocksData = [orderedBlocks.iloc[[250,251,255]], orderedBlocks.iloc[[1370, 1371]]]
	
	blocksDataOri = [i.copy() for i in blocksData]

	
	blocksCoords = {}
	for i in range(len(blocksData)):		
		bData = blocksData[i].copy()
		aMin = min(bData[['aStart','aEnd']].min())
		aMax = max(bData[['aStart','aEnd']].max())
		bMin = min(bData[['bStart','bEnd']].min())
		bMax = max(bData[['bStart','bEnd']].max())
		blocksCoords[i] = [aMin,aMax, bMin, bMax]		
	
	keyOrder = sorted(blocksCoords.keys(), key = lambda x : blocksCoords[x][0])
	
	gapsRemoved = []
#	startPosition = blocksCoords[keyOrder[0]][0]
	maxEnd = blocksCoords[keyOrder[0]][1]
	for i in range(1,len(keyOrder)):
		if blocksCoords[keyOrder[i]][0] > maxEnd:
			gapsRemoved.append(blocksCoords[keyOrder[i]][0] - blocksCoords[keyOrder[i-1]][1])
		else:
			gapsRemoved.append(0)
			
	for i in range(1, len(keyOrder)):
		leftShift = sum(gapsRemoved[:i])
		rightShift = max(0,log(leftShift,1.015))
		blocksData[keyOrder[i]]['aStart'] = blocksData[keyOrder[i]]['aStart'] -  leftShift + rightShift
		blocksData[keyOrder[i]]['aEnd'] = blocksData[keyOrder[i]]['aEnd'] -  leftShift + rightShift
	
	dataLimits = [[min(i[['aStart','aEnd']].min()), max(i[['aStart','aEnd']].max()),min(i[['bStart','bEnd']].min()), max(i[['bStart','bEnd']].max())] for i in blocksData]
	aMin = min(unlist([x[0:2] for x in dataLimits]))
	aMax = max(unlist([x[0:2] for x in dataLimits]))
	bMin = min(unlist([x[2:4] for x in dataLimits]))
	bMax = max(unlist([x[2:4] for x in dataLimits]))

	
#	
#	for i in range(1,len(blocksCoords)):
#		if blocksCoords[i][0] > blocksCoords[i-1][1] 
#	
	colList = ['aStart', 'aEnd', 'bStart', 'bEnd']	
	for bData in blocksData:
		for i in colList[:2]:
			bData[i] = (bData[i] - aMin)/(aMax - aMin)
			
		for i in colList[2:]:
			bData[i] = (bData[i] - bMin)/(bMax - bMin)
	
	colors = getColors(plt.cm.Dark2,len(blocksData))
	
	bLen = len(blocksData)
	
	for count in range(bLen):
		bData = blocksData[count]
		for i in range(bData.shape[0]):
			row = bData.iloc[i]
			plt.plot([row[0], row[1]], [0 - (0.1*count),0 - (0.1*count)], linewidth = 5, color = colors[count], path_effects = [pe.Stroke(linewidth = 10, foreground="k"),pe.Normal()])
			plt.plot([row[2], row[3]], [1 + 0.1*(bLen-1-count),1 + 0.1*(bLen-1-count)], linewidth = 5, color = colors[count],path_effects = [pe.Stroke(linewidth = 10, foreground="k"),pe.Normal()])
	plt.show()
	plt.gca().add_patch(patches.Rectangle((10,10), 30,5))
	plt.show()
	plt.plot()

	
