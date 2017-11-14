#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 12:00:19 2017

@author: goel
"""

import os
import sys
import pandas as pd
import numpy as np
from myUsefulFunctions import *
from synsearchFunctions import *
from itertools import cycle, product
from collections import Counter
from igraph import *
import datetime

def filterBreaks(breaks):
    while True:
        starts = breaks[:-1]
        ends = breaks[1:]
        lengths = np.array(subList(ends, starts))
        small = np.where(lengths < 100)[0]
        print(len(small))
        if len(small) == 0:
            break
        else:
            for i in small[::-1]:
                breaks[i] = int(breaks[i] + breaks[i+1])/2
                breaks = np.delete(breaks,i+1)
    return list(breaks)
    
def getMatchRegion(gen, regID):
    nodeID = adders[gen]+regID
    matches = np.where(regGMatrix.iloc[nodeID] == 1)[0]
    
    for genome in genNames:
        if genome != gen:
            matID = [i - adders[genome] for i in matches]
            for i in matID:
                if i >= 0 and i <= len(genomeRegions[genome]):
                    print(genome,genomeRegions[genome].iloc[i].tolist(), i)
    
    
#%%
dirNames = "an1_chr_c24_chr,col_an1_Chr,col_c24_Chr,col_cvi_Chr,col_eri_Chr,col_kyo_Chr,col_ler_Chr,cvi_chr_an1_chr,cvi_chr_c24_chr,cvi_chr_eri_chr,cvi_chr_kyo_chr,eri_chr_an1_chr,eri_chr_c24_chr,kyo_chr_an1_chr,kyo_chr_c24_chr,kyo_chr_eri_chr,ler_chr_an1_chr,ler_chr_c24_chr,ler_chr_cvi_chr,ler_chr_eri_chr,ler_chr_kyo_chr"

dirNames = dirNames.split(",")

rootDir = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/"

pathList = [rootDir + i for i in dirNames ]

genNames = ["col","ler","cvi","kyo","eri","an1","c24"]

genData = {i: [[ j for j in dirNames if i in j]] for i in genNames }

for key, values in genData.items():
#    print(values)
    genData[key].append([ j for i in values[0] for j in genNames if j != key and j in i])
    genData[key].append(["a" if j.index(key) == 0 else "b" for j in values[0]])

fileTypes = ("syn","inv","TL","invTL","dup","invDup")
fileNames = ["".join(i) for i in zip(cycle(["Chr1_"]),fileTypes,cycle(["Out.txt"]))]

dirFiles = {i:[rootDir + i + "/" + j for j in fileNames] for i in dirNames}

dirData = {i: getBlocksData(dirFiles[i], list(fileTypes)) for i in dirNames}

#lerDataSyn = lerData.loc[lerData.state == "syn"]
#cviDataSyn = cviData.loc[cviData.state == "syn"]
#kyoDataSyn = kyoData.loc[kyoData.state == "syn"]
#lerCviDataSyn = lerCviData.loc[lerCviData.state == "syn"]


###########################################################
## Make regions using ALL blocks
###########################################################


genDataList = {i:getDataList([dirData[j] for j in genData[i][0]], genData[i][1], genData[i][2]) for i in genNames}
genRegions = {i:getConservedRegions(genDataList[i], isSyn = False) for i in genNames}

#
#colDataList = getDataList([lerData, cviData],["ler","cvi"],["a","a"])
#colRegions = getConservedRegions(colDataList, isSyn = False)
#
#lerDataList = getDataList([lerData,lerCviData],["col","cvi"],["b","a"])
#lerRegions = getConservedRegions(lerDataList, isSyn = False)
#
#cviDataList = getDataList([cviData, lerCviData],["col","ler"],["b","b"])
#cviRegions = getConservedRegions(cviDataList, isSyn = False)
#
#genomeData = {"col":colDataList,"ler": lerDataList, "cvi" :cviDataList}
#genomeRegions = {"col":colRegions, "ler":lerRegions, "cvi":cviRegions}

genomeData = genDataList
genomeRegions = genRegions

#regG = Graph().as_undirected()
#regG.add_vertices(sum(map(len,genomeRegions.values())))

adders = [0]
for i in genNames[:-1]:
    adders.append(adders[-1] + len(genomeRegions[i]))
    print(adders)
adders = dict(zip(genNames,adders))

regGMatrix = np.zeros([sum(map(len,genomeRegions.values()))]*2, dtype = "bool")

 
for i in genNames:
    print(datetime.datetime.now())
    for j in genomeData[i].itertuples():
        aStart = j.start
        aEnd = j.end
        bStart = j.bStart
        bEnd = j.bEnd
        bGenome = j.bGenome
        state = j.state
        aRegionsIndices = np.intersect1d(np.where(genomeRegions[i].start >= aStart)[0],np.where(genomeRegions[i].end <= aEnd)[0])
        bRegionsIndices = np.intersect1d( np.where(genomeRegions[bGenome].start >= bStart)[0], np.where(genomeRegions[bGenome].end <= bEnd)[0])
#        aRegions = genomeRegions[i].iloc[aRegionsIndices]
#        bRegions = genomeRegions[bGenome].iloc[bRegionsIndices]        
        for a in aRegionsIndices:
            ast = genomeRegions[i].iat[a,0]
            aen = genomeRegions[i].iat[a,1]
            
            aStartPro = (ast - aStart)/(aEnd - aStart)
            aEndPro = (aen - aStart)/(aEnd - aStart)
            if state in ["syn","TL","dup"]:
                bst_ex = bStart + (aStartPro*(bEnd - bStart))
                ben_ex = bStart + (aEndPro*(bEnd - bStart))
            elif state in ["inv","invTL","invDup"]:
                bst_ex = bStart + ((1-aEndPro)*(bEnd - bStart))
                ben_ex = bStart + ((1-aStartPro)*(bEnd - bStart))
            
            for b in bRegionsIndices:
                bst = genomeRegions[bGenome].iat[b,0]
                ben = genomeRegions[bGenome].iat[b,1]
                
                if bst_ex > ben or bst > ben_ex:
                    continue
                else:
                    st = max(bst,bst_ex)
                    en = min(ben,ben_ex)
                    length = en - st
                    if length >= 0.75*min(ben-bst, ben_ex-bst_ex):
#                        if adders[i]+a == 2 and adders[bGenome]+b == adders["ler"]+2:
#                            print(i,bGenome, j, a, b)
                        regGMatrix[adders[i]+a,adders[bGenome]+b] = True
                        
#regG.simplify()
#regGMatrix = pd.DataFrame(regG.get_adjacency().data)
regGMatrix = pd.DataFrame(regGMatrix)
clqs = getCLQ(regGMatrix, list(map(lambda x: len(genomeRegions[x]), genNames)))


coords = []
for i in clqs:
    clqCoords = []
    for j in range(len(i)):
        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].start)
        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].end)
    coords.append(clqCoords)
        
allData = pd.DataFrame(coords)
clqsData = pd.DataFrame(clqs)

clqsData.columns = genNames
allData.columns = [i + j for i in genNames for j in ["_start","_end"]]

allData.to_csv("allData.txt", sep="\t", header = True, index = False)
clqsData.to_csv("clqsData.txt", sep="\t", header = True, index = False)








#clqs_filtered = []
#for i in clqs:
#    count = 0
#    for j in range(len(i)):
#        if genomeRegions[genNames[j]].iloc[i[j]].end - genomeRegions[genNames[j]].iloc[i[j]].start > 100:
#            count+=1
#    if count == 3:
#        clqs_filtered.append(i)
#    
#coords = []
#for i in clqs_filtered:
#    clqCoords = []
#    for j in range(len(i)):
#        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].start)
#        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].end)
#    coords.append(clqCoords)    
#
#allData_filtered = pd.DataFrame(coords)
#
#sizeData = []
#for i in clqs:
#    _d = []
#    for j in range(len(genNames)):
#        _d.append(genomeRegions[genNames[j]].iloc[i[j]].end - genomeRegions[genNames[j]].iloc[i[j]].start)
#    sizeData.append(_d)
#        
#    
#transData = []
#for i in allDataAll.itertuples(index = False):
#    a = np.intersect1d(np.where(allData[0] <= i[0])[0],np.where(allData[1] >= i[1])[0])
#    b = np.intersect1d(np.where(allData[2] <= i[2])[0],np.where(allData[3] >= i[3])[0])
#    c = np.intersect1d(np.where(allData[4] <= i[4])[0],np.where(allData[5] >= i[5])[0])
#    
##    if not len(intersect(a,b,c)) > 0:
#    if len(a) < 1 or len(b) < 1 or len(c) < 1:
#        transData.append(i)
#
#transData = pd.DataFrame(transData)
#coords = []
#for i in transData.itertuples(index = False):
#    if (i[1] - i[0] > 100) and (i[3] - i[2] > 100) and (i[5] - i[4] > 100):
#      coords.append(i)
#transData_filtered = pd.DataFrame(coords)   
