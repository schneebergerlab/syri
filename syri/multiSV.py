#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 11:12:56 2017

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
lerPath = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/col_ler_Chr/"
cviPath = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/col_cvi_Chr/"
lerCviPath = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/ler_chr_cvi_chr/"
kyoPath = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/col_kyo_Chr/"
lerKyoPath = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/ler_chr_kyo_chr/"
cviKyoPath = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/cvi_chr_kyo_chr/"

fileTypes = ("syn","inv","TL","invTL","dup","invDup")
fileNames = ["".join(i) for i in zip(cycle(["Chr1_"]),fileTypes,cycle(["Out.txt"]))]
        
lerFiles = [lerPath + i for i in fileNames]
cviFiles = [cviPath + i for i in fileNames]
kyoFiles = [kyoPath + i for i in fileNames]
lerCviFiles = [lerCviPath + i for i in fileNames]
lerKyoFiles = [lerKyoPath + i for i in fileNames]
cviKyoFiles = [cviKyoPath + i for i in fileNames]

lerData = getBlocksData(lerFiles, list(fileTypes))
cviData = getBlocksData(cviFiles, list(fileTypes))
kyoData = getBlocksData(kyoFiles, list(fileTypes))
lerCviData = getBlocksData(lerCviFiles, list(fileTypes))
lerKyoData = getBlocksData(lerKyoFiles, list(fileTypes))
cviKyoData = getBlocksData(cviKyoFiles, list(fileTypes))


lerDataSyn = lerData.loc[lerData.state == "syn"]
cviDataSyn = cviData.loc[cviData.state == "syn"]
kyoDataSyn = kyoData.loc[kyoData.state == "syn"]
lerCviDataSyn = lerCviData.loc[lerCviData.state == "syn"]


colBreaks = np.unique(unlist([lerDataSyn.aStart.tolist(), lerDataSyn.aEnd.tolist(), cviDataSyn.aStart.tolist(), cviDataSyn.aEnd.tolist()]))
lerBreaks = np.unique(unlist([lerDataSyn.bStart.tolist(),lerDataSyn.bEnd.tolist(), lerCviDataSyn.aStart.tolist(), lerCviDataSyn.aEnd.tolist()]))
cviBreaks = np.unique(unlist([cviDataSyn.bStart.tolist(), cviDataSyn.bEnd.tolist(),lerCviDataSyn.bStart.tolist(),lerCviDataSyn.bEnd.tolist()]))



colBreaks = filterBreaks(colBreaks)
lerBreaks = filterBreaks(lerBreaks)
cviBreaks = filterBreaks(cviBreaks)

colRegions = genSeqMap(colBreaks, "Chr1", "Col")
lerRegions = genSeqMap(lerBreaks,"Chr1","Ler")
cviRegions = genSeqMap(cviBreaks,"Chr1","Cvi")


colRIn = 0
lerRIn = 0

#dataList = lerDataSyn[["aStart","aEnd"]]
#bGenome = ["ler"]*len(dataList)
#
#tmpData = cviDataSyn[["aStart","aEnd"]]
#bGenome.extend(["cvi"]*len(tmpData))
#dataList = dataList.append(tmpData)
##
##tmpData = kyoDataSyn[["aStart","aEnd"]]
##bGenome.extend(["kyo"]*len(tmpData))
##dataList = dataList.append(tmpData)
#
#dataList["bGenome"] = bGenome
#
#dataList.sort_values(by = ["aStart","aEnd"], inplace = True)
#dat
#
#regions = getConservedRegions(dataList)




###########################################################
## Make regions using only Syn blocks
###########################################################
colDataList = getDataList([lerDataSyn, cviDataSyn],["ler","cvi"],["a","a"])
colRegions = getConservedRegions(colDataList)

lerDataList = getDataList([lerDataSyn,lerCviDataSyn],["col","cvi"],["b","a"])
lerRegions = getConservedRegions(lerDataList)

cviDataList = getDataList([cviDataSyn, lerCviDataSyn],["col","ler"],["b","b"])
cviRegions = getConservedRegions(cviDataList)

###########################################################
## Make regions using ALL blocks
###########################################################
colDataList = getDataList([lerData, cviData],["ler","cvi"],["a","a"])
colRegions = getConservedRegions(colDataList, isSyn = False)

lerDataList = getDataList([lerData,lerCviData],["col","cvi"],["b","a"])
lerRegions = getConservedRegions(lerDataList, isSyn = False)

cviDataList = getDataList([cviData, lerCviData],["col","ler"],["b","b"])
cviRegions = getConservedRegions(cviDataList, isSyn = False)





genomeData = {"col":colDataList,"ler": lerDataList, "cvi" :cviDataList}
genomeRegions = {"col":colRegions, "ler":lerRegions, "cvi":cviRegions}

#regG = Graph().as_undirected()
#regG.add_vertices(sum(map(len,genomeRegions.values())))

genNames = ["col","ler","cvi"]

adders = [0]
for i in genNames[:-1]:
    adders.append(sum(adders) + len(genomeRegions[i]))
adders = dict(zip(genNames,adders))



"""
for i in genNames:
    for region in genomeRegions[i].itertuples():
        regData = genomeData[i][(genomeData[i].start <= region.start) & (genomeData[i].end >= region.end)]
#        indices = np.intersect1d(np.where(genomeData[i].start <= region[0])[0], np.where(genomeData[i].end >= region[1])[0])
        for row in regData.itertuples(index=True):
            aStart = row.start
            aEnd = row.end
            bStart = row.bStart
            bEnd = row.bEnd
            bGenome = row.bGenome

            aStartPro = (region.start- aStart)/(aEnd - aStart)
            aEndPro = (region.end- aStart)/(aEnd - aStart)
            
            bRegStart = bStart + (aStartPro*(bEnd - bStart))
            bRegEnd = bStart + (aEndPro*(bEnd - bStart))
            
            length = bRegEnd - bRegStart
            bRegStart -= 0.01*length
            bRegEnd += 0.01*length
            
            matReg = genomeRegions[bGenome][(genomeRegions[bGenome].start >= bRegStart) & (genomeRegions[bGenome].end <= bRegEnd)]
            if len(matReg) > 0:
#                print(list(zip([adders[i] +region.Index]*len(matReg),adders[bGenome]+matReg.index.values)))
                regG.add_edges(zip([adders[i] +region.Index]*len(matReg),adders[bGenome]+matReg.index.values))
"""                   
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

clqs = getCLQ(regGMatrix, list(map(lambda x: len(genomeRegions[x]), genNames)))


coords = []
for i in clqs:
    clqCoords = []
    for j in range(len(i)):
        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].start)
        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].end)
    coords.append(clqCoords)
        
        
allData = pd.DataFrame(coords)
    

clqs_filtered = []
for i in clqs:
    count = 0
    for j in range(len(i)):
        if genomeRegions[genNames[j]].iloc[i[j]].end - genomeRegions[genNames[j]].iloc[i[j]].start > 100:
            count+=1
    if count == 3:
        clqs_filtered.append(i)
    
coords = []
for i in clqs_filtered:
    clqCoords = []
    for j in range(len(i)):
        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].start)
        clqCoords.append(genomeRegions[genNames[j]].iloc[i[j]].end)
    coords.append(clqCoords)    

allData_filtered = pd.DataFrame(coords)

sizeData = []
for i in clqs:
    _d = []
    for j in range(len(genNames)):
        _d.append(genomeRegions[genNames[j]].iloc[i[j]].end - genomeRegions[genNames[j]].iloc[i[j]].start)
    sizeData.append(_d)
        
    
transData = []
for i in allDataAll.itertuples(index = False):
    a = np.intersect1d(np.where(allData[0] <= i[0])[0],np.where(allData[1] >= i[1])[0])
    b = np.intersect1d(np.where(allData[2] <= i[2])[0],np.where(allData[3] >= i[3])[0])
    c = np.intersect1d(np.where(allData[4] <= i[4])[0],np.where(allData[5] >= i[5])[0])
    
#    if not len(intersect(a,b,c)) > 0:
    if len(a) < 1 or len(b) < 1 or len(c) < 1:
        transData.append(i)

transData = pd.DataFrame(transData)
coords = []
for i in transData.itertuples(index = False):
    if (i[1] - i[0] > 100) and (i[3] - i[2] > 100) and (i[5] - i[4] > 100):
      coords.append(i)
transData_filtered = pd.DataFrame(coords)   

    
commonBases = {}
for i in dirData.keys():
    if "col" in i:
        colStart = dirData[i].aStart
        colEnd = dirData[i].aEnd
        colBases = np.unique(unlist([list(range(colStart[i], colEnd[i]+1)) for i in range(len(colStart))]))
        commonBases[i] = colBases
        print(i, len(colBases))
        
keys = [i for i in dirData.keys() if "col" in i]

bases = commonBases[keys[0]]

for i in keys[1:]:
    bases = np.intersect1d(bases, commonBases[i])


colStartGenome = genomeData["col"].start
colEndGenome = genomeData["col"].end
#colBasesGenome = [list(range(colStartGenome[i], colEndGenome[i])) for i in range(len(colStartGenome))]
colBasesGenome = np.array([],dtype = "uint8")
count = 0
num = 0
for i in range(0,len(colStartGenome),100):
    print(i)
#    bases = []
    
    bases = np.unique(unlist([list(range(colStartGenome[j], colEndGenome[j])) for j in range(i,min(i+100, len(colStartGenome)))]))
    colBasesGenome = np.append(colBasesGenome,np.setdiff1d(bases, colBasesGenome))
#    colBasesGenome = np.append(colBasesGenome, np.setdiff1d(range(colStartGenome[i], colEndGenome[i]), colBasesGenome))
#    colBasesGenome.append(list(range(colStartGenome[i], colEndGenome[i])))
    if count == 10:
        count = 0
        num+=1
        colBasesGenome = unlist(colBasesGenome)
        print(num)
    else:
        count+=1

