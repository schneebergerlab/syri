#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 14:33:37 2017

@author: goel
"""
import os
import sys
import numpy as np
import pandas as pd
from synsearchFunctions import *
from myUsefulFunctions import *


def readSVData(cwdPath, uniChromo, ctxout):
    annoCoords = pd.DataFrame()
    for chromo in uniChromo:
        for fileType in ["syn","inv", "TL", "invTL"]:
            fileData = pd.read_table(cwdPath+chromo+"_"+fileType+"Out.txt", header = None, dtype = object)
            annoIndices = np.where(fileData[0] =="#")[0]
            coordsData = fileData.loc[fileData[0] !="#"].copy()
            coordsData = coordsData[[0,1,2,3]].astype(dtype = "int64")
            annoIndices = np.append(annoIndices,len(fileData))
            repCount = annoIndices[1:] - annoIndices[:-1] - 1
            reps = np.repeat(range(len(annoIndices)-1), repCount)
            coordsData["group"] = reps
            coordsData["aChr"] = chromo
            coordsData["bChr"] = chromo
            coordsData["state"] = fileType
            annoCoords = annoCoords.append(coordsData.copy())
                                   
    fileData = pd.read_table(cwdPath+ctxout, header = None, names = list(range(11)), dtype = object, sep ="\t")
    annoIndices = np.where(fileData[0] =="#")[0]
    states = list(fileData[9].loc[annoIndices])
    coordsData = fileData.loc[fileData[0] !="#"].copy()
    annoIndices = np.append(annoIndices,len(fileData))
    repCount = annoIndices[1:] - annoIndices[:-1] - 1
    reps = np.repeat(range(len(annoIndices)-1), repCount)
    stateReps = np.repeat(states, repCount)
    
    coordsData1 = coordsData[[0,1,2,3]].astype(dtype = "int64")
    coordsData1["aChr"] = coordsData[9]
    coordsData1["bChr"] = coordsData[10]
    coordsData1["group"] = reps
    coordsData1["state"] = stateReps
    coordsData1 = coordsData1[[0,1,2,3,"group","aChr","bChr","state"]]
    coordsData1 = coordsData1.loc[coordsData1["state"].isin(["translocation","invTranslocation"])]
    coordsData1.loc[coordsData1.state == "translocation","state"] = "ctx"
    coordsData1.loc[coordsData1.state == "invTranslocation","state"] = "invCtx"
    annoCoords = annoCoords.append(coordsData1)
    annoCoords.columns = ["aStart","aEnd","bStart","bEnd","group","aChr","bChr","state"]
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))
    invIndices = annoCoords.state == "invCtx"
    annoCoords.loc[invIndices, "bStart"] = annoCoords.bStart + annoCoords.bEnd
    annoCoords.loc[invIndices, "bEnd"] = annoCoords.bStart - annoCoords.bEnd
    annoCoords.loc[invIndices, "bStart"] = annoCoords.bStart - annoCoords.bEnd
    return annoCoords
        
#def getmumSNPIn(allAlignments):
#    mumSNPIn = allAlignments.copy()
#    invIndices = mumSNPIn.state.isin(["invCtx"])
#    mumSNPIn.loc[invIndices,"bStart"] = mumSNPIn.bStart + mumSNPIn.bEnd
#    mumSNPIn.loc[invIndices,"bEnd"] = mumSNPIn.bStart - mumSNPIn.bEnd
#    mumSNPIn.loc[invIndices,"bStart"] = mumSNPIn.bStart - mumSNPIn.bEnd
#    
#    mumSNPIn = mumSNPIn[["aStart","aEnd","bStart","bEnd","aChr","bChr"]].copy()
#    return mumSNPIn
    

#%%
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



def getNotAligned(cwdPath, uniChromo, ctxout):    
    annoCoords = pd.DataFrame()
    for chromo in uniChromo:
        for fileType in ["syn","inv", "TL", "invTL","dup", "invDup"]:
            fileData = pd.read_table(cwdPath+chromo+"_"+fileType+"Out.txt", header = None, dtype = object)
            coordsData = fileData.loc[fileData[0] == "#"].copy()
            coordsData = coordsData[[1,2,4,5]].astype(dtype="int64")
            coordsData["aChr"] = chromo
            coordsData["bChr"] = chromo
            coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
            annoCoords = annoCoords.append(coordsData.copy())

    fileData = pd.read_table(cwdPath+ctxout, header = None, names = list(range(11)), dtype = object, sep ="\t")
    coordsData = fileData.loc[fileData[0] == "#"]
    coordsData = coordsData[[1,2,3,4,7,8]].copy()
    coordsData[[1,2,3,4]] = coordsData[[1,2,3,4]].astype(dtype="int64")
    coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
    
    annoCoords = annoCoords.append(coordsData.copy())
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))
  
    fout = open(cwdPath + "notAligned.txt","w")
    
    fout.write("#Reference Genome\n")
    df = annoCoords[["aStart","aEnd","aChr"]].copy()
    df.sort_values(["aChr", "aStart", "aEnd"], inplace = True)
    for chrom in uniChromo:
        chromData = df.loc[df.aChr == chrom]
        maxEnd = chromData.iloc[0,1]
        for row in chromData.itertuples(index = False):
            if row.aStart > maxEnd+1:
                fout.write("\t".join([str(maxEnd+1),
                                      str(row.aStart - 1),
                                      chrom]) + "\n")
            if row.aEnd > maxEnd:
                maxEnd = row.aEnd
    
    fout.write("\n")
    
    fout.write("#Query Genome\n")
    df = annoCoords[["bStart","bEnd","bChr"]].copy()
    df.sort_values(["bChr", "bStart", "bEnd"], inplace = True)
    for chrom in uniChromo:
        chromData = df.loc[df.bChr == chrom]
        maxEnd = chromData.iloc[0,1]
        for row in chromData.itertuples(index = False):
            if row.bStart > maxEnd+1:
                fout.write("\t".join([str(maxEnd+1),
                                      str(row.bStart - 1),
                                      chrom]) + "\n")
            if row.bEnd > maxEnd:
                maxEnd = row.bEnd
                
    fout.close()
    return None
    
     
             
         
#    overlaps = (np.array(df.aEnd[:-1]) - np.array(df.aStart[1:]) < -1)
#    sameChromo = np.array(df.aChr[:-1]) == np.array(df.aChr[1:])
#    outIndices = overlaps & sameChromo
    fout.write("#Reference Genome\n")
    for i in range(len(outIndices)):
        if outIndices[i] == True:
            fout.write("\t".join([str(df.iloc[i,1] + 1),
                                  str(df.iloc[i+1,0] - 1),
                                  df.iloc[i,2]]) + "\n")
    
    fout.write("\n")
    
    df = allAlignments[["bStart","bEnd","bChr"]].copy()
    invIndices = df.bStart > df.bEnd
    
    df.loc[invIndices, "bStart"] = df.loc[invIndices, "bStart"] + df.loc[invIndices, "bEnd"]
    df.loc[invIndices, "bEnd"] = df.loc[invIndices, "bStart"] - df.loc[invIndices, "bEnd"]
    df.loc[invIndices, "bStart"] = df.loc[invIndices, "bStart"] - df.loc[invIndices, "bEnd"]
    
    df.sort_values(["bChr", "bStart", "bEnd"], inplace = True)
    overlaps = (np.array(df.bEnd[:-1]) - np.array(df.bStart[1:]) < -1)
    sameChromo = np.array(df.bChr[:-1]) == np.array(df.bChr[1:])
    outIndices = overlaps & sameChromo
    fout.write("#Query Genome\n")
    for i in range(len(outIndices)):
        if outIndices[i] == True:
            fout.write("\t".join([str(df.iloc[i,1] + 1),
                                  str(df.iloc[i+1,0] - 1),
                                  df.iloc[i,2]]) + "\n")
    
    fout.close()
    return None
    
    
            
            
    
    
                        
#%%                


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python getSVs.py <path_to_coords_file>")
        sys.exit()
    if len(sys.argv) == 2:
        fileLocation = sys.argv[1]   
        fileLocation = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/col_ler_Chr/out_m_i90_l100.coords"
    
    cwdPath = os.getcwd()+"/"
    cwdPath = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/testRuns/col_ler_Chr/"
    coords = pd.read_table(fileLocation, header = None) 
    coords.columns  =   ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr"]
    aChromo = set(coords["aChr"])
    bChromo = set(coords["bChr"])
    uniChromo = list(aChromo if len(aChromo) < len(bChromo) else bChromo)
    uniChromo.sort()
    ctxout = "ctxOut.txt"
    
    if ctxout not in os.listdir(cwdPath):
        print("No ctx Out file in directory. Exiting.")
        sys.exit()
    allAlignments = readSVData(cwdPath, uniChromo, ctxout)
    mumSNPIn = allAlignments[["aStart","aEnd","bStart","bEnd","aChr","bChr"]].copy()
    mumSNPIn.to_csv(cwdPath+"mumSNPIn.txt",sep="\t",header = False, index = False)
    getSV(cwdPath, allAlignments)
    
    getNotAligned(cwdPath, uniChromo, ctxout)
