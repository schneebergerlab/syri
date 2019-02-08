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
from syri.pyxFiles.function cimport getAllLongestPaths, getConnectivityGraph

cimport numpy as np
cimport cython

np.random.seed(1)

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

