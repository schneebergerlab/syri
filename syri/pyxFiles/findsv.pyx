# cython: language_level = 3
# distutils: language = c++

import numpy as np
from syri.scripts.func import *
from igraph import *
from scipy.stats import *
import pandas as pd
import logging
from collections import deque
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
        except pd.errors.EmptyDataError:
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
        # annoCoords = annoCoords.append(coordsData.copy())
        annoCoords = pd.concat([annoCoords, coordsData.copy()])

    fileData = None

    try:
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header = None, dtype = object)
    except pd.errors.ParserError as e:
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header=None, dtype = object, engine ="python")
    except pd.errors.EmptyDataError:
        print("ctxOut.txt is empty. Skipping analysing it.")
    except Exception as e:
        print("ERROR: while trying to read ", fileType, "Out.txt", e)

    if fileData is not None:
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
        # annoCoords = annoCoords.append(coordsData1)
        annoCoords = pd.concat([annoCoords, coordsData1])

    annoCoords.columns = ["aStart","aEnd","bStart","bEnd","group","aChr","bChr","state"]
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))
    return annoCoords
# END


def getSV(cwdPath, allAlignments, prefix, offset):
    """
    inverted regions are output in reverse as well
    :param cwdPath:
    :param allAlignments:
    :param prefix:
    :param offset:
    :return: Output coordinate info:
        In reference: HDR/DEL coordinates consists of the affected bases, INS coordinate consists of 1bp upstream coordinate. This ensure that the inserted sequence (as in VCF) matches the query genome sequence
        In query: HDR/INS coordinates consists of the affected bases, DEL coordinate consists of 1bp upstream coordinate (towards the 5' end for directed alignment and towards the 3' end for inverted alignments)
    """
    logger = logging.getLogger("getSV")
    offset = -abs(offset)
    fout = open(cwdPath + prefix + "sv.txt", "w")
    allAlignments["id"] = allAlignments.group.astype(
        "str") + 'Chr' + allAlignments.aChr + 'Chr' + allAlignments.bChr + allAlignments.state
    allBlocks = pd.unique(allAlignments.id)

    for i in allBlocks:
        blocksAlign = allAlignments.loc[allAlignments.id == i].copy()
        if len(pd.unique(blocksAlign["aChr"])) > 1 or len(pd.unique(blocksAlign["aChr"])) > 1:
            logger.error("More than one chromosome found for a SR\n"+blocksAlign.to_string())
            sys.exit()
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

            if offset <= m <= 0:  ## No significant overlap in reference genome
                if offset <= n <= 0:  ## No significant overlap in query genome
                    continue
                elif n > 0:
                    # s = str(max(blocksAlign.iat[j, 1], blocksAlign.iat[j + 1, 0]))
                    # Use the base upstream to the insertion
                    s = str(blocksAlign.iat[j, 1])
                    if ordered:
                        fout.write("\t".join(["INS",
                                              s,
                                              s,
                                              str(blocksAlign.iat[j, 3] + 1),
                                              str(blocksAlign.iat[j + 1, 2] - 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        fout.write("\t".join(["INS",
                                              s,
                                              s,
                                              str(blocksAlign.iat[j, 3] - 1),
                                              str(blocksAlign.iat[j + 1, 2] + 1),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

                elif n < offset:
                    if ordered:
                        j_prop = abs(n) / (blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2])
                        j1_prop = abs(n) / (blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2])
                        sCoord = round(
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
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
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

            elif m == 1:
                if offset <= n <= 0:
                    # e = str(min(blocksAlign.iat[j, 3], blocksAlign.iat[j + 1, 2]))
                    # Using the base upstream to the deletion in the direction of alignment
                    e = str(blocksAlign.iat[j, 3])
                    if ordered:
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              e,
                                              e,
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        # e = str(max(blocksAlign.iat[j, 3], blocksAlign.iat[j + 1, 2]))
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              e,
                                              e,
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
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
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
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
                        fout.write("\t".join(["CPL",
                                              str(sCoord),
                                              str(eCoord),
                                              str(blocksAlign.iat[j + 1, 2]),
                                              str(blocksAlign.iat[j, 3]),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")

            elif m > 1:
                if offset <= n <= 0:
                    # e = str(min(blocksAlign.iat[j, 3], blocksAlign.iat[j + 1, 2]))
                    # Using the base upstream to the deletion in the direction of alignment
                    e = str(blocksAlign.iat[j, 3])
                    if ordered:
                        # e = str(min(blocksAlign.iat[j, 3], blocksAlign.iat[j + 1, 2]))
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              e,
                                              e,
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        # e = str(max(blocksAlign.iat[j, 3], blocksAlign.iat[j + 1, 2]))
                        fout.write("\t".join(["DEL",
                                              str(blocksAlign.iat[j, 1] + 1),
                                              str(blocksAlign.iat[j + 1, 0] - 1),
                                              e,
                                              e,
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
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
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
                            blocksAlign.iat[j, 1] - j_prop * (blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                        eCoord = round(blocksAlign.iat[j + 1, 0] + j1_prop * (
                                    blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
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
                            blocksAlign.iat[j, 3] - j_prop * (blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2]))
                        eCoord = round(blocksAlign.iat[j + 1, 2] + j1_prop * (
                                    blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2]))
                        fout.write("\t".join(["CPG",
                                              str(blocksAlign.iat[j + 1, 0]),
                                              str(blocksAlign.iat[j, 1]),
                                              str(sCoord),
                                              str(eCoord),
                                              blocksAlign.iat[0, 5],
                                              blocksAlign.iat[0, 6]]) + "\n")
                    else:
                        sCoord = round(blocksAlign.iat[j, 3] + j_prop * (blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3]))
                        eCoord = round(blocksAlign.iat[j + 1, 2] - j1_prop * (blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3]))
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
                                        blocksAlign.iat[j, 3] - blocksAlign.iat[j, 2]))
                            eCoord = round(blocksAlign.iat[j + 1, 2] + j1_prop * (
                                        blocksAlign.iat[j + 1, 3] - blocksAlign.iat[j + 1, 2]))
                            fout.write("\t".join(["TDM",
                                                  str(blocksAlign.iat[j + 1, 0]),
                                                  str(blocksAlign.iat[j, 1]),
                                                  str(sCoord),
                                                  str(eCoord),
                                                  blocksAlign.iat[0, 5],
                                                  blocksAlign.iat[0, 6]]) + "\n")
                        else:
                            sCoord = round(blocksAlign.iat[j, 3] + j_prop * (
                                        blocksAlign.iat[j, 2] - blocksAlign.iat[j, 3]))
                            eCoord = round(blocksAlign.iat[j + 1, 2] - j1_prop * (
                                        blocksAlign.iat[j + 1, 2] - blocksAlign.iat[j + 1, 3]))
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
                                        blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                            eCoord = round(blocksAlign.iat[j + 1, 0] + k1_prop * (
                                        blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
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
                                        blocksAlign.iat[j, 1] - blocksAlign.iat[j, 0]))
                            eCoord = round(blocksAlign.iat[j + 1, 0] + k1_prop * (
                                        blocksAlign.iat[j + 1, 1] - blocksAlign.iat[j + 1, 0]))
                            fout.write("\t".join(["TDM",
                                                  str(sCoord),
                                                  str(eCoord),
                                                  str(blocksAlign.iat[j + 1, 2]),
                                                  str(blocksAlign.iat[j, 3]),
                                                  blocksAlign.iat[0, 5],
                                                  blocksAlign.iat[0, 6]]) + "\n")
    fout.close()
    return None
# END


def addsvseq(svfin: str, refname: str, qryname: str, chrlink: dict):
    """
    Add sequence for indels and hdr
    :param svfin: Complete path to the sv file generated by getSV. This file will be overwritten.
    :param refname: Complete path to the reference genome.
    :param qryname: Complete path to the query genome.
    :param hdrseq: Set True if sequences for HDRs are also required.
    :return: SV sequence. For reference genome, only forward strand sequence is added. For query genome, reverse completed sequence is added when the SR is on the reverse stranded (INV, INVTR, INVDP). Upstream matching base is not added.
    """
    # vt = {'INS', 'DEL', 'HDR'} if hdrseq else {'INS', 'DEL'}
    vt = {'INS', 'DEL', 'HDR'}
    svdata = deque()
    # Read sv.txt generated using getSV
    with open(svfin, 'r') as fin:
        for line in fin:
            svdata.append(line.strip().split())
    gseq = readfasta(refname)
    # Add reference genome sequence
    for sv in svdata:
        if sv[0] not in vt:
            sv.append('-')
            continue
        if sv[0] == 'INS':
            sv.append('-')
            continue
        sv.append(gseq[sv[5]][(int(sv[1])-1):int(sv[2])])
    # Add query genome sequence
    gseq = readfasta(qryname)
    if len(chrlink) > 0:
        for k, v in chrlink.items():
            gseq[chrlink[k]] = gseq[k]
            gseq.pop(k)
    for sv in svdata:
        if sv[0] not in vt:
            sv.append('-')
            continue
        if sv[0] == 'DEL':
            sv.append('-')
            continue
        # Check whether SV corresponds to inversion
        if int(sv[3]) < int(sv[4]):
            sv.append(gseq[sv[6]][(int(sv[3]) - 1):int(sv[4])])
        else:
            # reverse complement sequence from SVs in inverted regions
            sv.append(revcomp(gseq[sv[6]][(int(sv[4])-1):int(sv[3])]))
    with open(svfin, 'w') as fout:
        for sv in svdata:
            fout.write("\t".join(sv) + '\n')
# END


def getNotAligned(cwdPath, prefix, ref, qry, chrlink):
    logger = logging.getLogger("getNA")

    refSize = {chrid: len(seq) for chrid, seq in readfasta(ref).items()}
    qrySize = {chrid: len(seq) for chrid, seq in readfasta(qry).items()}
    # qrySize = {fasta.id: len(fasta.seq) for fasta in parse(qry,'fasta')}

    annoCoords = pd.DataFrame()
    for fileType in ["syn","inv", "TL", "invTL","dup", "invDup"]:
        try:
            fileData = pd.read_table(cwdPath+prefix+fileType+"Out.txt", header=None, dtype = object)
        except pd.errors.ParserError as e:
            fileData = pd.read_table(cwdPath+prefix+fileType+"Out.txt", header=None, dtype = object, engine ="python")
        except pd.errors.EmptyDataError:
            print(fileType, "Out.txt is empty. Skipping analysing it.")
            continue
        except Exception as e:
            print("ERROR: while trying to read ", fileType, "Out.txt", e)
            continue
        coordsData = fileData.loc[fileData[0] == "#",[2,3,6,7,1,5]].copy()
        coordsData[[2,3,6,7]] = coordsData[[2,3,6,7]].astype(dtype="int64")
        coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
        # annoCoords = annoCoords.append(coordsData.copy())
        annoCoords = pd.concat([annoCoords, coordsData.copy()])

    try:
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header = None, dtype = object)
    except pd.errors.ParserError as e:
        fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header=None, dtype = object, engine ="python")
    except pd.errors.EmptyDataError:
        print("ctxOut.txt is empty. Skipping analysing it.")
    except Exception as e:
        print("ERROR: while trying to read ", fileType, "Out.txt", e)
#    fileData = pd.read_table(cwdPath+prefix+"ctxOut.txt", header = None, names = list(range(11)), dtype = object, sep ="\t")
    coordsData = fileData.loc[fileData[0] == "#"]
    coordsData = coordsData[[2,3,6,7,1,5]].copy()
    coordsData[[2,3,6,7]] = coordsData[[2,3,6,7]].astype(dtype="int64")
    coordsData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]

    # annoCoords = annoCoords.append(coordsData.copy())
    annoCoords = pd.concat([annoCoords, coordsData.copy()])
    annoCoords.sort_values(by = ["aChr", "aStart","aEnd","bChr", "bStart","bEnd"], inplace = True)
    annoCoords.index = range(len(annoCoords))

    for i in annoCoords.bChr:
        if i not in qrySize.keys():
            for k,v in chrlink.items():
                if v == i:
                    qrySize[i] = qrySize.pop(k)

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
# END
