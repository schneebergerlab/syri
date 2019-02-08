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
from Bio.Alphabet import generic_dna
from re import findall

cimport numpy as np
cimport cython

np.random.seed(1)

from syri.pyxFiles.synsearchFunctions import readCoords

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
    logger = logging.getLogger("ShV")
    cwdpath = args.dir
    prefix = args.prefix


    allAlignments = readSRData(cwdpath, prefix, args.all)
    allAlignments["id"] = allAlignments.group.astype("str") + allAlignments.aChr + allAlignments.bChr + allAlignments.state
    allBlocks = pd.unique(allAlignments.id)

    if not args.cigar:
        logger.debug("finding short variation using MUMmer alignments")
        nc = args.nCores

        ## Could filter SNPs based no show-snp 'buff', but removed this feature
        # buff = args.buff
        buff = 0
        sspath = args.sspath
        delta = args.delta.name

        if delta not in os.listdir(cwdpath):
            logger.error("Delta file missing")
            sys.exit()

        blocklists = [allBlocks[_i:(_i+nc)] for _i in range(0, len(allBlocks), nc)]

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

    else:
        logger.debug("finding short variation using CIGAR string")
        coordsfin = args.infile.name
        chrmatch = args.chrmatch
        coords = readCoords(coordsfin, chrmatch, cigar=True)
        refg = {fasta.id:fasta.seq for fasta in parse(args.ref.name, 'fasta', generic_dna)}
        qryg = {fasta.id:fasta.seq for fasta in parse(args.qry.name, 'fasta', generic_dna)}

        with open('snps.txt', 'w') as fout:
            for b in allBlocks:
                block = allAlignments.loc[allAlignments.id == b].copy()
                fout.write("\t".join(["#",
                      str(block[["aStart", "aEnd"]].min().min()),
                      str(block[["aStart", "aEnd"]].max().max()),
                      str(block[["bStart", "bEnd"]].min().min()),
                      str(block[["bStart", "bEnd"]].max().max()),
                      pd.unique(block["aChr"])[0],
                      pd.unique(block["bChr"])[0]])+ "\n")
                for row in block.itertuples(index=False):
                    if not 'inv' in row.id:
                        cg = coords.loc[(coords.aStart == row.aStart) &
                                        (coords.aEnd == row.aEnd) &
                                        (coords.bStart == row.bStart) &
                        (coords.bEnd == row.bEnd) &
                        (coords.aChr == row.aChr) &
                        (coords.bChr == row.bChr), 'cigar']
                        brks = findall("(\d+)([IDX=])?", cg.iloc[0])

                        # chech for validity of CIGAR string
                        cuts = np.unique([i[1] for i in brks])
                        for i in cuts:
                            if i not in ['X','=','I','D']:
                                logger.error('Invalid CIGAR string. Only (X/=/I/D) operators are allowed')
                                sys.exit()

                        refseq = refg[row.aChr][(row.aStart-1):row.aEnd]
                        qryseq = qryg[row.bChr][(row.bStart-1):row.bEnd]

                        posa = 0                    # number of bases covered in genome a
                        posb = 0                    # number of bases covered in genome b
                        for i in brks:
                            if i[1] == '=':
                                posa += int(i[0])
                                posb += int(i[0])
                            elif i[1] == 'X':
                                for j in range(int(i[0])):
                                    out = [row.aStart+posa+j, refseq[posa+j], qryseq[posb+j], row.bStart+posb+j, 0, 0, 0, 0, 1, 1, row.aChr, row.bChr]
                                    fout.write("\t".join(list(map(str, out))) + '\n')
                                posa += int(i[0])
                                posb += int(i[0])
                            elif i[1] == 'I':
                                for j in range(int(i[0])):
                                    out = [row.aStart+posa-1, '.', qryseq[posb+j], row.bStart+posb+j, 0, 0, 0, 0, 1, 1, row.aChr, row.bChr]
                                    fout.write("\t".join(list(map(str, out))) + '\n')
                                posb += int(i[0])
                            elif i[1] == 'D':
                                for j in range(int(i[0])):
                                    out = [row.aStart+posa+j, refseq[posa+j], '.', row.bStart+posb-1, 0, 0, 0, 0, 1, 1,row.aChr, row.bChr]
                                    fout.write("\t".join(list(map(str, out))) + '\n')
                                posa += int(i[0])
                    else:
                        cg = coords.loc[(coords.aStart == row.aStart) &
                                        (coords.aEnd == row.aEnd) &
                                        (coords.bStart == row.bStart) &
                                        (coords.bEnd == row.bEnd) &
                                        (coords.aChr == row.aChr) &
                                        (coords.bChr == row.bChr), 'cigar']
                        brks = findall("(\d+)([IDX=])?", cg.iloc[0])

                        # chech for validity of CIGAR string
                        cuts = np.unique([i[1] for i in brks])
                        for i in cuts:
                            if i not in ['X','=','I','D']:
                                logger.error('Invalid CIGAR string. Only (X/=/I/D) operators are allowed')
                                sys.exit()

                        refseq = refg[row.aChr][(row.aStart-1):row.aEnd]
                        qryseq = qryg[row.bChr][(row.bEnd-1):row.bStart].reverse_complement()
                    # with open('snps.txt', 'w') as fout:
                        posa = 0                    # number of bases covered in genome a
                        posb = 0                    # number of bases covered in genome b
                        for i in brks:
                            if i[1] == '=':
                                posa += int(i[0])
                                posb += int(i[0])
                            elif i[1] == 'X':
                                for j in range(int(i[0])):
                                    out = [row.aStart+posa+j, refseq[posa+j], qryseq[posb+j], row.bStart-posb-j, 0, 0, 0, 0, 1, 1, row.aChr, row.bChr]
                                    fout.write("\t".join(list(map(str, out))) + '\n')
                                posa += int(i[0])
                                posb += int(i[0])
                            elif i[1] == 'I':
                                for j in range(int(i[0])):
                                    out = [row.aStart+posa-1, '.', qryseq[posb+j], row.bStart-posb-j, 0, 0, 0, 0, 1, 1, row.aChr, row.bChr]
                                    fout.write("\t".join(list(map(str, out))) + '\n')
                                posb += int(i[0])
                            elif i[1] == 'D':
                                for j in range(int(i[0])):
                                    out = [row.aStart+posa+j, refseq[posa+j], '.', row.bStart-posb+1, 0, 0, 0, 0, 1, 1,row.aChr, row.bChr]
                                    fout.write("\t".join(list(map(str, out))) + '\n')
                                posa += int(i[0])

        return

def minimapToTSV(finpath, lenCut, idenCut):
    '''
    This function transforms mimimap2 output file generated with parameters --eqx -cx asm* to a .tsv format suitable to
    be used by SyRI.
    '''
    with open(finpath, 'r') as fin:
        with open('coords.txt', 'w') as fout:
            for line in fin:
                line = line.strip().split()
                # add one to convert to 1-based positioning
                line = [int(line[7])+1, int(line[8]), int(line[2])+1, int(line[3]), int(line[8])-int(line[7]), int(line[3])-int(line[2]), format((int(line[9])/int(line[10]))*100, '.2f'), 1, 1 if line[4] == '+' else -1, line[5], line[0], line[-1].split(":")[-1]]
                if line[4] > lenCut and line[5] > lenCut and float(line[6]) > idenCut:
                    fout.write("\t".join(list(map(str, line))) + "\n")
