#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 14:17:12 2018

script to find the annotation of a region from syri's output

@author: goel
"""
import os
import argparse
import sys

def readData(args, cwdPath):
    annoCoords = pd.DataFrame()
    for f in args.f:
        fin = f.name
        try:
            findf = pd.read_table(cwdPath+fin, header=None, dtype=object)
        except pd.io.common.EmptyDataError:
            print(fin, "Out.txt is empty. Skipping analysing it.")
            continue
        except Exception as e:
            print("ERROR: while trying to read ", fin, "Out.txt", e)
            continue
        annoData = findf.loc[findf[0] == "#"].copy().drop([0, 4], axis=1).dropna(axis=1, how="all")
        annoData["state"] = fin
        annoData = annoData.reindex(columns=["state"]+annoData.columns[:-1].tolist())
        # annoCoords = annoCoords.append(annoData.copy(), sort=True)
        annoCoords = pd.concat([annoCoords, annoData], sort=True)
    annoCoords[[2, 3, 6, 7]] = annoCoords[[2, 3, 6, 7]].astype("int")
    colnames = ["SVtype", "refChr", "refStart", "refEnd", "qryChr", "qryStart", "qryEnd", "ctxType", "Genome"]
    annoCoords.columns = colnames[:annoCoords.shape[1]]
    annoCoords.sort_values(["refChr", "refStart", "refEnd", "qryChr", "qryStart", "qryEnd"], inplace=True)
    annoCoords.index = range(len(annoCoords))
    return annoCoords
# END


def getRegion(args, cwdPath):
    data = readData(args, cwdPath)
    if not args.q:
        regions = data.loc[(data.refChr == args.chr) & (data.refStart <= args.end) & (data.refEnd >= args.start)].copy()
    elif args.q:
        regions = data.loc[(data.qryChr == args.chr) & (data.qryStart <= args.end) & (data.qryEnd >= args.start)].copy()
        regions.sort_values(["qryChr", "qryStart", "qryEnd", "refChr", "refStart", "refEnd"], inplace=True)
    print(regions.to_csv(sep="\t", index=False, header=False))
# END


def readsyriout(f, vars):
    """
    Reads syri output and select SR annotations
    """
    from collections import OrderedDict, deque
    import logging
    import pandas as pd
    import numpy as np
    # Reads syri.out. Select: achr, astart, aend, bchr, bstart, bend, srtype
    logger = logging.getLogger("readsyriout")
    syri_regs = deque()
    skipvartype = ['CPG', 'CPL', 'DEL', 'DUPAL', 'HDR', 'INS', 'INVAL', 'INVDPAL', 'INVTRAL', 'NOTAL', 'SNP', 'TDM', 'TRANSAL']
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            # TODO: DECIDE WHETHER TO HAVE STATIC VARS OR FLEXIBLE ANNOTATION
            if l[10] in vars:
                syri_regs.append(l)
            else:
                if l[10] not in skipvartype:
                    skipvartype.append(l[10])
                    logger.warning("{} is not a valid annotation for alignments in file {}. Alignments should belong to the following classes {}. Skipping alignment.".format(l[10], f, vars))
    try:
        df = pd.DataFrame(list(syri_regs))[[0, 1, 2, 5, 6, 7, 10]]
    except KeyError:
        raise ImportError("Incomplete input file {}, syri.out file should have 11 columns.".format(f))
    df[[0, 5, 10]] = df[[0, 5, 10]].astype(str)
    try:
        df.loc[df[0] != '-', [1, 2]] = df.loc[df[0] != '-', [1, 2]].astype(int)
        df.loc[df[5] != '-', [6, 7]] = df.loc[df[5] != '-', [6, 7]].astype(int)
    except ValueError:
        raise ValueError("Non-numerical values used as genome coordinates in {}. Exiting".format(f))
    # chr ID map
    # chrid = []
    chrid_dict = OrderedDict()
    for i in np.unique(df[0]):
        if i == '-': continue
        # chrid.append((i, np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]))
        chrid_dict[i] = np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]
    df.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type']
    return df, chrid_dict
# END


def getregion2(args, cwdPath):
    syriout = args.syriout.name
    # anno = args.anno
    vtypes = ['SYN', 'SYNAL', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']
    vtypes = vtypes + args.anno if args.anno is not None else vtypes
    data, chrid_dict = readsyriout(syriout, vtypes)
    if not args.q:
        data = data.loc[data.achr != '-']
        data.loc[:, ['astart', 'aend']] = data.loc[:, ['astart', 'aend']].astype(int)
        regions = data.loc[(data.achr == args.chr) & (data.astart <= args.end) & (data.aend >= args.start)].copy()
    elif args.q:
        data = data.loc[data.bchr != '-']
        data.loc[:, ['bstart', 'bend']] = data.loc[:, ['bstart', 'bend']].astype(int)
        regions = data.loc[(data.bchr == args.chr) & (data.bstart <= args.end) & (data.bend >= args.start)].copy()
        regions.sort_values(["bchr", "bstart", "bend", "achr", "astart", "aend"], inplace=True)
    print(regions.to_csv(sep="\t", index=False, header=False))
# END


def getBed(args, cwdPath):
    data = readData(args, cwdPath)
    mar = args.m
    if not args.q:
        outD = data[["refChr", "refStart", "refEnd"]].copy()
        outD.refStart = outD.refStart - 1 - mar
        outD.loc[outD.refStart < 0, "refStart"] = 0
        outD.refEnd = outD.refEnd - 1 + mar
    elif args.q:
        outD = data[["qryChr", "qryStart", "qryEnd"]].copy()
        outD.qryStart = outD.qryStart - 1 - mar
        outD.loc[outD.qryStart < 0, "qryStart"] = 0
        outD.qryEnd = outD.qryEnd - 1 + mar
    if args.a:
        outD[["SVType", "ctxType", "genome"]] = data[["SVtype", "ctxType", "Genome"]]
    print(outD.to_csv(sep="\t", index=False, header=False))
# END

    
def getSnpAnno(args, cwdPath):
    fin = pd.read_table(args.fin, header=None)
    fin.columns = ["chr", "pos"]
    fin.sort_values(["chr", "pos"], inplace=True)
    data = readData(args, cwdPath)
    outData = deque()
    chroms = deque()
    poses = deque()
    if not args.q:
        chromo = ""
        for row in fin.itertuples():
            if chromo != row.chr:
                chrData = data.loc[(data.refChr==row.chr)].copy()
                chromo = row.chr
            regions = chrData.loc[(chrData.refStart <= row.pos) & (chrData.refEnd >= row.pos)].copy()
            chroms.extend([row.chr]*len(regions))
            poses.extend([row.pos]*len(regions))
            outData.append(outData.append(regions))
    elif args.q:
        chromo = ""
        for row in fin.itertuples():
            if chromo != row.chr:
                chrData = data.loc[(data.qryChr==row.chr)].copy()
                chromo = row.chr
            regions = chrData.loc[(data.qryStart <= row.pos) & (data.qryEnd >= row.pos)].copy()
            chroms.extend([row.chr]*len(regions))
            poses.extend([row.pos]*len(regions))
            outData.append(outData.append(regions))
    outData = pd.concat(list(outData))
    outData["inChr"] = list(chroms)
    outData["inPos"] = list(poses)
    print(outData.to_csv(sep="\t", index=False, header=False), end="")
# END


def syri2bedpe(args):
    input = args.input.name
    output = args.output.name
    f = args.f
    # TODO: Complete this function
    with open(input, 'r') as fin, open(output, 'w') as fout:
        for line in fin:
            line = line.strip().split()
# END

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    
    parser_region = subparsers.add_parser("region", help="get annotation for the region")
    parser_getbed = subparsers.add_parser("getbed", help="get bed file for the regions")
    parser_snpanno = subparsers.add_parser("snpanno", help="get annotation for SNPs list")
    parser_syri2bedpe = subparsers.add_parser("syri2bedpe", help="Converts the syri.out to BEDPE format")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    parser_syri2bedpe.set_defaults(func=syri2bedpe)
    parser_syri2bedpe.add_argument('input', help="syri.out file", type=argparse.FileType('r'))
    parser_syri2bedpe.add_argument('output', help="output BEDPE file", type=argparse.FileType('w'))
    parser_syri2bedpe.add_argument('-f', help="Annotation to filter out. Use multiple times to filter more than one annotation types.", type=str, action='append')


    parser_region.set_defaults(func=getregion2)
    parser_region.add_argument("chr", help="Chromosome ID", type=str)
    parser_region.add_argument("start", help="region start", type=int)
    parser_region.add_argument("end", help="region end", type=int)
    parser_region.add_argument("syriout", help="syri output file", type=argparse.FileType('r'))
    parser_region.add_argument("--anno", help="Extra syri annotation to consider", action='append', type=str)
    parser_region.add_argument("-p", help="prefix", type=str)
    parser_region.add_argument("-q", help="search region in query genome", action="store_true", default=False)
    parser_region.add_argument("-f", help="files to search", type=argparse.FileType('r'), nargs="+", default=None)
    
    parser_getbed.set_defaults(func=getBed)
    parser_getbed.add_argument("-q", help="search region in query genome", action="store_true", default=False)
    parser_getbed.add_argument("-f", help="files to search", type=argparse.FileType('r'), nargs="+", default=None)
    parser_getbed.add_argument("-m", help="margin to add on the regions", type=int, default=0)
    parser_getbed.add_argument("-a", help="output region annotation", action="store_true", default=False)
    
    parser_snpanno.set_defaults(func=getSnpAnno)
    parser_snpanno.add_argument("fin", help="input file in tab separated format. Column 1: Chr. Column 2: position", type=argparse.FileType('r'), default=None)
    parser_snpanno.add_argument("-q", help="search annotation in query genome", action="store_true", default=False)
    parser_snpanno.add_argument("-f", help="files to search. by default all files will be searched", type=argparse.FileType('r'), nargs="+", default=None)

    cwdpath = os.getcwd()+os.sep
    args = parser.parse_args()
    if args.func != getregion2:
        if "f" in args:
            if args.f == None:
                args.f = list(map(open, ['synOut.txt', 'invOut.txt', 'TLOut.txt', 'invTLOut.txt', 'dupOut.txt', 'invDupOut.txt', 'ctxOut.txt']))

    # import pandas as pd
    # from collections import deque
    args.func(args, cwdpath)

        
    
    
        
