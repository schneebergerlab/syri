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
            finD = pd.read_table(cwdPath+fin, header = None, dtype = object)
        except pd.io.common.EmptyDataError:
            print(fin, "Out.txt is empty. Skipping analysing it.")
            continue
        except Exception as e:
            print("ERROR: while trying to read ", fin, "Out.txt", e)
            continue
            
        annoData = finD.loc[finD[0] == "#"].copy().drop( [0,4],axis = 1).dropna(axis = 1 , how="all")
        annoData["state"] = fin
        annoData = annoData.reindex(columns = ["state"]+annoData.columns[:-1].tolist())
        annoCoords = annoCoords.append(annoData.copy())
    annoCoords[[2,3,6,7]] = annoCoords[[2,3,6,7]].astype("int")
    colNames = ["SVtype","refChr","refStart","refEnd","qryChr","qryStart","qryEnd","ctxType","Genome"]
    annoCoords.columns = colNames[:annoCoords.shape[1]]
    annoCoords.sort_values(["refChr","refStart","refEnd","qryChr","qryStart","qryEnd"],inplace=True)
    annoCoords.index = range(len(annoCoords))
    return annoCoords

def getRegion(args, cwdPath):
    data = readData(args,cwdPath)
    if not args.q:
        regions = data.loc[(data.refChr==args.chr) & (data.refStart <= args.end) & (data.refEnd >= args.start)].copy()
    elif args.q:
        regions = data.loc[(data.qryChr==args.chr) & (data.qryStart <= args.end) & (data.qryEnd >= args.start)].copy()
        regions.sort_values(["qryChr","qryStart","qryEnd","refChr","refStart","refEnd"],inplace=True)
    print(regions.to_csv(sep="\t",index=False, header=False))

def getBed(args, cwdPath):
    data = readData(args,cwdPath)
    mar = args.m
    if not args.q:
        outD = data[["refChr","refStart","refEnd"]].copy()
        outD.refStart = outD.refStart - 1 - mar
        outD.loc[outD.refStart < 0, "refStart"] = 0
        outD.refEnd = outD.refEnd - 1 + mar
    elif args.q:
        outD = data[["qryChr","qryStart","qryEnd"]].copy()
        outD.qryStart = outD.qryStart - 1 - mar
        outD.loc[outD.qryStart < 0, "qryStart"] = 0
        outD.qryEnd = outD.qryEnd - 1 + mar
    if args.a:
        outD[["SVType","ctxType","genome"]] = data[["SVtype","ctxType","Genome"]]
    print(outD.to_csv(sep="\t", index=False, header=False))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    
    parser_region = subparsers.add_parser("region", help="get annotation for the region")
    parser_getbed = subparsers.add_parser("getbed", help= "get bed file for the regions")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    parser_region.set_defaults(func=getRegion)
    parser_region.add_argument("chr", help="Chromosome ID", type=str)
    parser_region.add_argument("start", help="region start", type=int)
    parser_region.add_argument("end", help = "region end", type=int)
    parser_region.add_argument("-q",help="search region in query genome",action="store_true", default=False)
    parser_region.add_argument("-f",help="files to search",type=argparse.FileType('r'), nargs="+", default = None )
    
    parser_getbed.set_defaults(func=getBed)
    parser_getbed.add_argument("-q",help="search region in query genome",action="store_true", default=False)
    parser_getbed.add_argument("-f", help="files to search",type=argparse.FileType('r'), nargs="+", default = None)
    parser_getbed.add_argument("-m",help="margin to add on the regions", type= int, default = 0)
    parser_getbed.add_argument("-a",help="output region annotation", action="store_true", default = False)
    
    cwdPath = os.getcwd()+os.sep
    args = parser.parse_args()
    if "f" in args:
        if args.f == None:
            args.f = list(map(open,['synOut.txt', 'invOut.txt', 'TLOut.txt', 'invTLOut.txt', 'dupOut.txt', 'invDupOut.txt', 'ctxOut.txt']))
    
    import pandas as pd
    args.func(args, cwdPath)
    
    
        
    
    
        