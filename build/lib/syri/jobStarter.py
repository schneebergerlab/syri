#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 11:28:08 2017

@author: goel
"""

import pandas as pd
from subprocess import Popen


fileLocation = "/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/tests/nucmer_out/tair_ler/out_max_r.coords"
coords = pd.read_table(fileLocation, header = None) 

coords.columns  =   ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr"]
aChromo = set(coords["aChr"])
bChromo = set(coords["bChr"])
uniChromo = list(aChromo if len(aChromo) < len(bChromo) else bChromo)
uniChromo.sort()

# Set this as an argument
threshold = 50

for chromo in uniChromo:
    command = 'python -u  /netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/scripts/python/SynSearchPackage/SynSearch.py '+chromo
    args = ['bsub','-M','10000','-R',"'rusage[mem=100000]'","-n","2","-o","/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/scripts/python/SynSearchPackage/"+chromo+"_Synsearch_output.txt","-e","/netscratch/dep_coupland/grp_schneeberger/projects/SynSearch/scripts/python/SynSearchPackage/"+chromo+"_Synsearch_error.txt",command]
    Popen(args)
    