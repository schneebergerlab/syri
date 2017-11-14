# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%%
import pandas as pd
from collections import Counter
import numpy as np
from synsearchFunctions import *
from scipy.stats import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import *
from synsearchFunctions import *
import os
import sys
from SynSearch import SynSearch
from glob import glob



cwdPath = os.getcwd()+"/"
queryPath = cwdPath + sys.argv[1]
coordsFilePath = cwdPath + sys.argv[2]

scafSizeCutOff = 50000
alignSizeCutOff = 10000
nsLen = 500

print(cwdPath+"*.gp")

gpFileList = glob(cwdPath+"*.gp")
if len(gpFileList) == 0:
    sys.exit("NO mummerplot output file found")
else:
    print("Found info for generating",len(gpFileList),"chromosomes")
    
    
queryGenome = {fasta.id:fasta for fasta in SeqIO.parse(queryPath,'fasta')}
scafSize = {fasta.id: len(fasta.seq) for fasta in queryGenome.values()}
longScaf = getValues(list(scafSize.keys()),np.where(np.array(list(scafSize.values())) > scafSizeCutOff)[0])

coords = pd.read_table(coordsFilePath,header=None)
coordsFilter = coords.loc[coords[4] > alignSizeCutOff]
coordsFilter = coordsFilter.loc[coordsFilter[10].isin(longScaf)]
scafCountDict, scafSizeDict, scafChrDict = getscafDict(np.unique(coordsFilter[10]), coordsFilter, scafSize)


uniChrom = np.unique(coordsFilter[9])

NNs = "N"*nsLen

chrGpDict = {i:gp for i in uniChrom for gp in gpFileList if i in gp}

pseudoGenome = []
anchorInfo = []
for chromo in sorted(chrGpDict.keys()):
    gpFilePath = chrGpDict[chromo]
    chromoLen = 0
#    p=Popen(['mummerplot',"-f","-l","-r",chromo,"-p",chromo,"--postscript",'out_1_i90_l10000.delta'],cwd=cwdPath,stderr=PIPE,stdout=PIPE)
#    p.wait()
#    pOut = p.communicate()
#    print(pOut)
#    gpFilePath = cwdPath+chromo+".gp"
    orderedScafID, invertedID = orderFromMummerplot(gpFilePath)
    
    scaffoldList = []
    for i in orderedScafID:
        if i in scafChrDict.keys():
            if scafChrDict[i] == chromo:
                scaffoldList.append(i)
    chromSeq = []
    
    for i in scaffoldList:
        if i in invertedID:
            chromSeq.append(str(queryGenome[i].seq.reverse_complement()))
            anchorInfo.append([chromo,chromoLen+1,chromoLen+scafSize[i],i,scafSize[i],"-"])
        else:
            chromSeq.append(str(queryGenome[i].seq))
            anchorInfo.append([chromo,chromoLen+1,chromoLen+scafSize[i],i,scafSize[i],"+"])
        chromoLen = chromoLen + scafSize[i] + nsLen

    pseudoGenome.append(SeqRecord(seq=Seq((NNs).join(chromSeq)), id = chromo, description=""))

anchorInfo = pd.DataFrame(anchorInfo,columns = ["chrom","start","end","contig","length","orientation"])
anchorInfo.to_csv(cwdPath+"anchorInfo.txt", sep = "\t")
pseudoPath = cwdPath+"pseudoGenome.fasta"
SeqIO.write(pseudoGenome,pseudoPath,"fasta")


