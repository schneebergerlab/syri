# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:36:01 2017

@author: goel
"""
import numpy as np
import operator as op
from functools import reduce
from Bio import SeqIO
from os import remove

def unlist(nestedList):
    """Take a nested-list as input and return a 1d list of all elements in it"""
    
    outList = []
    for i in nestedList:
        if type(i) in (list, np.ndarray):
            outList.extend(unlist(i))
        else:
            outList.append(i)
    return(outList)

def getValues(l, index):
    """from list l get the values at indices specified by index"""
    return [l[i] for i in index]
			
def getColors(colorPalette, numOfCol):
	return([colorPalette(i/numOfCol) for i in range(numOfCol)])


def subList(lst1, lst2):
    return(list(map(op.sub,lst1, lst2)))
    
def intersect(*lists):
    return reduce(np.intersect1d,list(lists))
    
def extractSeq(filePath, seqID, start = 0, end = -1):
    querySeq = [fasta for fasta in SeqIO.parse(filePath,'fasta') if fasta.id == seqID][0]
    querySeq.seq = querySeq.seq[start:end+1]
    SeqIO.write(querySeq,seqID+"_"+str(start)+"_"+str(end)+".fasta","fasta")

def fileRemove(fName):
    try:
        remove(fName)
    except OSError as e:
        if e.errno != 2:    ## 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
            raise