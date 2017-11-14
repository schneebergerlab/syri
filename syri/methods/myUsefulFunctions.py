# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:36:01 2017

@author: goel
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import operator as op
from functools import reduce
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import importlib.util

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

def plotDensity(data):
    density = gaussian_kde(data)
    xs = np.linspace(min(data), max(data),1000)
    density.covariance_factor = lambda : .2
    density._compute_covariance()
    plt.plot(xs,density(xs))
    plt.show()

def subList(lst1, lst2):
    return(list(map(op.sub,lst1, lst2)))
    
def intersect(*lists):
    return reduce(np.intersect1d,list(lists))
    
def extractSeq(filePath, seqID, start = 0, end = -1):
    querySeq = [fasta for fasta in SeqIO.parse(filePath,'fasta') if fasta.id == seqID][0]
    print(type(start))
    querySeq.seq = querySeq.seq[start:end+1]
    SeqIO.write(querySeq,seqID+"_"+str(start)+"_"+str(end)+".fasta","fasta")
#    return str(querySeq[seqID][start:end+1])

def loadModule(moduleName, modulePath):
    spec = importlib.util.spec_from_file_location(moduleName, modulePath)
    foo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(foo)
    return foo
    