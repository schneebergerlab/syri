def maxsyn(coords, fix=False):
    '''
    Read alignments between chromosomes and decide chromosome order corresponding
    to maximum synteny between the two genomes.
    This is done by fusing the chromosomes to get larger ' putative ancestral chromosomes',
    while maximising syntenic regions.
    Optional: consider expanding it to multiple genomes
    :param coords: Alignments between the genomes (output of readCoords)
    :param fix: Fix the reference chromosomes. Only query chromosomes would be
    considered for fusing
    :return: Order of reference and query chromosome
    '''
    # from syri.synsearchFunctions import readCoords
    from hometools.classes import Namespace
    from collections import deque
    import logging
    import pandas as pd
    import numpy as np
    from itertools import product

    infin = '/srv/netscratch/dep_mercier/grp_schneeberger/projects/syri2/results/ancestral/thaliana_v_lyrata/out.paf'
    # coords = readCoords(infin, chrmatch, cwdpath, prefix, args, cigar = False)
    args = Namespace(ftype='P', f=True)
    coords = pd.read_table('/srv/netscratch/dep_mercier/grp_schneeberger/projects/syri2/results/ancestral/thaliana_v_lyrata/input_alignments.txt')
    achr = coords.aChr.unique()
    bchr = coords.bChr.unique()
    chrlink = deque()
    for a, b in product(achr[:5], bchr[:8]):
        chrcoords = coords.loc[(coords.aChr==a) & (coords.bChr==b)]
        fcoords = chrcoords.loc[chrcoords.bDir==1].copy()
        rcoords = chrcoords.loc[chrcoords.bDir==-1].copy()
        fcoords.iden = fcoords.iden.astype(float)
        fscore = sum((fcoords.aLen + fcoords.bLen) * fcoords.iden)/100000
        rscore = sum((rcoords.aLen + rcoords.bLen) * rcoords.iden)/100000
        try:
            chrlink.append([a, b, fscore/(fscore+rscore), rscore/(fscore+rscore), fscore, rscore])
        except ZeroDivisionError:
            chrlink.append([a, b, 0, 0, 0, 0])
        print(a, '\t', b)
        break



    return
# END