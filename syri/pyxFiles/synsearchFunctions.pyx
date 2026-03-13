# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import numpy as np

from syri.scripts.func import *
import sys
from collections import deque, defaultdict
# from scipy.stats import *
import pandas as pd
from multiprocessing import Pool
from functools import partial
import os
from gc import collect
import logging
import psutil

from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map as cpp_map
from libcpp.deque cimport deque as cpp_deq
from syri.pyxFiles.function cimport getmeblocks, getOverlapWithSynBlocks
cimport numpy as np
cimport cython

np.random.seed(1)


def samtocoords(f):
    from pandas import DataFrame
    from collections import deque
    logger = logging.getLogger('SAM reader')
    rc = {}        # Referece chromosomes
    rcs = {}        # Selected chromosomes
    al = deque()    # Individual alignment
    try:
        with open(f, 'r') as fin:
            for l in fin:
                if l[:3] == '@SQ':
                    c, s = 0, 0
                    for h in l.strip().split()[1:]:
                        h = h.split(':')
                        if h[0] == 'SN': c = h[1]
                        if h[0] == 'LN': s = int(h[1])
                    rcs[c] = s
                    continue
                elif l[0] == '@': continue

                l = l.split('\t')[:6]
                # if l[1] == '2064': break
                if l[2] == '*':
                    logger.warning(l[0]+ ' do not align with any reference sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')  # Skip rows corresponding to non-mapping sequences (contigs/scaffolds)
                    continue

                if 'M' in l[5]:
                    logger.error(f'Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: {l[5]}. If using minimap2 for alignment, then use the --eqx parameter.')
                    sys.exit()
                cgt = [[int(j[0]), j[1]] for j in [i.split(';') for i in l[5].replace('S', ';S,').replace('H', ';H,').replace('=', ';=,').replace('X', ';X,').replace('I', ';I,').replace('D', ';D,').split(',')[:-1]]]
                if len(cgt) > 2:
                    if True in [True if i[1] in ['S', 'H'] else False for i in cgt[1:-1]]:
                        logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                        sys.exit()

                bf = '{:012b}'.format(int(l[1]))

                rs = int(l[3])
                re = rs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

                if bf[7] == '0':    # forward alignment
                    if cgt[0][1] == '=':
                        qs = 1
                    elif cgt[0][1] in ['S', 'H']:
                        qs = cgt[0][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                elif bf[7] == '1':  # inverted alignment
                    if cgt[-1][1] == '=':
                        qs = 1
                    elif cgt[-1][1] in ['S', 'H']:
                        qs = cgt[-1][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                    qs, qe = qe, qs

                al.append([
                    rs,
                    re,
                    qs,
                    qe,
                    abs(re-rs) + 1,
                    abs(qs-qe) + 1,
                    format((sum([i[0] for i in cgt if i[1] == '=']) / sum(
                        [i[0] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])) * 100, '.2f'),
                    1,
                    1 if bf[7] == '0' else -1,
                    l[2],
                    l[0],
                    "".join([str(i[0])+i[1] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])
                ])
                rcs[l[2]] = 1
            rcs = list(rcs.keys())
            for k in list(rc.keys()):
                if k not in rcs: logger.warning(l[0]+ ' do not align with any query sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')
    except Exception as e:
        logger.error('Error in reading SAM file: ' + str(e))
        sys.exit()
    al = DataFrame(list(al))
    al[6] = al[6].astype('float')
    al.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
    al.index = range(len(al.index))
    return al
# END

def readSAMBAM(fin, type='B'):
    import pysam
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        if type == 'B':
            findata = pysam.AlignmentFile(fin,'rb')
        elif type == 'S':
            return samtocoords(fin)
        else:
            raise ValueError("Wrong parameter")
    except ValueError as e:
        logger.error("Error in opening BAM/SAM file. " + str(e))
        sys.exit()
    except OSError as e:
        logger.error("Error in reading input file." + str(e))
        sys.exit()
    except Exception as e:
        logger.error("Unexpected error in opening BAM/SAM file. " + str(e))
        sys.exit()

    try:
        qry_prim = {}
        ref_prim = {}
        cgdict = {1:'I', 2:'D', 7:'=', 8:'X'}
        coords = {}
        index = 0
        for aln in findata:
            index += 1
            ## Check whether every sequence has at least one primary alignment
            if aln.reference_name is not None:
                if aln.reference_name not in ref_prim.keys():
                    ref_prim[aln.reference_name] = False
            if aln.query_name not in qry_prim.keys():
                qry_prim[aln.query_name] = False
            if aln.reference_name is not None:
                if not ref_prim[aln.reference_name]:
                    if aln.flag < 256:
                        ref_prim[aln.reference_name] = True
            if not qry_prim[aln.query_name]:
                if aln.flag < 256:
                    qry_prim[aln.query_name] = True

            ## Pass non-alinging chromosomes
            if aln.cigarstring is None:
                logger.warning(aln.query_name + ' do not align with any reference chromosome and cannot be analysed')
                continue

            ## Check CIGAR:
            if False in [False if i[0] not in [1,2,4,5,7,8] else True for i in aln.cigartuples]:
                logger.error(f'Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: {aln.cigarstring}. If using minimap2 for alignment, then use the --eqx parameter.')
                sys.exit()
            if len(aln.cigartuples) > 2:
                if True in [True if i[0] in [4,5] else False for i in aln.cigartuples[1:-1]]:
                    logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                    sys.exit()

            ## Parse information from the aln object
            astart = aln.reference_start+1
            aend = aln.reference_end
            is_inv = True if np.binary_repr(aln.flag,12)[7] == '1' else False
            if not is_inv:
                if aln.cigartuples[0][0] in [4,5]:
                    bstart = aln.cigartuples[0][1]+1
                else:
                    bstart = 1
                bend = bstart + aln.query_alignment_length - 1
            else:
                if aln.cigartuples[-1][0] in [4,5]:
                    bend = aln.cigartuples[-1][1]+1
                else:
                    bend = 1
                bstart = bend + aln.query_alignment_length - 1
            alen = abs(aend - astart) + 1
            blen = abs(bend - bstart) + 1
            iden = format((sum([i[1] for i in aln.cigartuples if i[0] == 7])/sum([i[1] for i in aln.cigartuples if i[0] in [1,2,7,8]]))*100, '.2f')
            adir = 1
            bdir = -1 if is_inv else 1
            achr = aln.reference_name
            bchr = aln.query_name
            cg = "".join([str(i[1]) + cgdict[i[0]] for i in aln.cigartuples if i[0] not in [4,5]])
            coords[index] = [astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg]

        ## Give warning for chromosomes which do not have any primary alignment
        for k,v in ref_prim.items():
            if not v:
                logger.warning('No primary alignment found for reference sequence ' + k +'. This could mean that the entire chromosome '+ k +' is reapeated.')
        for k,v in qry_prim.items():
            if not v:
                logger.warning('No primary alignment found for query sequence ' + k +'. This could mean that the entire chromosome '+ k + ' is reapeated.')

        ## Return alignments
        coords = pd.DataFrame.from_dict(coords, orient= 'index')
        coords.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except Exception as e:
        logger.error("Error in reading BAM/SAM file. " + str(e))
        sys.exit()
# END

def readPAF(paf):
    coords = deque()
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        with open(paf, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                astart = int(line[7]) + 1
                aend = int(line[8])
                adir = 1
                bdir = 1 if line[4] == '+' else -1
                bstart = int(line[2]) + 1 if bdir == 1 else int(line[3])
                bend = int(line[3]) if bdir == 1 else int(line[2]) + 1
                alen = abs(aend - astart) + 1
                blen = abs(bend - bstart) + 1 if bdir == 1 else bstart - bend + 1
                cg = [i.split(":")[-1] for i in line[12:] if i[:2] == 'cg']
                if len(cg) != 1:
                    logger.error("CIGAR string is not present in PAF at line {}. Exiting.".format("\t".join(line)))
                    sys.exit()
                cg = cg[0]
                ## Check CIGAR:
                if not all([True if i[1] in {'I', 'D', 'H', 'S', 'X', '='} else False for i in cgtpl(cg)]):
                    logger.error(f'Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: {cg}. If using minimap2 for alignment, then use the --eqx parameter.')
                    sys.exit()
                if len(cgtpl(cg)) > 2:
                    if any([True if i[1] in {'H', 'S'} else False for i in cgtpl(cg)]):
                        logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + str(cg))
                        sys.exit()

                iden = round((sum([int(i[0]) for i in cgtpl(cg) if i[1] == '='])/sum([int(i[0]) for i in cgtpl(cg) if i[1] in {'=', 'X', 'D', 'I'}]))*100, 2)
                achr = line[5]
                bchr = line[0]
                coords.append([astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg])
        coords = pd.DataFrame(coords)
        coords.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except FileNotFoundError:
        logger.error("Cannot open {} file. Exiting".format(paf))
        sys.exit()
    except ValueError as e:
        logger.error("Error in reading PAF: {}. Exiting".format(e))
        sys.exit()
# END

def readCoords(coordsfin, chrmatch, cwdpath, prefix, args, cigar = False):
    logger = logging.getLogger('Reading Coords')
    logger.debug(args.ftype)
    chrlink = {}
    if args.ftype == 'T':
        logger.info("Reading input from .tsv file")
        try:
            coords = pd.read_table(coordsfin, header = None)
        except pd.errors.ParserError:
            coords = pd.read_table(coordsfin, header = None, engine = "python")
        except Exception as e:
            logger.error("Error in reading the alignment file. " + e)
            sys.exit()
    elif args.ftype == 'S':
        logger.info("Reading input from SAM file")
        try:
            coords = readSAMBAM(coordsfin, type='S')
        except Exception as e:
            logger.error("Error in reading the alignment file. " + e)
            sys.exit()
    elif args.ftype == 'B':
        logger.info("Reading input from BAM file")
        try:
            coords = readSAMBAM(coordsfin, type='B')
        except Exception as e:
            logger.error("Error in reading the alignment file" + e)
            sys.exit()
    elif args.ftype == 'P':
        logger.info("Reading input from PAF file")
        try:
            coords = readPAF(coordsfin)
        except Exception as e:
            logger.error("Error in reading the alignment file" + e)
            sys.exit()
    else:
        logger.error("Incorrect alignment file type specified.")
        sys.exit()

    if not cigar:
        if coords.shape[1] >= 12:
            coords = coords.iloc[:, 0:11]
        coords.columns = ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr"]
    else:
        if coords.shape[1] > 12:
            coords = coords.iloc[:, 0:12]
        coords.columns = ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']

    # Sanity check input file
    try:
        coords.aStart = coords.aStart.astype('int')
    except ValueError:
        logger.error('astart is not int')
        sys.exit()

    try:
        coords.aEnd = coords.aEnd.astype('int')
    except ValueError:
        logger.error('aend is not int')
        sys.exit()

    try:
        coords.bStart = coords.bStart.astype('int')
    except ValueError:
        logger.error('bstart is not int')
        sys.exit()

    try:
        coords.bEnd = coords.bEnd.astype('int')
    except ValueError:
        logger.error('abend is not int')
        sys.exit()

    try:
        coords.aLen = coords.aLen.astype('int')
    except ValueError:
        logger.error('alen is not int')
        sys.exit()

    try:
        coords.bLen = coords.bLen.astype('int')
    except ValueError:
        logger.error('blen is not int')
        sys.exit()

    try:
        coords.iden = coords.iden.astype('float')
    except ValueError:
        logger.error('iden is not float')
        sys.exit()

    try:
        coords.aDir = coords.aDir.astype('int')
    except ValueError:
        logger.error('aDir is not int')
        sys.exit()

    if any(coords.aDir != 1):
        logger.error('aDir can only have values 1')
        sys.exit()

    try:
        coords.bDir = coords.bDir.astype('int')
    except ValueError:
        logger.error('bDir is not int')
        sys.exit()

    for i in coords.bDir:
        if i not in [1,-1]:
            logger.error('bDir can only have values 1/-1')
            sys.exit()

    try:
        coords.aChr = coords.aChr.astype(str)
    except:
        logger.error('aChr is not string')
        sys.exit()

    try:
        coords.bChr = coords.bChr.astype(str)
    except:
        logger.error('bChr is not string')
        sys.exit()

    # Filter small alignments
    if args.f:
        logger.info('Filtering low-quality alignments (alignment quality < 90, alignment length < 100)')
        logger.debug('Number of alignments before filtering: {}'.format(coords.shape[0]))
        coords = coords.loc[coords.iden > 90]
        coords = coords.loc[(coords.aLen>100) & (coords.bLen>100)]
        logger.debug('Number of alignments after filtering: {}'.format(coords.shape[0]))

    ## check for bstart > bend when bdir is -1
    check = np.unique(coords.loc[coords.bDir == -1, 'bStart'] > coords.loc[coords.bDir == -1, 'bEnd'])
    if len(check) > 1:
        logger.error('Inconsistent start and end position for inverted alignment in query genome. For inverted alignments, either all bstart < bend or all bend > bstart')
        sys.exit()
    elif len(check) == 0:
        logger.info('No Inverted alignments present.')
    elif check[0] == True:
        pass
    else:
        logger.info('For inverted alignments, bstart was less than bend. Swapping them.')
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] + coords.loc[coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bEnd'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[coords.bDir == -1, 'bEnd']

    coords.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)

    ## Ensure that chromosome IDs are same for the two genomes.
    ## Either find best query match for every reference genome.
    ## Or if --no-chrmatch is set then remove non-matching chromosomes.
    if np.unique(coords.aChr).tolist() != np.unique(coords.bChr).tolist():
        logger.warning('Chromosomes IDs do not match.')
        if not chrmatch:
            if len(np.unique(coords.aChr)) != len(np.unique(coords.bChr)):
                logger.error("Unequal number of chromosomes in the genomes. Exiting")
                sys.exit()
            else:
                logger.warning("Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.")
                chromMaps = defaultdict(dict)
                for i in np.unique(coords.bChr):
                    for j in np.unique(coords.aChr):
                        a = np.array(coords.loc[(coords.bChr == i) & (coords.aChr == j), ["aStart", "aEnd"]])
                        a = mergeRanges(a)
                        chromMaps[j][i] = len(a) + (a[:, 1] - a[:, 0]).sum()

                assigned = []
                fout = open(cwdpath+prefix+"mapids.txt", "w")
                for chrom in np.unique(coords.aChr):
                    maxid = max(chromMaps[chrom].items(), key=lambda x: x[1])[0]
                    if maxid in assigned:
                        logger.error("{} in genome B is best match for two chromosomes in genome A. Cannot assign chromosomes automatically.".format(maxid))
                        fout.close()
                        fileRemove(cwdpath+prefix+"mapids.txt")
                        sys.exit()
                    assigned.append(maxid)
                    fout.write(chrom+"\t"+maxid+"\n")
                    logger.info("setting {} as {}".format(maxid, chrom))
                    coords.loc[coords.bChr == maxid, "bChr"] = chrom
                    chrlink[maxid] = chrom
                fout.close()
        else:
            logger.warning("--no-chrmatch is set. Not matching chromosomes automatically.")
            aChromo = set(coords["aChr"])
            bChromo = set(coords["bChr"])
            badChromo = list(aChromo - bChromo) + list(bChromo - aChromo)
            if len(badChromo) > 0:
                logger.warning(", ".join(badChromo) + " present in only one genome. Removing corresponding alignments")
            coords = coords.loc[~coords.aChr.isin(badChromo) & ~coords.bChr.isin(badChromo)]

    ## Check for presence of directed alignments
    achrs = np.unique(coords.aChr).tolist()
    for achr in achrs:
        if coords.loc[(coords.aChr==achr) & (coords.bChr==achr) & (coords.bDir == 1),].shape[0] == 0:
            hombchr = [k for k,v in chrlink.items() if v==achr]
            if len(hombchr) == 1:
                hombchr = hombchr[0]
            elif len(hombchr) == 0:
                hombchr = achr
            else:
                logger.error('Homologous chromosomes were not identified correctly. Try assigning the chromosome ids manually.')
                sys.exit()
            logger.warning('Reference chromosome ' + achr + ' do not have any directed alignments with its homologous chromosome in the query genome (' + hombchr + '). Filtering out all corresponding alignments.')
            coords = coords.loc[~(coords.aChr == achr)]
            coords = coords.loc[~(coords.bChr == achr)]

    ## Check for presence of too many inverted alignments
    for achr in achrs:
        dir_range = mergeRanges(np.array(coords.loc[(coords.aChr==achr) & (coords.bChr==achr) & (coords.bDir==1), ["aStart", "aEnd"]]))
        dir_len = len(dir_range) + (dir_range[:, 1] - dir_range[:, 0]).sum()
        inv_range = mergeRanges(np.array(coords.loc[(coords.aChr==achr) & (coords.bChr==achr) & (coords.bDir==-1), ["aStart", "aEnd"]]))
        inv_len = len(inv_range) + (inv_range[:, 1] - inv_range[:, 0]).sum()
        if inv_len > dir_len:
            hombchr = [k for k,v in chrlink.items() if v==achr]
            if len(hombchr) == 1:
                hombchr = hombchr[0]
            elif len(hombchr) == 0:
                hombchr = achr
            else:
                logger.error('Homologous chromosomes were not identified correctly. Try assigning the chromosome ids manually.')
                sys.exit()
            logger.warning('Reference chromosome ' + achr + ' has high fraction of inverted alignments with its homologous chromosome in the query genome (' + hombchr + '). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.')
    return coords, chrlink


def startSyri(args, coords):
    nCores = args.nCores
    bRT = args.bruteRunTime
    threshold = 50  ##args.threshold
    cwdPath = args.dir
    prefix = args.prefix
    tUC = args.TransUniCount
    tUP = args.TransUniPercent
    invgl = args.invgl
    tdgl = args.tdgl
    tdolp = args.tdolp

    logger = logging.getLogger("syri")
    logger.info("starting")
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    uniChromo = list(np.unique(coords.aChr))
    logger.info('Analysing chromosomes: {}'.format(uniChromo))
    # Identify intra-chromosomal events (synteny, inversions, intra-trans, intra-dup) for each chromosome as a separate
    # process in parallel
    with Pool(processes = nCores) as pool:
        p = pool.map(partial(syri,threshold=threshold,coords=coords, cwdPath= cwdPath, bRT = bRT, prefix = prefix, tUC=tUC, tUP=tUP, invgl=invgl, tdgl=tdgl, tdolp=tdolp), uniChromo)
    if p != [None]*len(uniChromo):
        sys.exit()
    # for chromo in uniChromo:
    #     print(chromo)
    #     syri(chromo,threshold=threshold,coords=coords, cwdPath= cwdPath, bRT = bRT, prefix = prefix, tUC=tUC, tUP=tUP, invgl=invgl, tdgl=tdgl, tdolp=tdolp)

    # Merge output of all chromosomes
    mergeOutputFiles(uniChromo,cwdPath, prefix)

    #Identify cross-chromosomal events in all chromosomes simultaneously
    from syri.tdfunc import getCTX
    getCTX(coords, cwdPath, uniChromo, threshold, bRT, prefix, tUC, tUP, nCores, tdgl, tdolp)

    # Recalculate syntenic blocks by considering the blocks introduced by CX events
    outSyn(cwdPath, threshold, prefix)
    return 'Finished'


def syri(chromo, threshold, coords, cwdPath, bRT, prefix, tUC, tUP, invgl, tdgl, tdolp):
    logger = logging.getLogger("syri."+chromo)
    coordsData = coords[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == 1)]
    logger.info(chromo+" " + str(coordsData.shape))
    logger.info("Identifying Synteny for chromosome " + chromo)
    df = apply_TS(coordsData.aStart.values,coordsData.aEnd.values,coordsData.bStart.values,coordsData.bEnd.values, threshold)
    blocks = [alignmentBlock(i, df[i], coordsData.iloc[i]) for i in df.keys()]
    for block in blocks:
        i = 0
        while(i < len(block.children)):
            block.children = list(set(block.children) - set(blocks[block.children[i]].children))
            i+=1
        block.children.sort()

        for child in block.children:
            blocks[child].addParent(block.id)

        scores = [blocks[parent].score for parent in block.parents]
        if len(scores) > 0:
            block.bestParent(block.parents[scores.index(max(scores))], max(scores))
    synPath = getSynPath(blocks)
    synData = coordsData.iloc[synPath].copy()
    del(coordsData, blocks, df)
    collect()

    ##########################################################################
    #   Finding Inversions
    ##########################################################################
    logger.info("Identifying Inversions for chromosome " + chromo)

    from syri.inversions import getInversions
    invertedCoordsOri, profitable, bestInvPath, invData, synInInv, badSyn = getInversions(coords,chromo, threshold, synData, tUC, tUP, invgl)

    ##########################################################
    #### Identify Translocation and duplications
    ##########################################################
    logger.info("Identifying translocation and duplication for chromosome " + chromo)

    # Import functions
    from syri.tdfunc import blocksdata, makeTransGroupList, transBlock, getBestClusterSubset, getTransClasses, getDupGenome, getTransCluster

    chromBlocks = coords[(coords.aChr == chromo) & (coords.bChr == chromo)]
    inPlaceIndices = sorted(list(synData.index.values) + list(invData.index.values))
    inPlaceBlocks = chromBlocks[chromBlocks.index.isin(sorted(list(synData.index.values)))].copy()

    for i in bestInvPath:
        invPos = profitable[i].invPos
        invBlockData = invertedCoordsOri.iloc[invPos]
        invCoord = [invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]]
        invCoord.append(invCoord[1] - invCoord[0])
        invCoord.append(invCoord[3] - invCoord[2])
        invCoord.append(sum((invBlockData.aLen+invBlockData.bLen)*invBlockData.iden)/(invCoord[-2] + invCoord[-1]))
        invCoord.extend([1,-1,chromo,chromo])
        for j in range(profitable[i].neighbours[0]+1, profitable[i].neighbours[1]):
            inPlaceBlocks = inPlaceBlocks[inPlaceBlocks.index != synData.iloc[j].name]
            try:
                inPlaceIndices.remove(synData.iloc[j].name)
            except:
                pass
        # inPlaceBlocks = inPlaceBlocks.append(pd.Series(invCoord, index = inPlaceBlocks.columns, name = invPos[0]))
        inPlaceBlocks = pd.concat([inPlaceBlocks, pd.DataFrame(dict(zip(inPlaceBlocks.columns, invCoord)), index=[0])])

    inPlaceBlocks.sort_values(["aChr","aStart","aEnd","bChr","bStart","bEnd"], inplace = True)
    inPlaceBlocks.index = range(inPlaceBlocks.shape[0])
    outPlaceBlocks = chromBlocks[~chromBlocks.index.isin(inPlaceIndices)]

    logger.debug("Translocations : Filtered inplace and outplace alignments" + chromo)

    ## Create connectivity tree for directed and inverted blocks
    ## find all translocations which don't have large gaps between its alignments
    ## and are not overlappign with the syntenic blocks
    ## Merge directed and inverted blocks

    transBlocks, invTransBlocks, allTransBlocks, allTransIndexOrder = blocksdata(outPlaceBlocks, inPlaceBlocks, threshold, tUC, tUP, chromo, tdgl)

    logger.debug("Translocations : found blocks" + chromo)


    allTransGenomeAGroups = makeTransGroupList(allTransBlocks, "aStart", "aEnd", threshold)
    allTransGenomeBGroups = makeTransGroupList(allTransBlocks, "bStart", "bEnd", threshold)

    allTransGroupIndices = {}
    for i in range(len(allTransGenomeAGroups)):
        for block in allTransGenomeAGroups[i].member:
            allTransGroupIndices[block] = [i]
    for i in range(len(allTransGenomeBGroups)):
        for block in allTransGenomeBGroups[i].member:
            allTransGroupIndices[block].append(i)

    logger.debug("Translocations : getting clusters " + chromo)
    allTransCluster = getTransCluster(allTransGroupIndices, {i:allTransGenomeAGroups[i].member for i in range(len(allTransGenomeAGroups))}, {i:allTransGenomeBGroups[i].member for i in range(len(allTransGenomeBGroups))})

    allTransClusterIndices = dict()
    for i in range(len(allTransCluster)):
        allTransClusterIndices.update(dict.fromkeys(allTransCluster[i], i))

    logger.debug("Translocations : making blocks data " + chromo)
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    if len(allTransBlocks) > 0:
        auni = getOverlapWithSynBlocks(np.array(allTransBlocks.aStart),
                                       np.array(allTransBlocks.aEnd),
                                       np.array([chromo]*allTransBlocks.shape[0]),
                                       np.array(inPlaceBlocks.aStart),
                                       np.array(inPlaceBlocks.aEnd),
                                       np.array([chromo]*inPlaceBlocks.shape[0]),
                                       threshold,
                                       allTransBlocks.shape[0],
                                       tUC,
                                       tUP)
        sortedInPlace = inPlaceBlocks.sort_values(["bStart","bEnd"])
        buni = getOverlapWithSynBlocks(np.array(allTransBlocks.bStart), np.array(allTransBlocks.bEnd), np.array([chromo]*allTransBlocks.shape[0]), np.array(sortedInPlace.bStart), np.array(sortedInPlace.bEnd), np.array([chromo]*inPlaceBlocks.shape[0]), threshold, allTransBlocks.shape[0], tUC, tUP)

    genomeGroupLengths = ([len(i.member) for i in allTransGenomeAGroups], [len(i.member) for i in allTransGenomeBGroups])

    allTransBlocksData = [transBlock(i) for i in range(allTransBlocks.shape[0])]
    count = 0
    for row in allTransBlocks.itertuples(index = False):
        allTransBlocksData[count].aStart = row.aStart
        allTransBlocksData[count].aEnd = row.aEnd
        allTransBlocksData[count].bStart = row.bStart
        allTransBlocksData[count].bEnd = row.bEnd
        allTransBlocksData[count].dir = row.dir
        allTransBlocksData[count].transClusterIndex = allTransClusterIndices[count]
        allTransBlocksData[count].transGroupIndices = allTransGroupIndices[count]
        allTransBlocksData[count].aUni = auni[count]
        allTransBlocksData[count].bUni = buni[count]
        if genomeGroupLengths[0][allTransBlocksData[count].transGroupIndices[0]] == 1:
            allTransBlocksData[count].genomeAUni = True
        if genomeGroupLengths[1][allTransBlocksData[count].transGroupIndices[1]] == 1:
            allTransBlocksData[count].genomeBUni = True
        if (allTransBlocksData[count].aUni and allTransBlocksData[count].genomeAUni) or (allTransBlocksData[count].bUni and allTransBlocksData[count].genomeBUni):
            allTransBlocksData[count].setStatus(1)
        count+=1

    logger.debug("Translocations : finished making blocks data on" + chromo)
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    aUni = np.array([allTransBlocksData[i].aUni for i in range(allTransBlocks.shape[0])], dtype="int")
    bUni = np.array([allTransBlocksData[i].bUni for i in range(allTransBlocks.shape[0])], dtype="int")
    status = np.array([allTransBlocksData[i].status for i in range(allTransBlocks.shape[0])], dtype="int")
    aIndex = np.array([allTransBlocksData[i].transGroupIndices[0] for i in range(allTransBlocks.shape[0])], dtype="int")
    bIndex = np.array([allTransBlocksData[i].transGroupIndices[1] for i in range(allTransBlocks.shape[0])], dtype="int")

    ## get sorted values. sorted based on genome coordinate in allTransBlocks
    aGroups = {}
    for i in range(len(allTransGenomeAGroups)):
        aGroups[i] = allTransBlocks.iloc[allTransGenomeAGroups[i].member].sort_values(['aStart','aEnd']).index.values
    bGroups = {}
    for i in range(len(allTransGenomeBGroups)):
        bGroups[i] = allTransBlocks.iloc[allTransGenomeBGroups[i].member].sort_values(['bStart', 'bEnd']).index.values
    clstrsize = np.array([len(allTransCluster[i.transClusterIndex]) for i in allTransBlocksData], dtype = 'int')

    if len(allTransBlocks) > 0:
        out = getmeblocks(np.array(allTransBlocks.aStart),
                          np.array(allTransBlocks.aEnd),
                          np.array(allTransBlocks.bStart),
                          np.array(allTransBlocks.bEnd),
                          threshold,
                          aUni,
                          bUni,
                          status,
                          aIndex,
                          bIndex,
                          aGroups,
                          bGroups,
                          clstrsize,
                          tdolp)

        for i in range(len(out[0])):
            if out[0][i]:
                allTransCluster[allTransClusterIndices[i]].remove(i)

        for i in out[1].keys():
            if len(out[1][i]) > 0:
                allTransBlocksData[i].addMEBlock(list(out[1][i]))

        for i in out[2].keys():
            allTransBlocksData[i].setMEList(list(out[2][i][0]),list(out[2][i][1]))

        # del(aUni, bUni, status, aIndex, bIndex, aGroups, bGroups, out)
        # collect()

    logger.debug("Translocations : finding solutions "+ chromo)
    clusterSolutions = []
    for i in range(len(allTransCluster)):
        if len(allTransCluster[i]) > 0:
            if len(allTransCluster[i]) > 10000:
                clusterSolutions.append(getBestClusterSubset(allTransCluster[i], allTransBlocksData, bRT, tdolp, chromo, aGroups, bGroups, threshold))
            else:
                clusterSolutions.append(getBestClusterSubset(allTransCluster[i], allTransBlocksData, bRT, tdolp, chromo))

    clusterSolutionBlocks = [i[1] for i in clusterSolutions]
    #clusterBlocks = unlist(clusterSolutionBlocks)

    logger.debug("Translocations : processing translocations " + chromo)

    garb = deque()
    for i in range(len(allTransBlocksData)):
        if not allTransBlocksData[i].aUni and not allTransBlocksData[i].bUni:
            garb.append(0)
        elif allTransBlocksData[i].status == 1:
            garb.append(0)
        elif not allTransBlocksData[i].aUni:
            garb.append(1)
        elif not allTransBlocksData[i].bUni:
            garb.append(2)
        else:
            garb.append(3)
    meclass = np.array(list(garb), np.uint16)

    transClasses = getTransClasses(clusterSolutionBlocks,
                                   allTransBlocksData,
                                   allTransGenomeAGroups,
                                   allTransGenomeBGroups,
                                   allTransBlocks.aStart.values.astype(np.uint),
                                   allTransBlocks.aEnd.values.astype(np.uint),
                                   allTransBlocks.bStart.values.astype(np.uint),
                                   allTransBlocks.bEnd.values.astype(np.uint),
                                   aIndex,
                                   bIndex,
                                   aGroups,
                                   bGroups,
                                   threshold,
                                   meclass,
                                   tdolp)
    dupData = allTransBlocks.iloc[transClasses["duplication"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    invDupData = allTransBlocks.iloc[transClasses["invDuplication"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    TLData = allTransBlocks.iloc[transClasses["translocation"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])
    invTLData = allTransBlocks.iloc[transClasses["invTranslocation"]].sort_values(by = ["aStart","aEnd","bStart","bEnd"])

    dupData = getDupGenome(dupData,
                           allTransBlocksData,
                           transClasses,
                           allTransBlocks.aStart.values.astype(np.uint),
                           allTransBlocks.aEnd.values.astype(np.uint),
                           allTransBlocks.bStart.values.astype(np.uint),
                           allTransBlocks.bEnd.values.astype(np.uint),
                           aIndex,
                           bIndex,
                           aGroups,
                           bGroups,
                           threshold,
                           meclass,
                           tdolp)
    invDupData = getDupGenome(invDupData,
                              allTransBlocksData,
                              transClasses,
                              allTransBlocks.aStart.values.astype(np.uint),
                              allTransBlocks.aEnd.values.astype(np.uint),
                              allTransBlocks.bStart.values.astype(np.uint),
                              allTransBlocks.bEnd.values.astype(np.uint),
                              aIndex,
                              bIndex,
                              aGroups,
                              bGroups,
                              threshold,
                              meclass,
                              tdolp)


    fout = open(cwdPath+prefix+chromo+"_invOut.txt","w")
    tempInvBlocks = []
    for i in bestInvPath:
        invPos = profitable[i].invPos
        tempInvBlocks.append([invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2]])
        fout.write("\t".join(map(str,["#",invertedCoordsOri.iat[invPos[0],0],invertedCoordsOri.iat[invPos[-1],1],"-",invertedCoordsOri.iat[invPos[-1],3],invertedCoordsOri.iat[invPos[0],2],"\n"])))
        for j in invPos:
            fout.write("\t".join(map(str,invertedCoordsOri.iloc[j][:4])))
            fout.write("\n")
    fout.close()


    ## Grouping Syn blocks : Final synblock identification is done after ctx identification.
    allBlocks, outClusters = groupSyn(tempInvBlocks, dupData, invDupData, invTLData, TLData, threshold, synData, badSyn)
    if outClusters == [[]]:
        logger.error(f"No syntenic region found for chromosome: {chromo}. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.")
        return -1

    orderedBlocks = outPlaceBlocks[outPlaceBlocks.bDir == 1]
    invertedBlocks = outPlaceBlocks[outPlaceBlocks.bDir == -1]

########################################################################################################################
    with open(cwdPath+prefix+chromo+"_synOut.txt","w") as fout:
        for i in outClusters:
            fout.write("\t".join(map(str,["#",allBlocks.at[i[0],"aStart"],allBlocks.at[i[-1],"aEnd"],"-",allBlocks.at[i[0],"bStart"],allBlocks.at[i[-1],"bEnd"],"\n"])))
            for j in i:
                fout.write("\t".join(map(str,allBlocks.loc[j][:-1])))
                if j in synInInv:
                    fout.write("\tSyn_in_Inv\n")
                else:
                    fout.write("\n")
########################################################################################################################

    with open(cwdPath+prefix+chromo+"_dupOut.txt","w") as fout:
        for i in dupData.index.values:
            fout.write("\t".join(map(str,["#",dupData.at[i,"aStart"],dupData.at[i,"aEnd"],"-",dupData.at[i,"bStart"],dupData.at[i,"bEnd"],"-", dupData.at[i,"dupGenomes"],"\n"])))
            for j in transBlocks[allTransIndexOrder[i]]:
                fout.write("\t".join(map(str,orderedBlocks.iloc[j][:4])))
                fout.write("\n")

########################################################################################################################

    with open(cwdPath+prefix+chromo+"_invDupOut.txt","w") as fout:
        for i in invDupData.index.values:
            fout.write("\t".join(map(str,["#",invDupData.at[i,"aStart"],invDupData.at[i,"aEnd"],"-",invDupData.at[i,"bStart"],invDupData.at[i,"bEnd"],"-", invDupData.at[i,"dupGenomes"],"\n"])))
            for j in invTransBlocks[allTransIndexOrder[i]]:
                fout.write("\t".join(map(str,invertedBlocks.iloc[j][:4])))
                fout.write("\n")

########################################################################################################################

    with open(cwdPath+prefix+chromo+"_TLOut.txt","w") as fout:
        for i in TLData.index.values:
            fout.write("\t".join(map(str,["#",TLData.at[i,"aStart"],TLData.at[i,"aEnd"],"-",TLData.at[i,"bStart"],TLData.at[i,"bEnd"],"\n"])))
            for j in transBlocks[allTransIndexOrder[i]]:
                fout.write("\t".join(map(str,orderedBlocks.iloc[j][:4])))
                fout.write("\n")

########################################################################################################################

    with open(cwdPath+prefix+chromo+"_invTLOut.txt","w") as fout:
        for i in invTLData.index.values:
            fout.write("\t".join(map(str,["#",invTLData.at[i,"aStart"],invTLData.at[i,"aEnd"],"-",invTLData.at[i,"bStart"],invTLData.at[i,"bEnd"],"\n"])))
            for j in invTransBlocks[allTransIndexOrder[i]]:
                fout.write("\t".join(map(str,invertedBlocks.iloc[j][:4])))
                fout.write("\n")

    return
# END
########################################################################################################################


cpdef apply_TS(long[:] astart, long[:] aend, long[:] bstart, long[:] bend, int threshold, int mxgap = 100000000000):
    cdef:
        Py_ssize_t                              i, j,  n = len(astart)
        cpp_map[long, cpp_deq[long]]            df
        cpp_map[long, cpp_deq[long]].iterator   mapit
    for i in range(<Py_ssize_t> n):
        for j in range(<Py_ssize_t> i+1, <Py_ssize_t> n):
            if (astart[j] - aend[i]) < mxgap:        # Select only alignments with small gaps
                if (astart[j] - astart[i]) > threshold:
                    if (aend[j] - aend[i]) > threshold:
                        if (bstart[j] - bend[i]) < mxgap:        # Select only alignments with small gaps
                            if (bstart[j] - bstart[i]) > threshold:
                                if (bend[j] - bend[i]) > threshold:
                                    df[i].push_back(j)
    out = {}
    for i in range(n):
        if df.count(i)== 1:
            out[i] = [df[i][j] for j in range(<Py_ssize_t> df[i].size())]
        else:
            out[i] = []
    return out


def getSynPath(blocks):
    cdef list synPath = []
    scores = [block.score for block in blocks]
    cdef int lastBlock = scores.index(max(scores))
    while blocks[lastBlock].bestParentID != -1:
        synPath.append(lastBlock)
        lastBlock = blocks[lastBlock].bestParentID        
    synPath.append(lastBlock)
    return(synPath[::-1])


def outSyn(cwdPath, threshold, prefix):
    ctxAnnoDict = {"duplication":"dupCtx",
                   "invDuplication":"invDupCtx",
                   "translocation":"TLCtx",
                   "invTranslocation":"invTLCtx"}
    reCoords =  pd.DataFrame()

    synData = []
    with open(cwdPath+prefix+"synOut.txt","r") as fin:
        for line in fin:
            line = line.strip().split("\t")
            if line[0] == "#":
                chromo = line[1]
                continue
            if len(line) == 4:
                synData.append(list(map(int,line[:4]))+[chromo,chromo])
            elif len(line) == 5:
                synData.append(list(map(int,line[:4]))+[chromo,chromo] + [line[4]])

    synData = pd.DataFrame(synData)
    if len(synData.columns) == 6:
        synData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"]
    else:
        synData.columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr","isinInv"]
    synData["class"] = "syn"

    for i in ["invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt","ctxOut.txt"]:
        data = []
        with open(cwdPath+prefix+i,"r") as fin:
            if i != "ctxOut.txt":
                for line in fin:
                    line = line.strip().split("\t")
                    if line[0] == "#":
                        data.append(list(map(int,getValues(line,[2,3,6,7]))) + [line[1],line[5]])
                data = pd.DataFrame(data, columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr"], dtype=object)
                data["class"] = i.split("Out.txt")[0]
                if len(data)>0:
                    # reCoords = reCoords.append(data)
                    reCoords = pd.concat([reCoords, data])
            else:
                for line in fin:
                    line = line.strip().split("\t")
                    if line[0] == "#":
                        data.append(list(map(int,getValues(line,[2,3,6,7]))) + [line[1],line[5],ctxAnnoDict[line[8]]])
                data = pd.DataFrame(data, columns = ["aStart","aEnd","bStart","bEnd","aChr","bChr","class"], dtype=object)
                if len(data)>0:
                    # reCoords = reCoords.append(data)
                    reCoords = pd.concat([reCoords, data])

    # allBlocks = synData[["aStart","aEnd","bStart","bEnd","aChr","bChr","class"]].append(reCoords)
    allBlocks = pd.concat([synData[["aStart","aEnd","bStart","bEnd","aChr","bChr","class"]], reCoords])
    allBlocks.index = range(allBlocks.shape[0])
    allBlocks.sort_values(["aChr","aStart","aEnd","bChr","bStart","bEnd"], inplace= True)
    synLocs = {np.where(allBlocks.index.values == i)[0][0]:i for i in range(synData.shape[0])}

    allBlocks.index = range(allBlocks.shape[0])

    aClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():
        if len(currentCluster) == 0:
            if row["class"] != "syn":
                continue
            elif row["class"] == "syn":
                curChr = row["aChr"]
                currentCluster.append(index)
        elif row["class"] == "syn":
            if row["aChr"] == curChr:
                currentCluster.append(index)
            else:
                aClusters.append(currentCluster)
                currentCluster = [index]
                curChr = row["aChr"]

        elif row["class"] in ["TL", "inv","invTL","TLCtx","invTLCtx"]:
            aClusters.append(currentCluster)
            currentCluster = []
            curChr = ""
        else:
            if row["aEnd"] < allBlocks.loc[currentCluster[-1]]["aEnd"] + threshold:
                continue
            else:
                allClasses = allBlocks["class"][index:]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if max(row["aStart"],allBlocks.loc[currentCluster[-1]]["aEnd"]) > allBlocks.loc[nextSyn]["aStart"] - threshold:
                        continue
                    else:
                        aClusters.append(currentCluster)
                        currentCluster = []
                else:
                    aClusters.append(currentCluster)
                    currentCluster = []
    aClusters.append(currentCluster)

    allBlocks.sort_values(["bChr","bStart","bEnd","bChr","aStart","aEnd"],inplace = True)
    bClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():
        if len(currentCluster) == 0:
            if row["class"] != "syn":
                continue
            elif row["class"] == "syn":
                curChr = row["bChr"]
                currentCluster.append(index)
        elif row["class"] == "syn":
            if row["bChr"] == curChr:
                currentCluster.append(index)
            else:
                bClusters.append(currentCluster)
                currentCluster = [index]
                curChr = row["bChr"]
        elif row["class"] in ["TL", "inv","invTL","TLCtx","invTLCtx"]:
            bClusters.append(currentCluster)
            currentCluster = []
            curChr = ""
        else:
            if row["bEnd"] < allBlocks.loc[currentCluster[-1]]["bEnd"] + threshold:
                continue
            else:
                allClasses = allBlocks["class"][list(allBlocks.index.values).index(index):]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if max(row["bStart"], allBlocks.loc[currentCluster[-1]]["bEnd"]) > allBlocks.loc[nextSyn]["bStart"] - threshold:
                        continue
                    else:
                        bClusters.append(currentCluster)
                        currentCluster = []
                else:
                    bClusters.append(currentCluster)
                    currentCluster = []
    bClusters.append(currentCluster)
    allBlocks.sort_values(["aChr","aStart","aEnd","bChr", "bStart","bEnd"],inplace = True)

    outClusters = []
    aIndex = 0
    bIndex = 0
    currentCluster = []
    for i in unlist(aClusters):
        if i in aClusters[aIndex] and i in bClusters[bIndex]:
            currentCluster.append(i)
        else:
            if i not in aClusters[aIndex]:
                aIndex+=1
            if i not in bClusters[bIndex]:
                bIndex+=1
            outClusters.append(currentCluster)
            currentCluster = [i]
    outClusters.append(currentCluster)

    hasSynInInv = "isinInv" in synData.columns

    with open(cwdPath+prefix+"synOut.txt","w", encoding="utf-8") as fout:
        for i in outClusters:
            fout.write("\t".join(map(str,["#",allBlocks.at[i[0],"aChr"],allBlocks.at[i[0],"aStart"],allBlocks.at[i[-1],"aEnd"],"-",allBlocks.at[i[0],"aChr"],allBlocks.at[i[0],"bStart"],allBlocks.at[i[-1],"bEnd"]])) +"\n")
            for j in i:
                fout.write("\t".join(map(str,allBlocks.loc[j][0:4])))
                if hasSynInInv and synData.loc[synLocs[j]]["isinInv"] == "Syn_in_Inv":
                    fout.write("\tSyn_in_Inv\n")
                else:
                    fout.write("\n")
    return None

        
def groupSyn(tempInvBlocks, dupData, invDupData, invTLData, TLData, threshold, synData, badSyn):
    
    synData = synData.drop(synData.index.values[badSyn])
    allBlocks = synData[["aStart","aEnd","bStart","bEnd"]].copy()
    allBlocks["class"] = "syn"
    
    tempInvBlocks = pd.DataFrame(tempInvBlocks,columns =["aStart","aEnd","bStart","bEnd"], dtype= object)
    tempInvBlocks["class"] = "inv"
    
    tempDupData = dupData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempDupData["class"] = "dup"
    
    tempInvDupData = invDupData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempInvDupData["class"] = "invDup"
    
    tempInvTLData = invTLData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempInvTLData["class"] = "invTL"
    
    tempTLData = TLData[["aStart","aEnd","bStart","bEnd"]].copy()
    tempTLData["class"] = "TL"
    
    allBlocks = pd.concat([allBlocks,tempInvBlocks, tempInvDupData, tempInvTLData, tempTLData, tempDupData])
    allBlocks.index = range(allBlocks.shape[0])
    
    """
    Take data of all blocks and create groups of syntenic blocks from syntenic alignments
    """
    
    allBlocks.sort_values(["aStart","aEnd","bStart","bEnd"],inplace = True)
    
    aClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():        
        if len(currentCluster) == 0 and row["class"] != "syn":
            continue
        
        if row["class"] == "syn":
            currentCluster.append(index)
        elif row["class"] in ["TL", "inv","invTL"]:
            aClusters.append(currentCluster)
            currentCluster = []
        else:
            if row["aEnd"] < allBlocks.loc[currentCluster[-1]]["aEnd"] + threshold:
                continue
            else:
                # allClasses = allBlocks["class"][index:]
                allClasses = allBlocks.iloc[index:]["class"]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if row["aStart"] > allBlocks.loc[nextSyn]["aStart"] - threshold:
                        continue
                    else:
                        aClusters.append(currentCluster)
                        currentCluster = []
                else:
                    aClusters.append(currentCluster)
                    currentCluster = []
    aClusters.append(currentCluster)
    
    allBlocks.sort_values(["bStart","bEnd","aStart","aEnd"],inplace = True)
    
    bClusters = []
    currentCluster = []
    for index, row in allBlocks.iterrows():
        
        if len(currentCluster) == 0 and row["class"] != "syn":
            continue
        
        if row["class"] == "syn":
            currentCluster.append(index)
        elif row["class"] in ["TL", "inv","invTL"]:
            bClusters.append(currentCluster)
            currentCluster = []
        else:
            if row["bEnd"] < allBlocks.loc[currentCluster[-1]]["bEnd"] + threshold:
                continue
            else:
                # allClasses = allBlocks["class"][index:]
                allClasses = allBlocks.iloc[index:]["class"]
                if len(np.where(allClasses=="syn")[0]) > 0:
                    nextSyn = allClasses.index[np.where(allClasses=="syn")[0][0]]
                    if row["bStart"] > allBlocks.loc[nextSyn]["bStart"] - threshold:
                        continue
                    else:
                        bClusters.append(currentCluster)
                        currentCluster = []
                else:
                    bClusters.append(currentCluster)
                    currentCluster = []
    bClusters.append(currentCluster)
    allBlocks.sort_values(["aStart","aEnd","bStart","bEnd"],inplace = True)
    
    outClusters = []
    aIndex = 0 
    bIndex = 0
    currentCluster = []
    for i in range(synData.shape[0]):
        if i in aClusters[aIndex] and i in bClusters[bIndex]:
            currentCluster.append(i)
        else:
            if i not in aClusters[aIndex]:
                aIndex+=1
            if i not in bClusters[bIndex]:
                bIndex+=1
            outClusters.append(currentCluster)
            currentCluster = [i]
    outClusters.append(currentCluster)
    return (allBlocks, outClusters)

def mergeOutputFiles(uniChromo,path,prefix):
    def addData(fName,anno, chromo):
        fPath = open(path+prefix+chromo+"_"+anno+"Out.txt","r")
        for line in fPath.readlines():
            line = line.strip().split("\t")
            if line[0] == "#":
                fName.write("\t".join(unlist([line[0], chromo, line[1:4], chromo, line[4:]])) + "\n")
            else:
                fName.write("\t".join(line) + "\n")
        fPath.close()
        fileRemove(path+prefix+chromo+"_"+anno+"Out.txt")
                
    fSyn = open(path+prefix+"synOut.txt","w")
    fInv = open(path+prefix+"invOut.txt","w")
    fTL = open(path+prefix+"TLOut.txt","w")
    fInvTL = open(path+prefix+"invTLOut.txt","w")
    fDup = open(path+prefix+"dupOut.txt","w")
    fInvDup = open(path+prefix+"invDupOut.txt","w")
    
    files = [fSyn, fInv, fTL, fInvTL, fDup, fInvDup]
    classes = ["syn","inv","TL","invTL","dup","invDup"]
    
    for chromo in uniChromo:
        for i in range(len(classes)):
            addData(files[i], classes[i], chromo)
            
    for f in files:
        f.close()


class alignmentBlock:
    def __init__(self, id, children, data):
        self.id = id
        self.children = list(children)
        self.parents = []
        self.score = (data.aLen + data.bLen) * data.iden
        self.bestParentID = -1
    
    def addParent(self, parentID):
        self.parents.append(parentID)
        
    def bestParent(self,parentID,maxScore):
        self.bestParentID = parentID
        self.score = self.score + maxScore

