import logging
import numpy as np
import pandas as pd

def readSAMBAM(fin, type='B'):
    import pysam
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        if type == 'B':
            findata = pysam.AlignmentFile(fin,'rb')
        elif type == 'S':
            findata = pysam.AlignmentFile(fin,'r')
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
                logger.error("Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(aln.cigarstring))
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
        coords = pd.DataFrame.from_dict(coords, orient='index')
        coords.sort_values([9, 0, 1, 2, 3, 10], inplace=True, ascending=True)
        return coords
    except Exception as e:
        logger.error("Error in reading BAM/SAM file. " + str(e))
        sys.exit()

finpath = '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/syri_test_genomes/analysis2/out.sam'

finpath = '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/col_ler_chr1.sam'
finpath1 = '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/col_ler_allchr.sam'

coords = readSAMBAM(finpath)
coords1 = readSAMBAM(finpath1)



lenCut = 0
idenCut = 0

minimapToTSV(finpath, lenCut, idenCut)



from collections import deque

def samtocoords(f):
    from pandas import DataFrame
    from collections import deque
    al = deque()
    with open(f, 'r') as fin:
        for l in fin:
            if l[0] == '@':
                continue
            l = l.split('\t')
            # TODO: Add check for non-mapping sequences
            # TODO Add check for correct CIGAR string
            cgt = [[int(j[0]), j[1]] for j in [i.split(';') for i in l[5].replace('S', ';S,').replace('H', ';H,').replace('=', ';=,').replace('X', ';X,').replace('I', ';I,').replace('D', ';D,').split(',')[:-1]]]

            bf = '{:012b}'.format(int(l[1]))

            rs = int(l[3])
            re = rs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

            # print(bf)
            if bf[7] == '0':
                if cgt[0][1] == '=':
                    qs = 1
                elif cgt[0][1] in ['S', 'H']:
                    qs = cgt[0][0] + 1
                else:
                    print('ERROR: CIGAR string starting with non-matching base')
                qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
            elif bf[7] == '1':
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
                re-rs+1,
                qe-qs+1,
                format((sum([i[0] for i in cgt if i[1] == '=']) / sum(
                    [i[0] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])) * 100, '.2f'),
                1,
                1 if bf[7] == '0' else -1,
                l[2],
                l[0],
                "".join([str(i[0])+i[1] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])
            ])
    al = DataFrame(al)
    al[6] = al[6].astype('float')
    # al = al.loc[al[6]>90]
    # al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] + al.loc[al[8] == -1, 3]
    # al.loc[al[8] == -1, 3] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    # al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    return al


fins = ['/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/ampril2/col_ler/out_r1.sam',
        '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/ampril2/col_ler/out_r2.sam',
        '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/ampril2/col_ler/out_r3.sam',
        '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/ampril2/col_ler/out_r4.sam',
        '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/ampril2/col_ler/out_r5.sam',
        '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/ampril2/col_ler/out_r6.sam',
        '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/ampril2/col_ler/out_r7.sam']

al = samtocoords(fins[0])
cs = deque()
for f in fins:
    cs.append(samtocoords(f))
