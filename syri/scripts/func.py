# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:36:01 2017
@author: goel
"""
import operator as op
from functools import reduce
from os import remove


def readfaidxbed(f):
    """
    Reads .faidx file from a genome assembly and returns a BED file consisting
    for entire chromosomes
    """
    from collections import deque
    import pybedtools as bt
    fabed = deque()
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            fabed.append([line[0], 1, int(line[1])])
    return list(fabed)
# END


def unlist(nestedList):
    """Take a nested-list as input and return a 1d list of all elements in it"""

    import numpy as np
    outList = []
    for i in nestedList:
        if type(i) in (list, np.ndarray):
            outList.extend(unlist(i))
        else:
            outList.append(i)
    return outList


def getValues(l, index):
    """from list l get the values at indices specified by index"""
    return [l[i] for i in index]


def getColors(colorPalette, numOfCol):
    return [colorPalette(i/numOfCol) for i in range(numOfCol)]


def subList(lst1, lst2):
    return list(map(op.sub, lst1, lst2))


def intersect(*lists):
    import numpy as np
    return reduce(np.intersect1d, list(lists))


def fileRemove(fName):
    try:
        remove(fName)
    except OSError as e:
        if e.errno != 2:    # 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
            raise


def mergeRanges(ranges):
    """
    Take a 2D numpy array, with each row as a range and return merged ranges i.e. ranges which are overlapping would be
    combined as well.
    :param ranges:
    :return:
    """
    from collections import deque
    import numpy as np
    if len(ranges) < 2:
        return ranges
    ranges = ranges[ranges[:, 0].argsort()]
    for i in ranges:
        if i[0] > i[1]:
            garb = i[0]
            i[0] = i[1]
            i[1] = garb
    min_value = ranges[0, 0]
    max_value = ranges[0, 1]
    out_range = deque()
    for i in ranges[1:]:
        if i[0] > max_value:
            out_range.append([min_value, max_value])
            min_value = i[0]
            max_value = i[1]
        elif i[1] > max_value:
            max_value = i[1]
    out_range.append([min_value, max_value])
    return np.array(out_range)


def revcomp(seq):
    assert type(seq) == str
    old = 'ACGTRYKMBDHVacgtrykmbdhv'
    rev = 'TGCAYRMKVHDBtgcayrmkvhdb'
    tab = str.maketrans(old, rev)
    return seq.translate(tab)[::-1]


def readfasta(f):
    from gzip import open as gzopen
    from gzip import BadGzipFile
    from collections import deque
    import sys
    out = {}
    chrid = ''
    chrseq = deque()
    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzopen(f, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                            chrseq = deque()
                        else:
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip().decode())
        else:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    if chrid != '':
        out[chrid] = ''.join(chrseq)
    # TODO: add check for the validation of input fasta files
    return out
# END

def cgtpl(cg):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    return [i.split(';') for i in cg.split(',')[:-1]]
#end


def plotref(refidx, varpos, syrireg, figout, bw=100000):
    """
    Plots distribution of bed file on a genome annotated with SRs
    :param refidx: path to ref.faidx
    :param syrireg: path to syri regions in BEDPE format
    :param varpos: variant positions in BED
    :param figout: output figure path
    :return:
    """
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    import pybedtools as bt # TODO: remove this dependency
    from collections import deque, defaultdict
    import numpy as np
    import logging
    logger = logging.getLogger('plotref')
    def _readbed(vp, refbed):
        _chrs = set([r[0] for r in refbed])
        bincnt = defaultdict(deque)
        skipchrs = []
        curchr = ''
        pos = deque()
        added_chrs = list()
        with open(vp, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if len(line) < 3:
                    logger.warning("Incomplete information in BED file at line: {}. Skipping it.".format("\t".join(line)))
                    continue
                if line[0] not in _chrs:
                    if line[0] not in skipchrs:
                        logger.info("Chromosome in BED is not present in FASTA or not selected for plotting. Skipping it. BED line: {}".format("\t".join(line)))
                        skipchrs.append(line[0])
                    continue
                if curchr == '':
                    curchr = line[0]
                    pos.append([int(line[1]), int(line[2])])
                elif curchr == line[0]:
                    pos.append([int(line[1]), int(line[2])])
                else:
                    if line[0] in added_chrs:
                        logger.error("BED file: {} is not sorted. For plotting tracks, sorted bed file is required. Exiting.".format(vp))
                        sys.exit()
                    if len(pos) > 1:
                        rngs = mergeRanges(np.array(pos))
                    else:
                        rngs = pos
                    chrpos = np.array(list(set([i for r in rngs for i in range(r[0], r[1])])))
                    # Get bin breakpoints for the chromosome
                    bins = np.concatenate((np.arange(0, [r[2] for r in refbed if r[0] == curchr][0], bw), np.array([r[2] for r in refbed if r[0] == curchr])), axis=0)
                    binval = np.histogram(chrpos, bins)[0]
                    bincnt[curchr] = deque([((bins[i] + bins[i+1])/2, binval[i]/bw) for i in range(len(binval))])
                    added_chrs.append(curchr)
                    # Set the new chromosome
                    curchr = line[0]
                    pos = deque([[int(line[1]), int(line[2])]])
            if len(pos) > 1:
                rngs = mergeRanges(np.array(pos))
            else:
                rngs = pos
            chrpos = np.array(list(set([i for r in rngs for i in range(r[0], r[1])])))
            # Get bin breakpoints for the chromosome
            bins = np.concatenate((np.arange(0, [r[2] for r in refbed if r[0] == curchr][0], bw), np.array([r[2] for r in refbed if r[0] == curchr])), axis=0)
            binval = np.histogram(chrpos, bins)[0]
            bincnt[curchr] = deque([((bins[i] + bins[i+1])/2, binval[i]/bw) for i in range(len(binval))])
        return bincnt
    # END
    chr_height = 0.3
    th = 0.6    # Track height for plotting bedfile data
    refbed = readfaidxbed(refidx)
    # TODO: remove dependency on pybedtools
    # syrioutbed = bt.BedTool('/srv/netscratch/dep_mercier/grp_schneeberger/projects/read_vs_assembly/results/human/hg002/reads15kb/svcalls/syri/hg002.reads15kb.winnowmap.hap1syri.bedpe').sort()
    syrioutbed = bt.BedTool(syrireg).sort()
    bedbin = _readbed(varpos, refbed)
    fig = plt.figure(figsize=[12, 10])
    ax = fig.add_subplot()
    ax.set_ylim([0, len(refbed)+1])
    ax.set_xlim([0, max([r[2] for r in refbed])])
    colors = {'SYN': 'lightgrey', 'INV': '#FFA500', 'TRANS': '#9ACD32', 'DUP': '#00BBFF'}
    peakpos = defaultdict()
    for i in range(len(refbed)):
        patches = deque()
        r = refbed[i]
        y = len(refbed) - i
        ax.add_patch(Rectangle((r[1], y), r[2], chr_height, fill=False, linewidth=0.5))
        # ax.add_patch(Rectangle((r[1], y), r[2], chr_height, linewidth=0.5, color='black'))
        bed = syrioutbed.filter(lambda b: b.chrom == r[0]).saveas()
        for b in bed:
            patches.append(Rectangle((b.start, y), b.stop-b.start, chr_height, color=colors[b[6]], linewidth=0))
        ax.add_collection(PatchCollection(patches, match_original=True))
        chrpos = [k[0] for k in bedbin[r[0]]]
        tpos = [k[1] for k in bedbin[r[0]]]
        tposmax = max(tpos)
        y0 = len(refbed) - i + chr_height
        ypos = [(t*th/tposmax)+y0 for t in tpos]
        ax.fill_between(chrpos, ypos, y0, color='blue', lw=0.1, zorder=2)
        y_mid = y0 + (th/2)
        peakpos[refbed[i][0]] = [(chrpos[j], tpos[j]) for j in range(len(ypos)) if ypos[j] > y_mid]
    plt.tight_layout()
    fig.savefig(figout)
    plt.close()
    return peakpos
# END