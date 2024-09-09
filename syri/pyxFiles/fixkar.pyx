# cython: language_level = 3
# distutils: language = c++

from collections import defaultdict, namedtuple
from itertools import product
import sys
from warnings import warn
from syri import __version__
import pandas as pd
import numpy as np
from collections import defaultdict, deque, Counter
from syri.scripts.func import mergeRanges, readfasta, revcomp, readfasta_iterator, isgzip, load_lrs
from syri.synsearchFunctions import readSAMBAM, readPAF
import operator



def getdata(reflength, refid, refdir):
    lrs = load_lrs()
    refdata = {}
    # For each scaffold in genome A find the list of aligning scaffolds in genome B
    for rid in refid:
        qiddata = {}
        if len(reflength[rid]) == 0:
            refdata[rid] = qiddata
            continue
        # Use sliding window (size: 5) for ref scaffolds larger than 50 units,
        # scaffolds which are not major alignment in any window are removed.
        # print("Sliding Window")
        if len(reflength[rid]) > 50:
            wsize = 5
            scaflist = []
            for i in range(0, len(reflength[rid]) - wsize + 1):
                qids = deque()
                for j in range(i, i+5):
                    try:
                        qids.append(reflength[rid][j*10000])
                    except KeyError:
                        pass
                qids = Counter(qids)
                for k, v in qids.items():
                    if v >= 3:
                        scaflist.append(k)
            scaflist = list(np.unique(scaflist))
            if len(scaflist) == 0:
                print('Error: For contig ', rid, ' no matching contig was found. This could be a result of incorrect assembly or extensive repeats or ', rid, ' could be a novel region. Pseudo-genome cannot be generated. Try removing this contig.')
                sys.exit()
        else:
            scaflist = list(np.unique(list(reflength[rid].values())))

        qids = np.unique(list(reflength[rid].values()))
        qind = dict(zip(qids, list(range(len(qids)))))
        ks = deque()
        qlist = deque()
        for k in sorted(list(reflength[rid].keys())):
            if reflength[rid][k] in scaflist:
                ks.append(k)
                qlist.append(qind[reflength[rid][k]])
        ks = list(ks)
        qlist = list(qlist)
        try:
            selected = lrs(qlist)
        except Exception as e:
            raise Exception('Error in running lrs: {}'.format(e))
            sys.exit()

        q = reflength[rid][ks[selected[0]]]
        start = ks[selected[0]]
        end = ks[selected[0]]
        for s in selected[1:]:
            if reflength[rid][ks[s]] == q:
                end = ks[s]
            else:
                qiddata[q] = {'s': start,
                              'e': end,
                              'l': end - start + 10000
                              }
                q = reflength[rid][ks[s]]
                start = ks[s]
                end = ks[s]
        qiddata[q] = {'s': start,
                      'e': end,
                      'l': end - start + 10000
                      }

        for qid in qiddata.keys():
            qiddata[qid]['d'] = refdir[rid][qid]
            refdata[rid] = qiddata
    return refdata
# END


def remove_multialign(refdata, qrydata, fixref):
    for rid in refdata.keys():
        if len(refdata[rid].keys()) == 0:
            continue

        if fixref:
            sqid = min(refdata[rid].keys(), key=lambda x: refdata[rid][x]['s'])
            eqid = max(refdata[rid].keys(), key=lambda x: refdata[rid][x]['e'])
            qids = list(refdata[rid].keys())
            qids.remove(sqid)
            try:
                qids.remove(eqid)
            except ValueError:
                pass
        else:
            qids = list(refdata[rid].keys())

        for qid in qids:
            if len(qrydata[qid]) > 1:
                if qrydata[qid][rid]['l'] > sum([qrydata[qid][i]['l'] for i in qrydata[qid].keys() if i != rid]):
                    ids = list(qrydata[qid].keys())
                    for i in ids:
                        if i == rid:
                            continue
                        try:
                            qrydata[qid].pop(i)
                        except KeyError:
                            pass
                        try:
                            refdata[i].pop(qid)
                        except KeyError:
                            pass
                else:
                    try:
                        refdata[rid].pop(qid)
                    except KeyError:
                        pass
                    try:
                        qrydata[qid].pop(rid)
                    except KeyError:
                        pass

    for qid in qrydata.keys():
        if not fixref:
            continue

        if len(qrydata[qid].keys()) == 0:
            continue

        srid = min(qrydata[qid].keys(), key=lambda x: qrydata[qid][x]['s'])
        erid = max(qrydata[qid].keys(), key=lambda x: qrydata[qid][x]['e'])
        rids = list(qrydata[qid].keys())
        rids.remove(srid)
        try:
            rids.remove(erid)
        except ValueError:
            pass

        for rid in rids:
            if len(refdata[rid]) > 1:
                if refdata[rid][qid]['l'] > sum([refdata[rid][i]['l'] for i in refdata[rid].keys() if i != qid]):
                    ids = list(refdata[rid].keys())
                    for i in ids:
                        if i == qid:
                            continue
                        try:
                            refdata[rid].pop(i)
                        except KeyError:
                            pass
                        try:
                            qrydata[i].pop(rid)
                        except KeyError:
                            pass
                else:
                    try:
                        refdata[rid].pop(qid)
                    except KeyError:
                        pass
                    try:
                        qrydata[qid].pop(rid)
                    except KeyError:
                        pass

    empty = [i for i in refdata.keys() if len(refdata[i]) == 0]
    for i in empty:
        refdata.pop(i)

    empty = [i for i in qrydata.keys() if len(qrydata[i]) == 0]
    for i in empty:
        qrydata.pop(i)
    return refdata, qrydata
# END


def scaf(incoords, refsize, qrysize, args):
    from gzip import open as gzopen
    from collections import OrderedDict

    refin   = args.ref.name
    qryfin   = args.qry.name
    # nn      = "N"*args.ncount
    # out     = args.out
    fixref   = args.fixref
    F       = args.ftype
    # Read coords and genome size

    # if F == 'T':    coords = pd.read_table(args.coords.name, header=None)
    # elif F == 'B':  coords = readSAMBAM(args.coords.name, type='B')
    # elif F == 'S':  coords = readSAMBAM(args.coords.name, type='S')
    # elif F == 'P':  coords = readPAF(args.coords.name)
    # coords = coords[list(range(11))]
    # coords.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr"]
    coords = incoords.copy()
    # coords.aChr = "ref"+coords.aChr
    coords.aChr =  "ref"+coords.aChr
    coords.bChr =  "qry"+coords.bChr
    # coords.bChr = "qry"+coords.bChr
    # refsize = {("ref"+id): len(seq) for id, seq in readfasta(refin).items()}
    # qrysize = {("qry"+id): len(seq) for id, seq in readfasta(qryin).items()}
    # fin = gzopen(refin, 'r') if isgzip(refin) else open(refin, 'r')
    # refsize = {f'ref{c}': len(s) for c, s in readfasta_iterator(refin)}
    # fin = gzopen(qryfin, 'r') if isgzip(qryfin) else open(qryfin, 'r')
    # qrysize = {f'qry{c}': len(s) for c, s in readfasta_iterator(qryfin)}
    refsize = {f"ref{k}": v for k, v in refsize.items()}
    qrysize = {f"qry{k}": v for k, v in qrysize.items()}

    reflength = defaultdict(dict)
    for i in np.unique(coords.aChr):
        for r in range(0, refsize[i], 10000):
            a = coords.loc[(coords.aChr == i) & (coords.aStart < (r + 10000)) & (coords.aEnd > r)]
            if a.shape[0] == 0:
                continue
            dd = dict()
            for j in np.unique(a.bChr):
                b = a.loc[a.bChr == j]
                if b.shape[0] == 0:
                    continue
                temp = mergeRanges(np.array(b[["aStart", "aEnd"]]))
                temp[temp[:, 0] < r, 0] = r
                temp[temp[:, 1] > (r+10000), 1] = r+10000
                dd[j] = (temp[:, 1] - temp[:, 0]).sum()
            reflength[i][r] = max(dd.items(), key=operator.itemgetter(1))[0]

    refdir = defaultdict(dict)
    for rid in reflength.keys():
        for qid in np.unique(list(reflength[rid].values())):
            c = coords.loc[(coords.aChr == rid) & (coords.bChr == qid)]
            cdir = c.loc[c.bDir == 1]
            cinv = c.loc[c.bDir == -1]
            cdm = mergeRanges(np.array(cdir[["aStart", "aEnd"]]))
            cim = mergeRanges(np.array(cinv[["aStart", "aEnd"]]))
            refdir[rid][qid] = 1 if (cdm[:, 1] - cdm[:, 0]).sum() >= (cim[:, 1] - cim[:, 0]).sum() else -1

    qrylength = defaultdict(dict)
    loci = coords.bDir == -1

    # Invert coordinates for inverted alignments
    coords.loc[loci, "bStart"] = coords.loc[loci, "bStart"] + coords.loc[loci, "bEnd"]
    coords.loc[loci, "bEnd"] = coords.loc[loci, "bStart"] - coords.loc[loci, "bEnd"]
    coords.loc[loci, "bStart"] = coords.loc[loci, "bStart"] - coords.loc[loci, "bEnd"]

    for i in np.unique(coords.bChr):
        for r in range(0, qrysize[i], 10000):
            a = coords.loc[(coords.bChr == i) & (coords.bStart < (r + 10000)) & (coords.bEnd > r)]
            if a.shape[0] == 0:
                continue
            dd = dict()
            for j in np.unique(a.aChr):
                b = a.loc[a.aChr == j]
                if b.shape[0] == 0:
                    continue
                temp = mergeRanges(np.array(b[["bStart", "bEnd"]]))
                temp[temp[:, 0] < r, 0] = r
                temp[temp[:, 1] > (r+10000), 1] = r+10000
                dd[j] = (temp[:, 1] - temp[:, 0]).sum()
            qrylength[i][r] = max(dd.items(), key=operator.itemgetter(1))[0]

    qrydir = defaultdict(dict)
    for qid in qrylength.keys():
        for rid in np.unique(list(qrylength[qid].values())):
            c = coords.loc[(coords.aChr == rid) & (coords.bChr == qid)]
            cdir = c.loc[c.bDir == 1]
            cinv = c.loc[c.bDir == -1]
            cdm = mergeRanges(np.array(cdir[["bStart", "bEnd"]]))
            cim = mergeRanges(np.array(cinv[["bStart", "bEnd"]]))
            qrydir[qid][rid] = 1 if (cdm[:, 1] - cdm[:, 0]).sum() >= (cim[:, 1] - cim[:, 0]).sum() else -1

    # Revert coordinates back to original
    coords.loc[loci, "bStart"] = coords.loc[loci, "bStart"] + coords.loc[loci, "bEnd"]
    coords.loc[loci, "bEnd"] = coords.loc[loci, "bStart"] - coords.loc[loci, "bEnd"]
    coords.loc[loci, "bStart"] = coords.loc[loci, "bStart"] - coords.loc[loci, "bEnd"]

    refid = np.unique(coords.aChr)
    qryid = np.unique(coords.bChr)

    # Get mapping data between scaffolds
    refdata = getdata(reflength, refid, refdir)
    qrydata = getdata(qrylength, qryid, qrydir)

    # remove scaffolds (Y) from the list of matched scaffolds for scaffold (X) if
    # X is not in the list of matched scaffolds for Y.
    for rid in refdata.keys():
        qids = list(refdata[rid].keys())
        for qid in qids:
            if rid not in qrydata[qid].keys():
                refdata[rid].pop(qid)

    for qid in qrydata.keys():
        rids = list(qrydata[qid].keys())
        for rid in rids:
            if qid not in refdata[rid].keys():
                qrydata[qid].pop(rid)

    # Scaffolds aligning in the middle cannot align with other scaffolds
    # Check if the middle scaffold is aligning with multiple scaffolds
    # If yes, then decide whether to remove it or not based on length of alignments
    refdata, qrydata = remove_multialign(refdata, qrydata, fixref)

    if not fixref:
        for qid in qrydata.keys():
            if len(qrydata[qid].keys()) > 1:
                print("branching found for query scaffold", qid)

    # Get the order in which to access the chromosomes
    chrs = deque()
    rpass = deque()
    qpass = deque()
    for rid in refdata:
        if rid not in rpass:
            rpass.append(rid)
            rgroup = deque([rid])
            qgroup = deque()
            rtogo = deque()
            qtogo = deque(list(refdata[rid].keys()))
            while len(rtogo) != 0 or len(qtogo) != 0:
                for r in rtogo:
                    if r not in rpass:
                        rpass.append(r)
                        rgroup.append(r)
                        qtogo.extend(list(refdata[r].keys()))
                rtogo = deque()

                for q in qtogo:
                    if q not in qpass:
                        qpass.append(q)
                        qgroup.append(q)
                        rtogo.extend(list(qrydata[q].keys()))
                qtogo = deque()
            chrs.append((rgroup, qgroup))

    chrid = 1
    # fref = open(out+".ref.fasta", "w")
    # fqry = open(out+'.qry.fasta', 'w')
    # fout = open(out + ".anno", 'w')
    # size = 60
    chrom_orders = deque()
    for chr in chrs:
        loci = {}
        count = 0

        # <editor-fold desc="Initiate nodes for reference and query genome">
        # Initiate nodes for qry genome
        rqids = {}
        for rid in chr[0]:
            pos = 0
            for qid in refdata[rid].keys():
                loci[count] = [rid, pos]
                count += 1
                pos += 1
            loci[count] = [rid, pos]
            count += 1

            # Ordered qids
            rqids[rid] = sorted(refdata[rid].keys(), key=lambda x: refdata[rid][x]['s'])

        # Initiate nodes for qry genome
        qrids = {}
        for qid in chr[1]:
            pos = 0
            for rid in qrydata[qid].keys():
                loci[count] = [qid, pos]
                count += 1
                pos += 1
            loci[count] = [qid, pos]
            count += 1

            # Ordered rids
            qrids[qid] = sorted(qrydata[qid].keys(), key=lambda x: qrydata[qid][x]['s'])
        # </editor-fold>

        # <editor-fold desc="Add neighbouring nodes">
        # Add upper self-neighbour
        for k, v in loci.items():
            if k-1 in loci.keys():
                if v[0] == loci[k-1][0]:
                    if v[1]-1 == loci[k-1][1]:
                        loci[k].append(k-1)
                    else:
                        loci[k].append("-")
                else:
                    loci[k].append("-")
            else:
                loci[k].append("-")

        # Add lower self-neighbour
        for k, v in loci.items():
            if k+1 in loci.keys():
                if v[0] == loci[k+1][0]:
                    if v[1]+1 == loci[k+1][1]:
                        loci[k].append(k+1)
                    else:
                        loci[k].append("-")
                else:
                    loci[k].append("-")
            else:
                loci[k].append("-")
        # </editor-fold>

        rhead = {}
        qhead = {}
        # get heads for rids/qids
        for k, v in loci.items():
            if v[2] == "-":
                if "ref" in v[0]:
                    rhead[v[0]] = k
                elif "qry" in v[0]:
                    qhead[v[0]] = k

        # Add empty list in which alignment neighbours will be added
        for v in loci.items():
            v[1].append([])

        # <editor-fold desc="Add alignment neighbours">
        for rid in chr[0]:
            rh = rhead[rid]
            for qid in rqids[rid]:
                # get rid index in qrid
                qh = qhead[qid] + qrids[qid].index(rid)
                if refdata[rid][qid]['d'] == 1:
                    loci[rh][4].append(qh)
                    loci[rh+1][4].append(qh+1)
                elif refdata[rid][qid]['d'] == -1:
                    loci[rh][4].append(qh+1)
                    loci[rh+1][4].append(qh)
                rh = loci[rh][3]

        for qid in chr[1]:
            qh = qhead[qid]
            for rid in qrids[qid]:
                rh = rhead[rid] + rqids[rid].index(qid)
                if qrydata[qid][rid]['d'] == 1:
                    loci[qh][4].append(rh)
                    loci[qh+1][4].append(rh+1)
                elif qrydata[qid][rid]['d'] == -1:
                    loci[qh][4].append(rh+1)
                    loci[qh+1][4].append(rh)
                qh = loci[qh][3]
        # </editor-fold>

        endlist = []
        for k, v in loci.items():
            if v[2] == "-" or v[3] == "-":
                if loci[v[4][0]][2] == "-" or loci[v[4][0]][3] == '-':
                    if "ref" in v[0]:
                        endlist.append([k, v[4][0]])
                    elif "qry" in v[0]:
                        endlist.append([v[4][0], k])

        if len(endlist) == 0:
            print("Following scaffolds form circular chromosome. Filtered out.")
            print(chr)
            continue
        else:
            endlist = np.unique(endlist, axis=0)

        # <editor-fold desc="Find all linear paths between the connected nodes">
        paths = []
        for e in endlist:
            ends = e.copy()
            rdir = 1 if loci[ends[0]][2] == "-" else -1
            qdir = 1 if loci[ends[1]][2] == "-" else -1

            rout = [ends[0]]
            qout = [ends[1]]

            routstack = []
            qoutstack = []
            rdirstack = []
            endtogo = []
            while True:
                if rdir != 0:
                    ends[0] = loci[ends[0]][3] if rdir == 1 else loci[ends[0]][2]
                if qdir != 0:
                    ends[1] = loci[ends[1]][3] if qdir == 1 else loci[ends[1]][2]

                if rdir == 0:
                    ends[0] = [i for i in [loci[ends[0]][2], loci[ends[0]][3]] if i in loci[ends[1]][4]][0]
                if qdir == 0:
                    ends[1] = [i for i in [loci[ends[1]][2], loci[ends[1]][3]] if i in loci[ends[0]][4]][0]

                rout.append(ends[0])
                qout.append(ends[1])

                # Check whether the nodes/edges match
                if ends[1] not in loci[ends[0]][4] or ends[0] not in loci[ends[1]][4]:
                    print("error in traversing")
                    print(chr)
                    break

                # Check if there is branching. i.e. more than one possible directions more traversing
                # add current out-path to stack
                if len(loci[ends[0]][4]) > 1 and len(loci[ends[1]][4]) > 1:
                    routstack.append(rout.copy())
                    qoutstack.append(qout.copy())
                    endtogo.append([ends[0], [i for i in loci[ends[0]][4] if i != ends[1]][0]])
                    rdirstack.append(rdir)
                    ends[0] = [i for i in loci[ends[1]][4] if i != ends[0]][0]
                    rout.append(ends[0])
                    if len(loci[ends[0]][4]) == 1:
                        rdir = 1 if loci[ends[0]][2] == "-" else -1
                    else:
                        rdir = 0
                    continue

                # When end node is reached: go back to previous branch or end loop
                if len(loci[ends[0]][4]) == 1 and len(loci[ends[1]][4]) == 1:
                    paths.append([rout, qout])
                    if len(endtogo) > 0:
                        ends = endtogo.pop()
                        rout = routstack.pop()
                        qout = qoutstack.pop()
                        qout.append(ends[1])
                        rdir = rdirstack.pop()
                        if len(loci[ends[1]][4]) == 1:
                            qdir = 1 if loci[ends[1]][2] == "-" else -1
                        else:
                            qdir = 0
                        continue
                    else:
                        break

                # Switch alignment if end of current alignment is obtained
                if loci[ends[0]][2] == "-" or loci[ends[0]][3] == "-":
                    ends[0] = [i for i in loci[ends[1]][4] if i != ends[0]][0]
                    rout.append(ends[0])
                    if len(loci[ends[0]][4]) == 1:
                        rdir = 1 if loci[ends[0]][2] == "-" else -1
                    else:
                        rdir = 0
                    continue

                if loci[ends[1]][2] == "-" or loci[ends[1]][3] == "-":
                    ends[1] = [i for i in loci[ends[0]][4] if i != ends[1]][0]
                    qout.append(ends[1])
                    if len(loci[ends[1]][4]) == 1:
                        qdir = 1 if loci[ends[1]][2] == "-" else -1
                    else:
                        qdir = 0
                    continue
        # </editor-fold>

        # Can't recall why the below check was added.
        if any([len(Counter(Counter(i[0]).values())) > 1 or len(Counter(Counter(i[1]).values())) > 1 for i in paths]):
            continue

        if len(paths) == 0:
            print('Could not assemble ', chr[0][0])
            continue

        # Remove inverse paths
        unipaths = []
        for path in paths:
            if (path[0] not in [i[0][::-1] for i in unipaths]) and (path[1] not in [i[1][::-1] for i in unipaths]):
                unipaths.append(path.copy())

        bestpath = defaultdict()
        for path in unipaths:
            rids = pd.unique([loci[i][0] for i in path[0]])
            qids = pd.unique([loci[i][0] for i in path[1]])
            score = np.mean([np.sum([refsize[i] for i in rids]), np.sum([qrysize[i] for i in qids])])

            if len(bestpath) == 0:
                bestpath = {'rpath': path[0],
                            'qpath': path[1],
                            'rids': rids,
                            'qids': qids,
                            'score': score}
            elif score > bestpath['score']:
                bestpath = {'rpath': path[0],
                            'qpath': path[1],
                            'rids': rids,
                            'qids': qids,
                            'score': score}

        # <editor-fold desc="Get directions of the selected chromosomes">
        rdirs = {}
        last = ""
        for i in bestpath['rpath']:
            if loci[i][0] not in rdirs:
                if len(loci[i][4]) == 2:
                    last = loci[i][0]
                else:
                    if loci[i][2] == '-':
                        if last != '':
                            rdirs[loci[i][0]] = -1
                            last = ''
                        else:
                            rdirs[loci[i][0]] = 1
                            last = ''
                    elif loci[i][3] == '-':
                        if last != '':
                            rdirs[loci[i][0]] = 1
                            last = ''
                        else:
                            rdirs[loci[i][0]] = -1
                            last = ''

        qdirs = {}
        last = ''
        for i in bestpath['qpath']:
            if loci[i][0] not in qdirs:
                if len(loci[i][4]) == 2:
                    last = loci[i][0]
                else:
                    if loci[i][2] == '-':
                        if last != '':
                            qdirs[loci[i][0]] = -1
                            last = ''
                        else:
                            qdirs[loci[i][0]] = 1
                            last = ''
                    elif loci[i][3] == '-':
                        if last != '':
                            qdirs[loci[i][0]] = 1
                            last = ''
                        else:
                            qdirs[loci[i][0]] = -1
                            last = ''
        # </editor-fold>
        chrom_orders.append((bestpath, rdirs, qdirs))

        # if fixref:
        #     fref.write(">Pseudochrom{}\n".format(chrid))
        #     spacer = ""
        #     for rid in bestpath['rids']:
        #         fref.write(spacer)
        #         spacer = nn
        #         if rdirs[rid] == 1:
        #             fref.write([str(seq) for id, seq in readfasta(refin).items() if ("ref"+id) == rid][0])
        #         elif rdirs[rid] == -1:
        #             fref.write([revcomp(seq) for id, seq in readfasta(refin).items() if ("ref"+id) == rid][0])
        #     fref.write("\n")
        #
        # if fixref:
        #     fqry.write(">Pseudochrom{}\n".format(chrid))
        # else:
        #     fqry.write(">"+chr[0][0].replace('ref', '')+'\n')
        # spacer = ""
        # for qid in bestpath['qids']:
        #     fqry.write(spacer)
        #     spacer = nn
        #     if qdirs[qid] == 1:
        #         fqry.write([seq for id, seq in readfasta(qryin).items() if ("qry"+id) == qid][0])
        #     elif qdirs[qid] == -1:
        #         fqry.write([revcomp(seq) for id, seq in readfasta(qryin).items() if ("qry"+id) == qid][0])
        # fqry.write("\n")
        # if fixref:
        #     fout.write(">Pseudochrom{}\n".format(chrid))
        # else:
        #     fout.write(">"+chr[0][0].replace('ref', '')+'\n')
        # fout.write("\t".join(list(bestpath['rids'])).replace('ref', '') + '\n')
        # fout.write("\t".join([str(rdirs[i]) for i in bestpath['rids']]) + '\n')
        # fout.write("\t".join(list(bestpath['qids'])).replace('qry', '') + '\n')
        # fout.write("\t".join([str(qdirs[i]) for i in bestpath['qids']]) + '\n')
        chrid += 1
    # fref.close()
    # fqry.close()
    # fout.close()

    cleanout = [(c[0]['rids'], c[0]['qids'], c[1], c[2]) for c in chrom_orders]
    for i in range(len(cleanout)):
        c0 = OrderedDict({c.replace("ref", ""): cleanout[i][2][c] for c in cleanout[i][0]})
        c1 = OrderedDict({c.replace("qry", ""): cleanout[i][3][c] for c in cleanout[i][1]})
        # c1 = [c.replace("ref", "") for c in cleanout[i][0]]
        # c1 = [c.replace("qry", "") for c in cleanout[i][1]]
        # c2 = {k.replace("ref", ""): v for k, v in cleanout[i][2].items()}
        # c3 = {k.replace("qry", ""): v for k, v in cleanout[i][3].items()}
        cleanout[i] = [c0, c1]
    return cleanout
# END


# # if __name__ == "__main__":
# def main():
#     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     parser.add_argument("coords", help='Alignment coordinates in a tsv format', type=argparse.FileType("r"))
#     parser.add_argument("ref", help='Assembly of genome A in multi-fasta format', type=argparse.FileType("r"))
#     parser.add_argument("qry", help='Assembly of genome B in multi-fasta format', type=argparse.FileType("r"))
# #    parser.add_argument("-l", dest="length", help='length of the range', type=int, default=10000)
#     parser.add_argument("-n", dest='ncount', help="number of N's to be inserted", type=int, default=500)
#     parser.add_argument('-o', dest='out', help="output file prefix", default="out", type=str)
#     parser.add_argument('-noref', dest='noref', help="Use this parameter when no assembly is at chromosome level", default=False, action='store_true')
#     parser.add_argument('-F', dest="ftype", help="Input coords type. T: Table, S: SAM, B: BAM, P:PAF", default="T", choices=['T', 'S', 'B', 'P'])
#     parser.add_argument('--version', action='version', version='{version}'.format(version=__version__))
#     args = parser.parse_args()
#     scaf(args)
#     print('Finished chroder')
