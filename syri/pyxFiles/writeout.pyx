import numpy as np
from syri.bin.func.myUsefulFunctions import *
import sys
import time
from igraph import *
from collections import Counter, deque, defaultdict
from scipy.stats import *
from datetime import datetime, date
import pandas as pd
import os
from Bio.SeqIO import parse
import logging

np.random.seed(1)


##################################################################
### Generate formatted output
##################################################################

def getsrtable(cwdpath, prefix):
    import re

    files = ["synOut.txt", "invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt", "ctxOut.txt"]
    entries = defaultdict()

    re_dup = re.compile("dup", re.I)
    align_num = 1
    sr_num = 0
    for file in files:
        sr_type = file.split("O")[0]
        if sr_type == "syn":
            sr_type = "SYN"
        elif sr_type == "inv":
            sr_type = "INV"
        elif sr_type == "TL":
            sr_type = "TRANS"
        elif sr_type == "invTL":
            sr_type = "INVTR"
        elif sr_type == "dup":
            sr_type = "DUP"
        elif sr_type == "invDup":
            sr_type = "INVDP"

        with open(cwdpath + prefix + file, "r") as fin:
            for line in fin:
                line = line.strip().split("\t")
                if line[0] == "#":
                    sr_num += 1
                    if file == 'ctxOut.txt':
                        if line[8] == "translocation":
                            sr_type = "TRANS"
                        elif line[8] == "invTranslocation":
                            sr_type = "INVTR"
                        elif line[8] == "duplication":
                            sr_type = "DUP"
                        elif line[8] == "invDuplication":
                            sr_type = "INVDP"
                    entries[sr_type + str(sr_num)] = {
                        'achr': line[1],
                        'astart': line[2],
                        'aend': line[3],
                        'bchr': line[5],
                        'bstart': line[6],
                        'bend': line[7],
                        'vartype': sr_type,
                        'parent': "-",
                        'dupclass': "-",
                        'aseq': "-",
                        "bseq": "-"
                    }
                    if re_dup.search(file):
                        entries[sr_type + str(sr_num)]["dupclass"] = "copygain" if line[9] == "B" else "copyloss"
                    if file == "ctxOut.txt":
                        entries[sr_type + str(sr_num)]["vartype"] = sr_type
                        if re_dup.search(line[8]):
                            entries[sr_type + str(sr_num)]["dupclass"] = "copygain" if line[9] == "B" else "copyloss"
                else:
                    entries[sr_type + "AL" + str(align_num)] = {
                        'achr': entries[sr_type + str(sr_num)]['achr'],
                        'astart': line[0],
                        'aend': line[1],
                        'bchr': entries[sr_type + str(sr_num)]['bchr'],
                        'bstart': line[2],
                        'bend': line[3],
                        'vartype': sr_type + "AL",
                        'parent': sr_type + str(sr_num),
                        'dupclass': "-",
                        'aseq': "-",
                        "bseq": "-"
                    }
                    align_num += 1

    anno = pd.DataFrame.from_dict(entries, orient="index")
    anno.loc[:, ['astart', 'aend', 'bstart', 'bend']] = anno.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
    anno['id'] = anno.index.values
    anno = anno.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
    anno.sort_values(['achr', 'astart', 'aend'], inplace=True)
    return anno


def extractseq(_gen, _pos):
    chrs = defaultdict(dict)
    for fasta in parse(_gen, 'fasta'):
        if fasta.id in _pos.keys():
            chrs[fasta.id] = {_i:fasta.seq[_i-1] for _i in _pos[fasta.id]}
    return chrs


def getTSV(cwdpath, prefix, ref):
    """
    :param cwdpath: Path containing all input files
    :return: A TSV file containing genomic annotation for the entire genome
    """
    import pandas as pd
    import sys
    from collections import defaultdict
    logger = logging.getLogger("getTSV")
    logger.debug('cwdpath:' + cwdpath + ", prefix:" + prefix + ", ref:" + ref)

    logger.debug('Get SR anno')
    anno = getsrtable(cwdpath, prefix)
    logger.debug("Number of SR annotations: " + str(anno.shape[0]))

    logger.debug('Get SV data')
    hasSV = True
    if not os.path.isfile(cwdpath + prefix + "sv.txt"):
        hasSV = False
        logger.info(cwdpath + prefix + "sv.txt"+' cannot be opened. Cannot output SVs.')

    if hasSV:
        svdata = pd.read_table(cwdpath + prefix + "sv.txt", header=None)
        logger.debug("Number of SV annotations read from file: " + str(svdata.shape[0]))
        svdata.columns = ["vartype", "astart", 'aend', 'bstart', 'bend', 'achr', 'bchr']

        entries = deque()
        count = 1
        for row in svdata.itertuples(index=False):
            if row.vartype == "#":
                _parent = anno.loc[(anno.achr == row.achr) &
                                   (anno.astart == row.astart) &
                                   (anno.aend == row.aend) &
                                   (anno.bchr == row.bchr) &
                                   (anno.bstart == row.bstart) &
                                   (anno.bend == row.bend) &
                                   (anno.parent == "-"), "id"]
                if len(_parent) != 1:
                    logger.error("Error in finding parent for SV")
                    logger.error(row.to_string() + "\t" + _parent.to_string())
                    sys.exit()
                else:
                    _parent = _parent.to_string(index=False, header=False)
                continue
            entries.append({
                'achr': row.achr,
                'astart': row.astart,
                'aend': row.aend,
                'bchr': row.bchr,
                'bstart': row.bstart,
                'bend': row.bend,
                'vartype': row.vartype,
                'parent': _parent,
                'id': row.vartype + str(count),
                'dupclass': "-",
                'aseq': "-",
                "bseq": "-"
            })
            count += 1

        sv = pd.DataFrame.from_records(entries)
        if sv.shape[0] != 0:
            logger.debug("NO SV found in SV file. NO SV will be reported." + str(sv.shape[0]))
            sv.index = sv['id']
            sv.loc[:, ['astart', 'aend', 'bstart', 'bend']] = sv.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
            sv = sv.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
            sv.sort_values(['achr', 'astart', 'aend'], inplace=True)
        else:
            sv = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
    else:
        sv = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
    logger.debug("Number of SV annotations to output: " + str(sv.shape[0]))

    logger.debug('Get notal data')
    hasNotal = True
    if not os.path.isfile(cwdpath + prefix + "notAligned.txt"):
        hasNotal = False
        logger.info(cwdpath + prefix + "notAligned.txt"+' cannot be opened. Cannot output not aligned regions.')

    if hasNotal:
        isempty = False
        try:
            notal = pd.read_table(cwdpath + prefix + "notAligned.txt", header=None)
        except pd.errors.EmptyDataError as e:
            isempty = True
            logger.debug("NOTAL file is empty")
            notal = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass', 'selected'])

        if not isempty:
            logger.debug("Number of NOTAL annotations read from file: " + str(notal.shape[0]))

            notal.columns = ["gen", "start", "end", "chr"]
            notal[["start", "end"]] = notal[["start", "end"]].astype("int")
            entries = defaultdict()
            _c = 0
            for row in notal.itertuples(index=False):
                _c += 1
                if row.gen == "R":
                    entries["NOTAL" + str(_c)] = {
                        'achr': row.chr,
                        'astart': row.start,
                        'aend': row.end,
                        'bchr': "-",
                        'bstart': "-",
                        'bend': "-",
                        'vartype': "NOTAL",
                        'parent': "-",
                        'id': "NOTAL" + str(_c),
                        'dupclass': "-",
                        'aseq': "-",
                        "bseq": "-"
                    }
                elif row.gen == "Q":
                    entries["NOTAL" + str(_c)] = {
                        'achr': "-",
                        'astart': "-",
                        'aend': "-",
                        'bchr': row.chr,
                        'bstart': row.start,
                        'bend': row.end,
                        'vartype': "NOTAL",
                        'parent': "-",
                        'id': "NOTAL" + str(_c),
                        'dupclass': "-",
                        'aseq': "-",
                        'bseq': "-"
                    }
            notal = pd.DataFrame.from_dict(entries, orient="index")
            notal = notal.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
            notal.sort_values(['achr', 'astart', 'aend'], inplace=True)
            notal['selected'] = -1
    else:
        notal = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass', 'selected'])
    logger.debug("Number of NOTAL annotations to output: " + str(notal.shape[0]))


    def p_indel():
        vtype = "INS" if indel == 1 else "DEL"
        entries.append({
            'achr': _ac,
            'astart': _as,
            'aend': _ae,
            'bchr': _bc,
            'bstart': _bs,
            'bend': _be,
            'vartype': vtype,
            'parent': _p,
            'aseq': "-" if vtype == "INS" else _seq,
            'bseq': "-" if vtype == "DEL" else _seq,
            'id': vtype + str(count),
            'dupclass': "-"
        })

    entries = deque()

    logger.debug('Get SNP data')
    hasSNP = True
    if not os.path.isfile(cwdpath + prefix + "snps.txt"):
        hasSNP = False
        logger.info(cwdpath + prefix + "snps.txt"+' cannot be opened. Cannot output SNPs and short indels.')

    if hasSNP:
        with open(cwdpath + prefix + "snps.txt", "r") as fin:
            indel = 0                           # marker for indel status. 0 = no_indel, 1 = insertion, -1 = deletion
            _as = -1                            # astart
            _ae = -1
            _bs = -1
            _be = -1
            _ac = -1
            _bc = -1
            _p = -1
            _seq = ""

            i = 1
            for line in fin:
                i += 1
                line = line.strip().split("\t")
                try:
                    if line[0] == "#" and len(line) == 7:
                        if indel != 0:
                            p_indel()
                            indel = 0
                            _seq = ""
                        _parent = anno.loc[(anno.achr == line[5]) &
                                           (anno.astart == int(line[1])) &
                                           (anno.aend == int(line[2])) &
                                           (anno.bchr == line[6]) &
                                           (anno.bstart == int(line[3])) &
                                           (anno.bend == int(line[4])) &
                                           (anno.parent == "-"), "id"]
                        if len(_parent) != 1:
                            logger.error("Error in finding parent for SNP")
                            logger.error("\t".join(line) + "\t" + _parent.to_string())
                            sys.exit()
                        else:
                            _parent = _parent.to_string(index=False, header=False)
                        continue
                    elif line[1] != "." and line[2] != "." and len(line) == 12:
                        if indel != 0:
                            p_indel()
                            indel = 0
                            _seq = ""
                        count += 1
                        entries.append({
                            'achr': line[10],
                            'astart': int(line[0]),
                            'aend': int(line[0]),
                            'bchr': line[11],
                            'bstart': int(line[3]),
                            'bend': int(line[3]),
                            'vartype': "SNP",
                            'parent': _parent,
                            'aseq': line[1],
                            'bseq': line[2],
                            'id': "SNP" + str(count),
                            'dupclass': "-"
                        })
                    elif indel == 0:
                        count += 1
                        _as = int(line[0])
                        _ae = int(line[0])
                        _bs = int(line[3])
                        _be = int(line[3])
                        _ac = line[10]
                        _bc = line[11]
                        _p = _parent
                        indel = 1 if line[1] == "." else -1
                        _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                    elif indel == 1:
                        if int(line[0]) != _as or line[1] != "." or line[10] != _ac or int(line[3])!=(_be+1):
                            p_indel()
                            _seq = ""
                            count += 1
                            _as = int(line[0])
                            _ae = int(line[0])
                            _bs = int(line[3])
                            _be = int(line[3])
                            _ac = line[10]
                            _bc = line[11]
                            _p = _parent
                            indel = 1 if line[1] == "." else -1
                            _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                        else:
                            _be = int(line[3])
                            _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                    elif indel == -1:
                        if int(line[3]) != _bs or line[2] != "." or line[11] != _bc or int(line[0])!=(_ae+1):
                            p_indel()
                            _seq = ""
                            count += 1
                            _as = int(line[0])
                            _ae = int(line[0])
                            _bs = int(line[3])
                            _be = int(line[3])
                            _ac = line[10]
                            _bc = line[11]
                            _p = _parent
                            indel = 1 if line[1] == "." else -1
                            _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                        else:
                            _ae = int(line[0])
                            _seq = _seq + line[2] if line[1] == "." else _seq + line[1]
                except IndexError as _e:
                    logger.error(_e)
                    logger.error("\t".join(line))
                    sys.exit()

        snpdata = pd.DataFrame.from_records(entries)
        try:
            snpdata.loc[:, ['astart', 'aend', 'bstart', 'bend']] = snpdata.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
            snpdata.index = snpdata['id']
            snpdata = snpdata.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
            snpdata.sort_values(['achr', 'astart', 'aend'], inplace=True)
        except KeyError as e:
            snpdata = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
            logger.warning('No SNPs were found. This could be an error. Remove any old snps.txt and try re-running without --nosnp.')
    else:
        snpdata = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
    logger.debug("Number of SNPs annotations read: " + str(snpdata.shape[0]))

    positions = defaultdict()
    for _chr in snpdata.achr.unique():
        positions[_chr] = snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "INS"), "astart"].tolist() + (snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "DEL"), "astart"] - 1).tolist()
        # positions[chrom] = (snpdata.loc[(snpdata.achr == chrom) & (snpdata.vartype == "DEL"), "astart"] - 1).tolist()

    seq = extractseq(ref, positions)
    logger.debug('Fixing coords for insertions and deletions')

    for _chr in snpdata.achr.unique():
        ## Fix coordinates and sequence for insertions
        _indices = snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "INS")].index.values
        _seq = pd.Series([seq[_chr][_i] for _i in snpdata.loc[_indices, "astart"]], index=_indices)
        _dir = _indices[~snpdata.loc[_indices, "parent"].str.contains("INV")]
        _inv = _indices[snpdata.loc[_indices, "parent"].str.contains("INV")]

        # Add previous characted
        snpdata.loc[_indices, "aseq"] = _seq
        snpdata.loc[_dir, "bseq"] = _seq.loc[_dir] + snpdata.loc[_dir, "bseq"]
        snpdata.loc[_inv, "bseq"] = _seq.loc[_inv] + snpdata.loc[_inv, "bseq"].str[::-1]

        # Change coordinates
        snpdata.loc[_dir, "bstart"] = snpdata.loc[_dir, "bstart"] - 1
        snpdata.loc[_inv, "bend"] = snpdata.loc[_inv, "bstart"] + snpdata.loc[_inv, "bend"]
        snpdata.loc[_inv, "bstart"] = snpdata.loc[_inv, "bend"] - snpdata.loc[_inv, "bstart"] + 1
        snpdata.loc[_inv, "bend"] = snpdata.loc[_inv, "bend"] - snpdata.loc[_inv, "bstart"] + 1

        ## Fix coordinates and sequence for deletions
        _indices = snpdata.loc[(snpdata.achr == _chr) & (snpdata.vartype == "DEL")].index.values
        _seq = pd.Series([seq[_chr][_i-1] for _i in snpdata.loc[_indices, "astart"]], index=_indices)
        _dir = _indices[~snpdata.loc[_indices, "parent"].str.contains("INV")]
        _inv = _indices[snpdata.loc[_indices, "parent"].str.contains("INV")]

        # Add previous characted
        snpdata.loc[_indices, "aseq"] = _seq + snpdata.loc[_indices, "aseq"]
        snpdata.loc[_indices, "bseq"] = _seq

        # Change coordinates
        snpdata.loc[_indices, "astart"] = snpdata.loc[_indices, "astart"] - 1
        snpdata.loc[_inv, "bstart"] = snpdata.loc[_inv, "bstart"] + 1
        snpdata.loc[_inv, "bend"] = snpdata.loc[_inv, "bend"] + 1

    events = anno.loc[anno.parent == "-"]
    logger.debug('Starting output file generation')

    with open(cwdpath + prefix + "syri.out", "w") as fout:
        logger.debug("All annotation count. " + "SR anno: " + str(anno.shape[0]) + " SV anno: " + str(sv.shape[0]) + " ShV anno: " + str(snpdata.shape[0]) + " notal anno: " + str(notal.shape[0]))
        annogrp = anno.groupby('parent')
        svgrp = sv.groupby('parent')
        snpdatagrp = snpdata.groupby('parent')

        notA = notal.loc[notal.achr != "-"].copy()
        notA.loc[:, ["astart", "aend"]] = notA.loc[:, ["astart", "aend"]].astype("int")
        notB = notal.loc[notal.bchr != "-"].copy()
        notB.loc[:, ["bstart", "bend"]] = notB.loc[:, ["bstart", "bend"]].astype("int")
        row_old = -1
        for row in events.itertuples(index=False):
            if len(notA) > 0:
                if row_old != -1 and row_old.achr != row.achr:
                    _notA = notA.loc[(notA.achr == row_old.achr) & (notA.astart == row_old.aend+1) & (notA.selected != 1), notA.columns != 'selected']
                    notA.loc[(notA.achr == row_old.achr) & (notA.astart == row_old.aend + 1), 'selected'] = 1
                    if len(_notA) == 0:
                        pass
                    elif len(_notA) == 1:
                        fout.write("\t".join(list(map(str, _notA.iloc[0]))) + "\n")
                    else:
                        logger.error("too many notA regions")
                        sys.exit()
                _notA = notA.loc[(notA.achr == row.achr) & (notA.aend == row.astart-1) & (notA.selected != 1), notA.columns != 'selected']
                notA.loc[(notA.achr == row.achr) & (notA.aend == row.astart - 1), 'selected'] = 1
                if len(_notA) == 0:
                    pass
                elif len(_notA) == 1:
                    fout.write("\t".join(list(map(str, _notA.iloc[0]))) + "\n")
                else:
                    logger.error("too many notA regions")
                    sys.exit()
            fout.write("\t".join(list(map(str, row))) + "\n")
            row_old = row

            try:
                a = annogrp.get_group(row.id)
            except KeyError:
                a = pd.DataFrame()
            except Exception as e:
                logger.debug('Error in finding key for anno.' + e)

            try:
                b = svgrp.get_group(row.id)
            except KeyError as k:
                b = pd.DataFrame()
            except Exception as e:
                logger.debug('Error in finding key for anno.' + e)

            try:
                c = snpdatagrp.get_group(row.id)
            except KeyError:
                c = pd.DataFrame()
            except Exception as e:
                logger.debug('Error in finding key for anno.' + e)

            outdata = pd.concat([a, b, c])
            outdata.sort_values(["astart", "aend"], inplace=True)
            fout.write(outdata.to_csv(sep="\t", index=False, header=False))
        if len(notA) > 0:
            _notA = notA.loc[(notA.achr == row_old.achr) & (notA.astart == row_old.aend+1) & (notA.selected != 1), notA.columns != 'selected']
            notA.loc[(notA.achr == row_old.achr) & (notA.astart == row_old.aend + 1), 'selected'] = 1
            if len(_notA) == 0:
                pass
            elif len(_notA) == 1:
                fout.write("\t".join(list(map(str, _notA.iloc[0]))) + "\n")
            else:
                logger.error("too many notA regions")
                sys.exit()
        fout.write(notB.loc[:, notB.columns != 'selected'].to_csv(sep="\t", index=False, header=False))

    logger.debug('Remapping query genome ids')
    if os.path.isfile(cwdpath+prefix+"mapids.txt"):
        chroms = {}
        with open(cwdpath+prefix+"mapids.txt", "r") as m:
            for line in m:
                l = line.strip().split()
                chroms[l[0]] = l[1]
        lines = open(cwdpath + prefix + "syri.out", "r").readlines()
        with open(cwdpath + prefix + "syri.out", "w") as fout:                 ## CHANGE THIS BEFORE PUSHING
            for line in lines:
                line = line.strip().split('\t')
                if line[5] != "-":
                    line[5] = chroms[line[5]]
                fout.write("\t".join(line) + '\n')
    return 0


def getVCF(finname, foutname, cwdpath, prefix):
    """
    does not output notal in qry genome
    :param finname:
    foutname:
    :return:
    """
    data = pd.read_table(cwdpath + prefix + finname, header=None)
    data.columns = ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']
    data = data.loc[data['achr'] != "-"].copy()
    data.sort_values(['achr', 'astart', 'aend'], inplace=True)
    data.loc[:, ['astart', 'aend', 'bstart', 'bend']] = data.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype(str)

    with open(cwdpath + prefix + foutname, 'w') as fout:
        fout.write('##fileformat=VCFv4.3\n')
        fout.write('##fileDate=' + str(date.today()).replace('-', '') + '\n')
        fout.write('##source=syri\n')
        fout.write('##ALT=<ID=SYN,Description="Syntenic region">' + '\n')
        fout.write('##ALT=<ID=INV,Description="Inversion">' + '\n')
        fout.write('##ALT=<ID=TRANS,Description="Translocation">' + '\n')
        fout.write('##ALT=<ID=INVTR,Description="Inverted Translocation">' + '\n')
        fout.write('##ALT=<ID=DUP,Description="Duplication">' + '\n')
        fout.write('##ALT=<ID=INVDP,Description="Inverted Duplication">' + '\n')
        fout.write('##ALT=<ID=SYNAL,Description="Syntenic alignment">' + '\n')
        fout.write('##ALT=<ID=INVAL,Description="Inversion alignment">' + '\n')
        fout.write('##ALT=<ID=TRANSAL,Description="Translocation alignment">' + '\n')
        fout.write('##ALT=<ID=INVTRAL,Description="Inverted Translocation alignment">' + '\n')
        fout.write('##ALT=<ID=DUPAL,Description="Duplication alignment">' + '\n')
        fout.write('##ALT=<ID=INVDPAL,Description="Inverted Duplication alignment">' + '\n')
        fout.write('##ALT=<ID=HDR,Description="Highly diverged regions">' + '\n')
        fout.write('##ALT=<ID=INS,Description="Insertion in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=DEL,Description="Deletion in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=CPG,Description="Copy gain in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=CPL,Description="Copy loss in non-reference genome">' + '\n')
        fout.write('##ALT=<ID=SNP,Description="Single nucleotide polymorphism">' + '\n')
        fout.write('##ALT=<ID=TDM,Description="Tandem repeat">' + '\n')
        fout.write('##ALT=<ID=NOTAL,Description="Not Aligned region">' + '\n')
        fout.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position on reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=ChrB,Number=1,Type=String,Description="Chromoosme ID on the non-reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=StartB,Number=1,Type=Integer,Description="Start position on non-reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=EndB,Number=1,Type=Integer,Description="End position on non-reference genome">' + '\n')
        fout.write('##INFO=<ID=Parent,Number=1,Type=String,Description="ID of the parent SR">' + '\n')
        fout.write(
            '##INFO=<ID=VarType,Number=1,Type=String,Description="Start position on non-reference genome">' + '\n')
        fout.write(
            '##INFO=<ID=DupType,Number=1,Type=String,Description="Copy gain or loss in the non-reference genome">' + '\n')
        # fout.write('##INFO=<ID=NotAlGen,Number=1,Type=String,Description="Genome containing the not aligned region">' + '\n')

        fout.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']) + '\n')
        for line in data.itertuples(index=False):
            pos = [line[0], line[1], ".", 'N', '<' + line[10] + '>', '.', 'PASS']

            if line[10] in ["SYN", "INV", "TRANS", "INVTR", "DUP", "INVDUP"]:
                _info = ';'.join(['END='+line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent=.', 'VarType='+'SR', 'DupType='+line[11]])
                pos.append(_info)
                fout.write('\t'.join(pos) + '\n')

            elif line[10] == "NOTAL":
                if line[0] != "-":
                    _info = ';'.join(
                        ['END=' + line[2], 'ChrB=.', 'StartB=.', 'EndB=.', 'Parent=.', 'VarType=.', 'DupType=.'])
                    pos.append(_info)
                    fout.write('\t'.join(pos) + '\n')

            elif line[10] in ["SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDUPAL"]:
                _info = ';'.join(
                    ['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=.', 'DupType=.'])
                pos.append(_info)
                fout.write('\t'.join(pos) + '\n')

            elif line[10] in ['CPG', 'CPL', 'TDM', 'HDR']:
                _info = ";".join(['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=ShV', 'DupType=.'])
                pos.append(_info)
                fout.write('\t'.join(pos) + '\n')

            elif line[10] in ['SNP', 'DEL', 'INS']:
                if (line[3] == '-') != (line[4] == '-'):
                    sys.exit("Inconsistency in annotation type. Either need seq for both or for none.")
                elif line[3] == '-' and line[4] == '-':
                    _info = ";".join(['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=ShV', 'DupType=.'])
                    pos.append(_info)
                    fout.write('\t'.join(pos) + '\n')
                elif line[3] != '-' and line[4] != '-':
                    pos = [line[0], line[1], ".", line[3], line[4], '.', 'PASS']
                    _info = ";".join(['END=' + line[2], 'ChrB=' + line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent=' + line[9], 'VarType=ShV', 'DupType=.'])
                    pos.append(_info)
                    fout.write('\t'.join(pos) + '\n')
    return 0
