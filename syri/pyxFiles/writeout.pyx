# cython: language_level = 3
# distutils: language = c++

import numpy as np
from syri.scripts.func import *
from collections import deque, defaultdict
from datetime import date
import pandas as pd
import os
import logging
import sys

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
# END


def extractseq(gen: str, pos: defaultdict):
    chrs = defaultdict(dict)
    for chrid, seq in readfasta(gen).items():
        if chrid in pos.keys():
            chrs[chrid] = {i: seq[i-1] for i in pos[chrid]}
    return chrs
# END


def parsesvs(f: str, anno: pd.DataFrame, count: int, ref: str):
    """
    Read the SVs from sv.txt. Update the position and base to include upstream base/pos so that the SVs are suitable for writing in VCF format as well
    :param fin: path to sv.txt
    :param anno: SR annotation coordinates
    :param count: variation counter
    :param ref: variation counter
    :return: Dataframe containing all SVs with corrected reference pos (for indels) and sequence (after adding upstream base with correct orientiation)
    """
    # Define logger
    logger = logging.getLogger("parsesvs")
    logger.debug('f: {}'.format(f))
    if not os.path.isfile(f):
        logger.info('{f} cannot be opened. Cannot output SVs.'.format(f=f))
        sv = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
    else:
        # Read sv.txt
        svdata = pd.read_table(f, header=None)
        logger.debug("Number of SV annotations read from file: " + str(svdata.shape[0]))
        svdata.columns = ["vartype", "astart", 'aend', 'bstart', 'bend', 'achr', 'bchr', 'aseq', 'bseq']
        # Ensure that chromosome names are strings
        svdata[['achr', 'bchr']] = svdata[['achr', 'bchr']].astype('str')
        # Read data
        entries = deque()
        for row in svdata.itertuples(index=False):
            if row.vartype == "#":
                parent = anno.loc[(anno.achr == row.achr) &
                                   (anno.astart == row.astart) &
                                   (anno.aend == row.aend) &
                                   (anno.bchr == row.bchr) &
                                   (anno.bstart == row.bstart) &
                                   (anno.bend == row.bend) &
                                   (anno.parent == "-"), "id"]
                if len(parent) != 1:
                    logger.error("Error in finding parent for SV")
                    logger.error(row.to_string() + "\t" + parent.to_string())
                    sys.exit()
                else:
                    parent = parent.to_string(index=False, header=False)
                continue
            entries.append({
                'achr': row.achr,
                'astart': row.astart,
                'aend': row.aend,
                'bchr': row.bchr,
                'bstart': row.bstart,
                'bend': row.bend,
                'vartype': row.vartype,
                'parent': parent,
                'id': row.vartype + str(count),
                'dupclass': "-",
                'aseq': row.aseq,
                "bseq": row.bseq
            })
            count += 1
        sv = pd.DataFrame.from_records(entries)
        if sv.shape[0] != 0:
            logger.debug("SV found in SV file. Number of SV that will be reported." + str(sv.shape[0]))
            sv.index = sv['id']
            sv.loc[:, ['astart', 'aend', 'bstart', 'bend']] = sv.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
            sv = sv.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
            sv.sort_values(['achr', 'astart', 'aend'], inplace=True)
            # Get sequence at start position in reference sequence
            positions = defaultdict()
            for grp in sv.groupby('achr'):
                p = deque()
                p.extend(grp[1].loc[grp[1].vartype == "INS", "astart"].tolist())
                p.extend((grp[1].loc[grp[1].vartype == "DEL", "astart"] - 1).tolist())
                p.extend((grp[1].loc[grp[1].vartype == "HDR", "astart"] - 1).tolist())
                positions[grp[0]] = sorted(p)
            seq = extractseq(ref, positions)
            logger.debug('Fixing coords for insertions and deletions')
            # Update sequence and position for insertions
            indices = sv.loc[sv.vartype == "INS"].index.values
            s = pd.Series([seq[row[0]][row[1]] for row in sv.loc[indices].itertuples(index=False)], index=indices, dtype=str)
            sv.loc[indices, "aseq"] = s
            sv.loc[indices, "bseq"] = s + sv.loc[indices, 'bseq']
            d = indices[~sv.loc[indices, "parent"].str.contains("INV")]
            sv.loc[d, "bstart"] = sv.loc[d, "bstart"] - 1
            # For inverted annotations, select the position matching the selected reference base and swap the coordinates so that  in output qry start < end
            i = indices[sv.loc[indices, "parent"].str.contains("INV")]
            sv.loc[i, "bstart"] = sv.loc[i, "bstart"] + 1
            sv.loc[i, "bstart"] = sv.loc[i, "bstart"] + sv.loc[i, "bend"]
            sv.loc[i, "bend"] = sv.loc[i, "bstart"] - sv.loc[i, "bend"]
            sv.loc[i, "bstart"] = sv.loc[i, "bstart"] - sv.loc[i, "bend"]

            # Update sequence and positions for deletions
            indices = sv.loc[sv.vartype == "DEL"].index.values
            s = pd.Series([seq[row[0]][row[1]-1] for row in sv.loc[indices].itertuples(index=False)], index=indices, dtype=str)
            sv.loc[indices, "aseq"] = s + sv.loc[indices, 'aseq']
            sv.loc[indices, "bseq"] = s
            sv.loc[indices, "astart"] = sv.loc[indices, "astart"] - 1
            ## For deletions in inverted annotations, do not need to change coordinates

            # Update sequence and positions for HDRs
            indices = sv.loc[sv.vartype == "HDR"].index.values
            s = pd.Series([seq[row[0]][row[1]-1] for row in sv.loc[indices].itertuples(index=False)], index=indices)
            sv.loc[indices, "aseq"] = s + sv.loc[indices, 'aseq']
            sv.loc[indices, "bseq"] = s + sv.loc[indices, 'bseq']
            sv.loc[indices, "astart"] = sv.loc[indices, "astart"] - 1
            d = indices[~sv.loc[indices, "parent"].str.contains("INV")]
            sv.loc[d, "bstart"] = sv.loc[d, "bstart"] - 1
            # For inverted annotations, select the position matching the selected reference base and swap the coordinates so that  in output qry start < end
            i = indices[sv.loc[indices, "parent"].str.contains("INV")]
            sv.loc[i, "bstart"] = sv.loc[i, "bstart"] + 1
            sv.loc[i, "bstart"] = sv.loc[i, "bstart"] + sv.loc[i, "bend"]
            sv.loc[i, "bend"] = sv.loc[i, "bstart"] - sv.loc[i, "bend"]
            sv.loc[i, "bstart"] = sv.loc[i, "bstart"] - sv.loc[i, "bend"]
        else:
            logger.debug("NO SV found in SV file. NO SV will be reported.")
            sv = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
    logger.debug("Number of SV annotations to output: " + str(sv.shape[0]))
    return sv, count
# END

def parsesnps(f: str, anno: pd.DataFrame, count: int, ref: str):
    def p_indel():
        vtype = "INS" if indel == 1 else "DEL"
        entries.append({
            'achr': achr,
            'astart': astart,
            'aend': aend,
            'bchr': bchr,
            'bstart': bstart,
            'bend': bend,
            'vartype': vtype,
            'parent': p,
            'aseq': "-" if vtype == "INS" else seq,
            'bseq': "-" if vtype == "DEL" else seq,
            'id': vtype + str(count),
            'dupclass': "-"
        })
    # END
    logger = logging.getLogger("parsesnps")
    logger.debug('parsing SNPs and short indels for writing')
    logger.debug('f: {}'.format(f))
    if not os.path.isfile(f):
        logger.warning('{f} cannot be opened. Cannot output SNPs and short indels.'.format(f=f))
        snp = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
    else:
        entries = deque()
        with open(f, "r") as fin:
            indel = 0                           # marker for indel status. 0 = no_indel, 1 = insertion, -1 = deletion
            astart, aend, bstart, bend, achr, bchr, p = [-1] * 7
            seq = ""
            i = 1
            for line in fin:
                i += 1
                line = line.strip().split("\t")
                try:
                    if line[0] == "#" and len(line) == 7:
                        if indel != 0:
                            p_indel()
                            indel = 0
                            seq = ""
                        parent = anno.loc[(anno.achr == line[5]) &
                                          (anno.astart == int(line[1])) &
                                          (anno.aend == int(line[2])) &
                                          (anno.bchr == line[6]) &
                                          (anno.bstart == int(line[3])) &
                                          (anno.bend == int(line[4])) &
                                          (anno.parent == "-"), "id"]
                        if len(parent) != 1:
                            logger.error("Error in finding parent for SNP")
                            logger.error("\t".join(line) + "\t" + parent.to_string())
                            sys.exit()
                        else:
                            parent = parent.to_string(index=False, header=False)
                        continue
                    elif line[1] != "." and line[2] != "." and len(line) == 12:
                        if indel != 0:
                            p_indel()
                            indel = 0
                            seq = ""
                        count += 1
                        entries.append({
                            'achr': line[10],
                            'astart': int(line[0]),
                            'aend': int(line[0]),
                            'bchr': line[11],
                            'bstart': int(line[3]),
                            'bend': int(line[3]),
                            'vartype': "SNP",
                            'parent': parent,
                            'aseq': line[1],
                            'bseq': line[2],
                            'id': "SNP" + str(count),
                            'dupclass': "-"
                        })
                    elif indel == 0:
                        count += 1
                        astart, aend = int(line[0]), int(line[0])
                        bstart, bend = int(line[3]), int(line[3])
                        achr, bchr = line[10], line[11]
                        p = parent
                        indel = 1 if line[1] == "." else -1
                        seq = seq + line[2] if line[1] == "." else seq + line[1]
                    elif indel == 1:
                        n = bend-1 if 'INV' in parent else bend+1
                        if int(line[0]) != astart or line[1] != "." or line[10] != achr or int(line[3]) != n:
                            p_indel()
                            seq = ""
                            count += 1
                            astart, aend = int(line[0]), int(line[0])
                            bstart, bend = int(line[3]), int(line[3])
                            achr, bchr = line[10], line[11]
                            p = parent
                            indel = 1 if line[1] == "." else -1
                            seq = seq + line[2] if line[1] == "." else seq + line[1]
                        else:
                            bend = int(line[3])
                            seq = seq + line[2] if line[1] == "." else seq + line[1]
                    elif indel == -1:
                        if int(line[3]) != bstart or line[2] != "." or line[11] != bchr or int(line[0]) != (aend+1):
                            p_indel()
                            seq = ""
                            count += 1
                            astart, aend = int(line[0]), int(line[0])
                            bstart, bend = int(line[3]), int(line[3])
                            achr, bchr = line[10], line[11]
                            p = parent
                            indel = 1 if line[1] == "." else -1
                            seq = seq + line[2] if line[1] == "." else seq + line[1]
                        else:
                            aend = int(line[0])
                            seq = seq + line[2] if line[1] == "." else seq + line[1]
                except IndexError as _e:
                    logger.error(_e)
                    logger.error("\t".join(line))
                    sys.exit()

        snp = pd.DataFrame.from_records(entries)
        try:
            snp.loc[:, ['astart', 'aend', 'bstart', 'bend']] = snp.loc[:, ['astart', 'aend', 'bstart', 'bend']].astype('int')
            snp.index = snp['id']
            snp = snp.loc[:, ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
            snp.sort_values(['achr', 'astart', 'aend'], inplace=True)
        except KeyError as e:
            snp = pd.DataFrame(columns=['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass'])
            logger.warning('No SNPs were found. This could be an error. Remove any old snps.txt and try re-running without --nosnp.')
        logger.debug("Number of SNPs annotations read: " + str(snp.shape[0]))
        positions = defaultdict()
        for c in snp.achr.unique():
            positions[c] = snp.loc[(snp.achr == c) & (snp.vartype == "INS"), "astart"].tolist() + (snp.loc[(snp.achr == c) & (snp.vartype == "DEL"), "astart"] - 1).tolist()
        seq = extractseq(ref, positions)
        logger.debug('Fixing coords for insertions and deletions')

        # Update sequence and positions for deletions
        indices = snp.loc[snp.vartype == "INS"].index.values
        s = pd.Series([seq[row[0]][row[1]] for row in snp.loc[indices].itertuples(index=False)], index=indices)
        snp.loc[indices, "aseq"] = s
        snp.loc[indices, "bseq"] = s + snp.loc[indices, 'bseq']
        d = indices[~snp.loc[indices, "parent"].str.contains("INV")]
        snp.loc[d, "bstart"] = snp.loc[d, "bstart"] - 1
        # For inverted annotations, select the position matching the selected reference base and swap the coordinates so that  in output qry start < end
        i = indices[snp.loc[indices, "parent"].str.contains("INV")]
        snp.loc[i, "bstart"] = snp.loc[i, "bstart"] + 1
        snp.loc[i, "bstart"] = snp.loc[i, "bstart"] + snp.loc[i, "bend"]
        snp.loc[i, "bend"] = snp.loc[i, "bstart"] - snp.loc[i, "bend"]
        snp.loc[i, "bstart"] = snp.loc[i, "bstart"] - snp.loc[i, "bend"]

        # Update sequence and positions for deletions
        indices = snp.loc[snp.vartype == "DEL"].index.values
        s = pd.Series([seq[row[0]][row[1]-1] for row in snp.loc[indices].itertuples(index=False)], index=indices)
        snp.loc[indices, "aseq"] = s + snp.loc[indices, 'aseq']
        snp.loc[indices, "bseq"] = s
        snp.loc[indices, "astart"] = snp.loc[indices, "astart"] - 1
        ## For deletions in inverted annotations, do not need to change coordinates
    return snp, count
# END


def getTSV(cwdpath: str, prefix: str, ref: str, hdrseq: bool, maxs: int):
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
    count = 1
    # Read structure variants
    logger.debug('Get SV data')
    sv, count = parsesvs(cwdpath + prefix + "sv.txt", anno, count, ref)
    if not hdrseq:
        sv.loc[sv.vartype == 'HDR', 'aseq'] = '-'
        sv.loc[sv.vartype == 'HDR', 'bseq'] = '-'
    if maxs != -1:
        sv.loc[(sv.vartype.isin(['HDR', 'DEL'])) & (sv.aend - sv.astart > maxs), 'aseq'] = '-'
        sv.loc[(sv.vartype.isin(['HDR', 'DEL'])) & (sv.aend - sv.astart > maxs), 'bseq'] = '-'
        sv.loc[(sv.vartype.isin(['HDR', 'INS'])) & (sv.bend - sv.bstart > maxs), 'aseq'] = '-'
        sv.loc[(sv.vartype.isin(['HDR', 'INS'])) & (sv.bend - sv.bstart > maxs), 'bseq'] = '-'
    # Read not aligned regions
    logger.debug('Get notal data')
    hasNotal = True
    if not os.path.isfile(cwdpath + prefix + "notAligned.txt"):
        hasNotal = False
        logger.info(cwdpath + prefix + "notAligned.txt"+' cannot be opened. Cannot output not aligned regions.')
    if hasNotal:
        isempty = False
        try:
            notal = pd.read_table(cwdpath + prefix + "notAligned.txt", header=None)
            notal[3] = notal[3].astype('str')
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
    # Read short variants
    logger.debug('Get SNPs/indels data')
    snp, count = parsesnps(cwdpath + prefix + "snps.txt", anno, count, ref)
    if maxs != -1:
        snp.loc[(snp.vartype == 'DEL') & (snp.aend - snp.astart > maxs), 'aseq'] = '-'
        snp.loc[(snp.vartype == 'DEL') & (snp.aend - snp.astart > maxs), 'bseq'] = '-'
        snp.loc[(snp.vartype == 'INS') & (snp.bend - snp.bstart > maxs), 'aseq'] = '-'
        snp.loc[(snp.vartype == 'INS') & (snp.bend - snp.bstart > maxs), 'bseq'] = '-'

    events = anno.loc[anno.parent == "-"]
    logger.debug('Starting output file generation')
    with open(cwdpath + prefix + "syri.out", "w") as fout:
        logger.debug("All annotation count. " + "SR anno: " + str(anno.shape[0]) + " SV anno: " + str(sv.shape[0]) + " ShV anno: " + str(snp.shape[0]) + " notal anno: " + str(notal.shape[0]))
        annogrp = anno.groupby('parent')
        svgrp = sv.groupby('parent')
        snpgrp = snp.groupby('parent')

        notA = notal.loc[notal.achr != "-"].copy()
        # notA.loc[:, ["astart", "aend"]] = notA.loc[:, ["astart", "aend"]].astype("int")
        notA[["astart", "aend"]] = notA[["astart", "aend"]].astype("int")
        notB = notal.loc[notal.bchr != "-"].copy()
        notB[["bstart", "bend"]] = notB[["bstart", "bend"]].astype("int")
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
            # Update row_old when chromosome change or the max coordinate increases
            if row_old == -1: row_old = row
            elif row.achr != row_old.achr: row_old = row
            elif row_old.aend < row.aend: row_old = row
            # row_old = row

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
                c = snpgrp.get_group(row.id)
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
        with open(cwdpath + prefix + "syri.out", "w") as fout:
            for line in lines:
                line = line.strip().split('\t')
                if line[5] != "-":
                    line[5] = chroms[line[5]]
                fout.write("\t".join(line) + '\n')
    return 0
# END


def getVCF(finname, foutname, cwdpath, prefix, sname):
    """
    does not output notal in qry genome
    :param finname:
    foutname:
    :return:
    """
    logger = logging.getLogger("getVCF")
    data = pd.read_table(cwdpath + prefix + finname, header=None, keep_default_na=False, dtype=object)
    data.columns = ['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']
    data = data.loc[data['achr'] != "-"].copy()
    dtypes = {'achr': str,
              'astart': int,
              'aend': int,
              'aseq': str,
              'bseq': str,
              'bchr': str,
              'bstart': str,
              'bend': str,
              'id': str,
              'parent': str,
              'vartype': str,
              'dupclass': str}
    data = data.astype(dtypes)
    # Translation string for fixing IUPAC codes in sequence strings
    old = 'ACGTNacgtnRYSWKMBDHVryswkmbdhv'
    rev = 'ACGTNacgtnACCAGACAAAaccagacaaa'
    tab = str.maketrans(old, rev)
    try:
        data['achr'] = data['achr'].astype('int')
    except ValueError as ve:
        logger.debug('Chromosome values are sorted lexicographically.')
    data.sort_values(['achr', 'astart', 'aend'], inplace=True)
    data.loc[:, ['achr', 'astart', 'aend', 'bstart', 'bend']] = data.loc[:, ['achr', 'astart', 'aend', 'bstart', 'bend']].astype(str)
    chr_sizes = data.aend.astype(int).groupby(data.achr).max()
    with open(cwdpath + prefix + foutname, 'w') as fout:
        fout.write('##fileformat=VCFv4.3\n')
        fout.write('##fileDate=' + str(date.today()).replace('-', '') + '\n')
        fout.write('##source=syri\n')
        for cntg, ln in chr_sizes.items():
            fout.write('##contig=<ID=' + cntg + ',length=' + str(ln) + '>\n')
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
            '##INFO=<ID=VarType,Number=1,Type=String,Description="SR for syntenic regions, ShV for short variants, missing otherwise">' + '\n')
        fout.write(
            '##INFO=<ID=DupType,Number=1,Type=String,Description="Copy gain or loss in the non-reference genome">' + '\n')
        fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n')
        # fout.write('##INFO=<ID=NotAlGen,Number=1,Type=String,Description="Genome containing the not aligned region">' + '\n')

        fout.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sname]) + '\n')

        ## Format information
        frmt = '\t'.join(['GT', '1'])
        ## Iterate over each line from syri.out
        for line in data.itertuples(index=False):
            pos = [line[0], line[1], line[8], 'N', '<' + line[10] + '>', '.', 'PASS']

            if line[10] in ["SYN", "INV", "TRANS", "INVTR", "DUP", "INVDP"]:
                _info = ';'.join(['END='+line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent=.', 'VarType='+'SR', 'DupType='+line[11]])

            elif line[10] == "NOTAL":
                if line[0] != "-":
                    _info = ';'.join(['END=' + line[2], 'ChrB=.', 'StartB=.', 'EndB=.', 'Parent=.', 'VarType=.', 'DupType=.'])

            elif line[10] in ["SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDPAL"]:
                _info = ';'.join(['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=.', 'DupType=.'])

            elif line[10] in ['CPG', 'CPL', 'TDM']:
                _info = ";".join(['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=ShV', 'DupType=.'])

            elif line[10] in ['SNP', 'DEL', 'INS', 'HDR']:
                if (line[3] == '-') != (line[4] == '-'):
                    logger.error("Inconsistency in annotation type. Either need seq for both or for none.")
                elif line[3] == '-' and line[4] == '-':
                    _info = ";".join(['END=' + line[2], 'ChrB='+line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent='+line[9], 'VarType=ShV', 'DupType=.'])
                elif line[3] != '-' and line[4] != '-':
                    if line[3].translate(tab).upper() == line[4].translate(tab).upper(): continue
                    pos = [line[0], line[1], line[8], line[3].translate(tab), line[4].translate(tab), '.', 'PASS']
                    _info = ";".join(['END=' + line[2], 'ChrB=' + line[5], 'StartB='+line[6], 'EndB='+line[7], 'Parent=' + line[9], 'VarType=ShV', 'DupType=.'])

            pos.append(_info)
            pos.append(frmt)
            fout.write('\t'.join(pos) + '\n')
    return 0
# END


def getsum(finname, foutname, cwdpath, prefix):
    """
    Read syri.out file and output summary statistics
    """
    logger = logging.getLogger("getting summary")
    logger.debug('starting generating summary file')

    # save counts ref_length and qry_length
    stats = {'SYN': [0, 0, 0],
             'INV': [0, 0, 0],
             'TRANS': [0, 0, 0],
             'DUPA': [0, 0, 0],
             'DUPB': [0, 0, 0],
             'NOTALR': [0, 0, 0],
             'NOTALQ': [0, 0, 0],
             'SNP': [0, 0, 0],
             'INS': [0, 0, 0],
             'DEL': [0, 0, 0],
             'CPG': [0, 0, 0],
             'CPL': [0, 0, 0],
             'HDR': [0, 0, 0],
             'TDM': [0, 0, 0]}

    logger.debug('reading syri.out')

    try:
        with open(cwdpath + prefix + finname, 'r') as fin:
            for line in fin:
                line = line.strip().split('\t')
                if line[10] == 'SYN':
                    stats['SYN'][0] += 1
                    stats['SYN'][1] += abs(int(line[1]) - int(line[2])) + 1
                    stats['SYN'][2] += abs(int(line[6]) - int(line[7])) + 1
                if line[10] == 'INV':
                    stats['INV'][0] += 1
                    stats['INV'][1] += abs(int(line[1]) - int(line[2])) + 1
                    stats['INV'][2] += abs(int(line[6]) - int(line[7])) + 1

                if line[10] in ['TRANS', 'INVTR']:
                    stats['TRANS'][0] += 1
                    stats['TRANS'][1] += abs(int(line[1]) - int(line[2])) + 1
                    stats['TRANS'][2] += abs(int(line[6]) - int(line[7])) + 1

                ## Discuss DUP length stat
                if line[10] in ['DUP', 'INVDP']:
                    if line[11] == 'copygain':
                        stats['DUPB'][0] += 1
                        stats['DUPB'][2] += abs(int(line[6]) - int(line[7])) + 1
                    if line[11] == 'copyloss':
                        stats['DUPA'][0] += 1
                        stats['DUPA'][1] += abs(int(line[1]) - int(line[2])) + 1

                if line[10] == 'NOTAL':
                    try:
                        stats['NOTALR'][1] += abs(int(line[1]) - int(line[2])) + 1
                        stats['NOTALR'][0] += 1
                    except ValueError:
                        stats['NOTALQ'][2] += abs(int(line[6]) - int(line[7])) + 1
                        stats['NOTALQ'][0] += 1

                if line[10] == 'SNP':
                    stats['SNP'][0] += 1
                    stats['SNP'][1] += abs(int(line[1]) - int(line[2])) + 1
                    stats['SNP'][2] += abs(int(line[6]) - int(line[7])) + 1

                if line[10] == 'INS':
                    stats['INS'][0] += 1
                    stats['INS'][2] += abs(int(line[6]) - int(line[7]))
                if line[10] == 'DEL':
                    stats['DEL'][0] += 1
                    stats['DEL'][1] += abs(int(line[1]) - int(line[2]))

                # Discuss Copychange statistics
                if line[10] == 'CPG':
                    stats['CPG'][0] += 1
                    # stats['CPG'][1] += abs(int(line[1]) - int(line[2])) + 1
                    stats['CPG'][2] += abs(int(line[6]) - int(line[7])) + 1
                if line[10] == 'CPL':
                    stats['CPL'][0] += 1
                    stats['CPL'][1] += abs(int(line[1]) - int(line[2])) + 1
                    # stats['CPL'][2] += abs(int(line[6]) - int(line[7])) + 1

                if line[10] == 'HDR':
                    stats['HDR'][0] += 1
                    stats['HDR'][1] += abs(int(line[1]) - int(line[2])) + 1
                    stats['HDR'][2] += abs(int(line[6]) - int(line[7])) + 1
                if line[10] == 'TDM':
                    stats['TDM'][0] += 1
                    stats['TDM'][1] += abs(int(line[1]) - int(line[2])) + 1
                    stats['TDM'][2] += abs(int(line[6]) - int(line[7])) + 1
    except FileNotFoundError:
        logger.error(cwdpath + prefix + finname + 'not found.')
        sys.exit()
    except Exception as e:
        logger.error('Error in reading SyRI\'s output' + e)
        sys.exit()

    try:
        with open(cwdpath + prefix + foutname, 'w') as fout:
            fout.write('#Structural annotations\n')
            fout.write('#Variation_type\tCount\tLength_ref\tLength_qry\n')
            fout.write('{}\t{}\t{}\t{}\n'.format('Syntenic regions', stats['SYN'][0], stats['SYN'][1], stats['SYN'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Inversions', stats['INV'][0], stats['INV'][1], stats['INV'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Translocations', stats['TRANS'][0], stats['TRANS'][1], stats['TRANS'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Duplications (reference)', stats['DUPA'][0], stats['DUPA'][1], '-'))
            fout.write('{}\t{}\t{}\t{}\n'.format('Duplications (query)', stats['DUPB'][0], '-', stats['DUPB'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Not aligned (reference)', stats['NOTALR'][0], stats['NOTALR'][1], '-'))
            fout.write('{}\t{}\t{}\t{}\n'.format('Not aligned (query)', stats['NOTALQ'][0], '-', stats['NOTALQ'][2]))

            fout.write('\n\n#Sequence annotations\n')
            fout.write('#Variation_type\tCount\tLength_ref\tLength_qry\n')
            fout.write('{}\t{}\t{}\t{}\n'.format('SNPs', stats['SNP'][0], stats['SNP'][1], stats['SNP'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Insertions', stats['INS'][0], '-', stats['INS'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Deletions', stats['DEL'][0], stats['DEL'][1], '-'))
            fout.write('{}\t{}\t{}\t{}\n'.format('Copygains', stats['CPG'][0], '-', stats['CPG'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Copylosses', stats['CPL'][0], stats['CPL'][1], '-'))
            fout.write('{}\t{}\t{}\t{}\n'.format('Highly diverged', stats['HDR'][0], stats['HDR'][1], stats['HDR'][2]))
            fout.write('{}\t{}\t{}\t{}\n'.format('Tandem repeats', stats['TDM'][0], stats['TDM'][1], stats['TDM'][2]))

    except PermissionError:
        logger.error('Cannot create file' + cwdpath + prefix + foutname + '. Permission denied.')
        sys.exit()
    return 0
# END