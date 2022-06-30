#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:05:51 2017

@author: goel

Current VCF handling tools seems to be cater VCFs generated from read alignment.
VCF generated from assemblies, have other requirement.
This package consists of methods to manipulate VCFs generated from the comparison
of haplotype-resolved assemblies
"""

import argparse

import pysam.libcbcf


def setlogconfig(lg):
    # TODO: Setup one logging function for syri and other utilities in func.py
    import logging.config
    logging.config.dictConfig({
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'log_file': {
                'format': "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s",
            },
            'stdout': {
                'format': "%(name)s - %(levelname)s - %(message)s",
            },
        },
        'handlers': {
            'stdout': {
                'class': 'logging.StreamHandler',
                'formatter': 'stdout',
                'level': 'WARNING',
            },
        },
        'loggers': {
            '': {
                'level': lg,
                'handlers': ['stdout'],
                # 'handlers': ['stdout', 'log_file'],
            },
        },
    })
#END



def hapmerge(vcf1, vcf2, outvcf, sample):
    '''
    :param vcf1: Input VCF (also BCF?) for haplotype 1
    :param vcf2: Input VCF (also BCF?) for haplotype 2
    :param sample: Sample name to be added to VCF
    Merges VCF files from different haplotypes of a sample.
    Current limitations:
        1) Only merges two VCFs (i.e. for a diplid genome)
        2) Only works with SNPs and INDELs                     
    :return: VCF (also BCF?) file where the SNPs/INDELs are merged with GT tag. Var IDs and INFO fields would be deleted. FORMAT and sample would be added. Positions with more than 1 SNP/indels in any VCF would be filtered out in both VCFs.
    '''

    def readvcf(v):
        '''
        Read VCF file
        :param v: VariantFile object
        :return: VCF data for unique and list of duplicated positions
        '''
        cids = deque()      # Chromosome ids
        lastchr = ''
        vd = defaultdict(deque)
        for var in v.fetch():
            if var.chrom != lastchr:
                if var.chrom not in cids:
                    cids.append(var.chrom)
            vd[var.chrom].append((var.pos, var.ref, var.alts))
        gp = deque()            # garb pos
        v1d = defaultdict(dict)
        for k in vd:
            vd[k] = deque(set(vd[k]))
            garb = Counter([v[0] for v in vd[k]])
            garb = set([k for k, v in garb.items() if v > 1])
            for g in garb:
                gp.append((k, g))
            n = {v[0]: (v[1:]) for v in vd[k] if v[0] not in garb}     # new var list
            v1d[k] = n
        return v1d, gp, cids


    # TODO: Perform the same operation on the syri.out file as well
    from pysam import VariantFile, VariantRecord
    import logging
    from collections import deque, defaultdict, Counter

    setlogconfig(args.log)
    logger = logging.getLogger("hapmerge")

    v1 = VariantFile(vcf1, 'r')
    v2 = VariantFile(vcf2, 'r')

    # if v1.header != v2.header:
    #     logger.info(f"Input VCFs have different headers. Using the header from {vcf1} in the output file.")
    header = v1.header
    header.formats.add(id='GT', number=1, type='String', description="Genotype after merging haplotype calls")
    header.samples.add(sample)

    # Read first vcf, save chromosome order, position and variations
    v1dict, garbv1, chridsv1 = readvcf(v1)
    v2dict, garbv2, chridsv2 = readvcf(v2)
    # TODO: Add IDs as well (can be added as a semicolon separated list)

    # Remove garb positions
    garbv1.extend(garbv2)
    for v in garbv1:
        try: v1dict[v[0]].pop(v[1])
        except KeyError: pass
        try: v2dict[v[0]].pop(v[1])
        except KeyError: pass

    vo = VariantFile(outvcf, 'w', header=header)
    for k in chridsv1:
        print(k)
        chrom = k
        pos = sorted(set(v1dict[k]).union(v2dict[k]))
        for p in pos:
            count += 1
            ref1, ref2, alt1, alt2 = ['']*4
            try:
                ref1 = v1dict[k][p][0]
                alt1 = v1dict[k][p][1][0]
            except KeyError: pass
            try:
                ref2 = v2dict[k][p][0]
                alt2 = v2dict[k][p][1][0]
            except KeyError: pass
            if ref1 == '' and ref2 == '':
                logger.error(f"Both ref are empty for position {k}:{p}. Error. Skipping this position.")
                continue
            elif ref1 == '':
                # Heterozygous mutations present only in paternal chromosome
                ref, alt, gt = ref2, alt2, (0, 1)
            elif ref2 == '':
                # Heterozygous mutations present only in maternal chromosome
                ref, alt, gt = ref1, alt1, (1, 0)
            elif ref1 == ref2:
                ref = ref1
                # Homozygous mutation
                if alt1 == alt2:
                    alt = alt1
                    gt = (1, 1)
                # Multi-allelic variation involving only SNPs, insertions
                else:
                    alt = (alt1, alt2)
                    gt = (1, 2)
            else:
                # Mutli-allelic variation involving deletions
                loge = f"Incompatible reference alleles at {k}:{p}: ref1-{ref1}, ref2-{ref2}. Skipping it."
                if len(ref1) == len(ref2):
                    logger.error(loge)
                    continue
                if len(ref1) > len(ref2):
                    if ref1[:len(ref2)] != ref2:
                        logger.error(loge)
                        continue
                    ref, alt = ref1, (alt1, alt2+ref1[len(ref2):])
                else:
                    if ref1 != ref2[:len(ref1)]:
                        logger.error(loge)
                        continue
                    ref, alt = ref2, (alt1+ref2[len(ref1):], alt2)
                gt = (1, 2)
            alleles = (ref, alt) if type(alt) is not tuple else (ref, alt[0], alt[1])
            r = vo.new_record(contig=chrom, start=p-1, stop=p+len(ref)-1, alleles=alleles, id='.', filter='PASS')
            r.samples[sample]['GT'] = gt
            vo.write(r)
    vo.close()
# END




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # TODO: setup arguments for the hapmerge function
    subparsers = parser.add_subparsers()

    other.add_argument('--log', help='Log-level', choices=['DEBUG', 'INFO', 'WARN'], default='WARN', type=str)


    parser_region = subparsers.add_parser("region", help="get annotation for the region")
    parser_getbed = subparsers.add_parser("getbed", help="get bed file for the regions")
    parser_snpanno = subparsers.add_parser("snpanno", help="get annotation for SNPs list")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    parser_region.set_defaults(func=getRegion)
    parser_region.add_argument("chr", help="Chromosome ID", type=str)
    parser_region.add_argument("start", help="region start", type=int)
    parser_region.add_argument("end", help="region end", type=int)
    parser_region.add_argument("-q", help="search region in query genome", action="store_true", default=False)
    parser_region.add_argument("-f", help="files to search", type=argparse.FileType('r'), nargs="+", default=None)
