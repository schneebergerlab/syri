#!/usr/bin/env python3
import argparse

if __name__=='__main__':
    parser = argparse.ArgumentParser('')
    parser.add_argument('fa', help='name of fasta file', type=argparse.FileType('r'))
    parser.add_argument('seqnames', help='file containing names sequence to complement', type=argparse.FileType('r'))
    args = parser.parse_args()

    ids = [line.strip() for line in args.seqnames.readlines()]
    # from syri.scripts.func import readfasta, revcomp
    from func import readfasta, revcomp
    # from Bio.SeqIO import parse
    # from Bio.Seq import reverse_complement
    import re

    with open(args.fa.name+'.rc', 'w') as fout:
        for chrid, seq in readfasta(args.fa.name).items():
            fout.write('>'+chrid+'\n')
            if chrid in ids:
                fout.write(re.sub("(.{64})", "\\1\n", str(revcomp(seq)), 0, re.DOTALL))
            else:
                fout.write(re.sub("(.{64})", "\\1\n", str(seq), 0, re.DOTALL))
            fout.write('\n')
