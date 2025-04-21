#!/usr/bin/env python3
import argparse

def updatevcf(args):
    FIN = args.fin.name
    FOUT = args.fout.name

    with open(FIN) as fin, open(FOUT, 'w') as fout:
        for line in fin:
            if line[0] == '': continue
            if line[:2] == '##':
                fout.write(line)
                continue
            if line[:2] == '#C':
                fout.write('\t'.join(line.strip().split('\t')[:-1]) + '\n')
                continue
            line = line.strip().split('\t')
            vec = line[7].split(';')[1]
            assert vec[:8] == 'SUPP_VEC'
            vec = vec.split('=')[1]
            assert len(vec) == 2
            fout.write('\t'.join(line[:8] + ['GT', f'{vec[0]}|{vec[1]}\n']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fin', help='input vcf', type=argparse.FileType('r'))
    parser.add_argument('fout', help='output vcf', type=argparse.FileType('w'))
    args = parser.parse_args()
    updatevcf(args)



