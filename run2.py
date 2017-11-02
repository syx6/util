#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Pw @ 2017-11-02 16:44:42

import sys
import pysam

COMPLEMENT = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
}

def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def strandCheck(seq, strand):
    if strand == '-':
        resultSeq = complement(seq)
    else:
        resultSeq = seq
    return resultSeq

def main(file, fa):
    print('Reading hg19.bed')
    genome_fa = pysam.FastaFile(fa)
    with open(file, 'r') as f:
        for i in f:
            chrom = i.split('\t')[0]
            startPoint = int(i.split('\t')[1])
            endPoint = int(i.split('\t')[2])
            strand = i.strip().split('\t')[-1]
            cir = i.split('\t')[3]
            if startPoint <= endPoint:
                if endPoint - startPoint < 100:
                    midPoint = ((endPoint - startPoint) / 2) + startPoint
                    seq = genome_fa.fetch(chrom, midPoint, endPoint)
                    seq += genome_fa.fetch(chrom,startPoint, (midPoint-1))
                    resultSeq = strandCheck(seq, strand)
                else:
                    seq = genome_fa.fetch(chrom, (endPoint-49), endPoint)
                    seq += genome_fa.fetch(chrom, startPoint, (startPoint+49))
                    resultSeq = strandCheck(seq, strand)
            print(cir + '\t' + resultSeq)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
