#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Pw @ 2017-10-25 10:35:00

import sys
import scipy.stats as stats

def main(file1, file2):
    counts = {}
    cir_srptm = {}
    with open(file1, 'r') as f1:
        for i in f1:
            if i.startswith('sample'):
                pass
            else:
                gene = i.split('\t')[0]
                hcounts = [float(n) for n in i.strip().split('\t')[1:]]
                counts[gene] = hcounts
    with open(file2, 'r') as f2:
        for i in f2:
            cir = i.split('\t')[0].split('_')[0]
            srptm = [float(n) for n in i.strip().split('\t')[1:]]
            cir_srptm[cir] = srptm
    with open(file2, 'r') as f:
        for i in f:
            cir = i.split('\t')[0].split('_')[0]
            gene = i.split('\t')[0].split('_')[1]
            if gene in counts:
                if cir in cir_srptm:
                    linear_gene = counts[gene]
                    circrna = cir_srptm[cir]
                    print cir + '\t' + gene + '\t' + str(stats.pearsonr(linear_gene, circrna)[0]) + '\t' + str(stats.pearsonr(linear_gene, circrna)[1])
                    #cal correlation coefficient

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
