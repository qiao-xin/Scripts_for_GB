#!/usr/bin/python
# -*- coding: UTF-8 -*- 

# Xin Qiao, 11/16/2016
# Xin Qiao, 01/16/2019

import os
import sys
from sys import argv

import pandas as pd

from scipy import stats

#Loading expression values for each gene under different conditions
rpkm = pd.read_table(argv[1], sep="\t", header=0, index_col='GeneID')#, header=None)

#Loading gene pairs
pairs = pd.read_table(argv[2], sep="\t", header=None)

c = {}
for k in rpkm.index:
    #print k
    c[k] = 1

d = 0
infile = open(argv[2])
for line in infile.readlines():
    d = d + 1
#print d

a = []
h = 0
g = ''
sm = 0
outfile = open(argv[3], 'w')
f = 'Duplicate1\tDuplicate2\tpearsonr\tp-value\n'
outfile.write(f)
for j in range(0, d):
    for i in pairs.loc[j]:
        #print i
        if not (c.has_key(i)):
            continue
        sm = rpkm.loc[i].sum()
        if sm == 0:
            continue
        else:
            h = h + 1
            info = rpkm.loc[i]
            a.append(info)
            if h == 1:
                g = g + i + '\t'
                #print g
            else:
                g = g + i
                #print g
    alen = len(a)
    if alen < 2:
        g = ''
        a = []
        h = 0
        continue
    else:
        #print len(a), a[0], a[1]
        #print g
        r, p_value = stats.pearsonr(a[0], a[1])
        #print r,"\t",p_value
        b = '%s\t%s\t%s\n' % (g,r,p_value)
        outfile.write(b)
        a = []
        g = ''
        h = 0
outfile.close