#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Author: Xin Qiao
# Date: 01/02/2018
# Date: 10/12/2018
import sys
from sys import argv

try:
    argv[1]
    argv[2]
    argv[3]
except:
    print "Usage: python %s <Species_Abbr.synteny.blocks.ks.info> <components> <Species_Abbr>" % argv[0]
    sys.exit()

import os
import re
import pandas as pd
import numpy as np
from numpy import log
#from sklearn.mixture import GMM, GaussianMixture
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt





def gaussian(fi,component):
    data = pd.read_table(fi)
    sort_data = data.sort_values(['Average Ks'], ascending=True)
    aks = sort_data['Average Ks']
    aks = np.array(aks)
    
    line_distribution   = np.random.choice(a = aks, size = 100000)
    number_points       = len(line_distribution)
    
    gmm = GaussianMixture(n_components = component)
    gmm.fit(np.reshape(line_distribution, (number_points, 1)))
    gauss_mixt = np.array([p * norm.pdf(aks, mu, sd) for mu, sd, p in zip(gmm.means_.flatten(), np.sqrt(gmm.covariances_.flatten()), gmm.weights_)])
    gauss_mixt_t = np.sum(gauss_mixt, axis = 0)
    return sort_data,aks,gmm,gauss_mixt

def hist_plot(sort_data):
    histd = sort_data['Average Ks'].values
    l = len(histd)
    nb = int(10*log(l))
    grey  = [0.9, 0.9, 0.9]
    ax.hist(histd, bins = nb, normed = True, color = '#f1f1f1', edgecolor = '#B8BFC6', label = "Ks distribution")#Empirical

def gaussian_fit(aks,gauss_mixt,gmm):
    color=['#7bb3ff','#fb6a6a','#96cd69','#ffcc5c']
    for i in range(len(gauss_mixt)):
        ax.plot(aks, gauss_mixt[i], color=color[i], label = 'Component ' + str(i+1))
    ax.axes.spines['right'].set_color('none')
    ax.axes.spines['top'].set_color('none')
    ax.spines['left'].set_position(('axes', -0.01))
    ax.spines['bottom'].set_position(('axes', -0.01))
    
    plt.setp(ax.get_xticklabels(), fontsize=10, fontname="Helvetica")
    plt.setp(ax.get_yticklabels(), fontsize=10, fontname="Helvetica")
    
    ax.set_xlabel('Ks', fontdict)#Synonymous substitutions per site (Ks)
    ax.set_ylabel('The frequency of syntenic blocks', fontdict)
    
    j=0
    #sd means Standard deviation, the smaller sd is better; p means weights or the proportion of each component
    for mu, sd, p in zip(gmm.means_.flatten(), np.sqrt(gmm.covariances_.flatten()), gmm.weights_):
        j = j + 1
        print "Component %s" % j
        #print mu,sd,p
        print "%s\t%s\t%s" % (mu,sd,p)

fig = plt.figure(figsize=(6,4))
fontdict = {"color":"k", "size":10, 'fontname':'Helvetica'}
fontdict2 = {"color":"k", "size":9, 'fontname':'Helvetica'}

comp = int(float(argv[2]))
sort_data,aks,gmm,gauss_mixt = gaussian(argv[1],comp)

ax = fig.add_subplot(111)
hist_plot(sort_data)
gaussian_fit(aks,gauss_mixt,gmm)
plt.ylim(0,6)
plt.xlim(0,5)
F = ax.legend(loc='best', fontsize=10, frameon=False)
plt.setp(F.texts, family='Helvetica')

fgname = '%s.synteny.blocks.ks.distri.pdf' % argv[3]
plt.savefig(fgname)
