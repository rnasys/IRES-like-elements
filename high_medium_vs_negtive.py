#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2017-11-24
E-mail: fanxiaojuan@picb.ac.cn
Description: Correlation figure of high & medium library vs negative library
"""

"""
Input format:
#--------------------------------------------------------------------------------------------
motif    high_medium_count    high_medium_freq    negative_count    negative_freq    fold_change    z-score
GAACGT    35598    0.000130173806349    20470    0.000130197590204    0.999817324922    -0.0208286648758
CGTCGT    31113    0.000113773179306    19969    0.000127011024855    0.89577404352    -12.1455606045
CTTCTT    34894    0.000127599438135    20504    0.000130413844139    0.97841942301    -2.47956834884
"""

"""
Output format (.png):
#---------------------------------------------------------------------------------------------

"""

"""
Note: remove intron or not; strand or strandless
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from matplotlib import rcParams

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input file', dest='fnIn', default='F:/experiment/Yangyun/IRES/6-mer_add/add_pair_end_2bp/IRES_6-mer_high_medium_vs_negative_add.txt', type=str,help='input file')
    parser.add_argument('-o', '--output file', dest='fnOut', default='', type=str,help='output file')
    op=parser.parse_args()
    return op

def data_trim(IRES_peak):
    """
    Extract lines to draw scatter plot
    """
    data_1 = []
    data_2 = []
    for i in xrange(0,len(IRES_peak)):
        words = IRES_peak[i].strip().split('\t')
        ctl_coverage = float(words[2])
        IRES_coverage = float(words[4])
        #print ctl_coverage,IRES_coverage
        data_1.append(ctl_coverage)
        data_2.append(IRES_coverage)
    print pearsonr(data_1,data_2)
    return data_1,data_2

if __name__ == '__main__':
    op = createHelp()
    
    motif_file = open(op.fnIn).readlines()
    high_medium,negative = data_trim(motif_file)
    
    plt.figure(figsize = (5,4))
    ax=plt.subplot(111)
    plt.scatter(high_medium,negative,linewidth=3.0,color='b',alpha=0.7,s=1)
    ax.set_xlim(-0.0001,0.003)
    ax.set_ylim(-0.0001,0.003)
#     plt.xscale('log')
#     plt.yscale('log')
#     ax.tick_params(direction='out', length=6, width=1, labelsize=11)
#     #rcParams.update({'font.family': 'Arial'}) #convert font family of all words
#     ax.set_xticklabels(('$\mathregular{10^{-3}}$','$\mathregular{10^{-2}}$','$\mathregular{10^{-1}}$','$\mathregular{10^0}$',\
#                        '$\mathregular{10^1}$','$\mathregular{10^2}$','$\mathregular{10^3}$','$\mathregular{10^4}$'),\
#                        family='Arial', minor=False) #modify x tick label family
#     ax.set_yticklabels(('$\mathregular{10^{-3}}$','$\mathregular{10^{-2}}$','$\mathregular{10^{-1}}$','$\mathregular{10^0}$',\
#                        '$\mathregular{10^1}$','$\mathregular{10^2}$','$\mathregular{10^3}$',\
#                        '$\mathregular{10^4}$','$\mathregular{10^5}$'),\
#                        family='Arial', minor=False) #modify x tick label family
    ax.set_yticks([0,0.0015,0.003])
    ax.set_xticks([0,0.0015,0.003])
    plt.savefig("F:/experiment/Yangyun/IRES/6-mer_add/add_pair_end_2bp/IRES_6-mer_high_medium_vs_negative_add.png",dpi=800)
    plt.show()