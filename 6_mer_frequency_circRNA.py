#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2015-10-19
E-mail: fanxiaojuan@picb.ac.cn
Description: hexamer 
"""

import random
import argparse
import numpy as np
import itertools
import time

def definition(m6a_seq_,sequrnce_):
    global m6a_seq
    global sequence
    m6a_seq = m6a_seq_
    sequence = sequrnce_

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to shuffle sequence'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i_circRNA', '--input-file', dest='fnIn_cir', default='F:/experiment/Yangyun/IRES/CIRI2/HEK293_RNase_R_merge_full_circRNA_filter.fa', help='input file')
    parser.add_argument('-l', '--motif-length', dest='motif_len', default=6, help='input motif length')
    parser.add_argument('-o', '--output shuffle sequence', dest='fnOut', default='F:/experiment/Yangyun/IRES/6-mer_add/add_pair_end_2bp/enrich_motif_frequency/hexamer_frequency_circRNA_Yuan2016.txt', help='output file')
    op=parser.parse_args()
    return op

def seq_treatment(sequence_file):
    """
    Clean DNA seuqences.
    """
    sequence = []
    #print sequence_file
    i = 0
    while i<len(sequence_file):
        #print sequence_file[i]
        if '>' in sequence_file[i]:
            if 'N' not in sequence_file[i+1]:
                if 'n' not in sequence_file[i+1]:
                    #print sequence_file[i+1]
                    my_seq = []
                    for j in xrange(i+1,len(sequence_file)):
                        if '>' in sequence_file[j]:
                            i=j-1
                            break
                        else:
                            my_seq.append(sequence_file[j].upper().strip())
                    #print my_seq[0]
                    my_seq = ''.join(my_seq[0])
                    #print my_seq
                    my_final_seq = ''.join(str(my_seq) + my_seq[0:op.motif_len-1])
                    sequence.append(my_final_seq)
        i += 1
    return sequence

def motif_count(sequence,motif_dic):
    motif_frequency_dic = {}
    for i in xrange(0,len(sequence)):
        motif_count_dic = {}
        motif_count_dic = motif_count_dic.fromkeys(motif_dic.keys(),0)
        #print len(motif_count_dic)
        for n in xrange(0,len(sequence[i])-5):
            #print sequence[i]
            motif_count_dic[sequence[i][n:n+6]]=motif_count_dic.get(sequence[i][n:n+6],0)+1
            #print motif_count_dic
        for key in motif_count_dic:
            motif_frequency_dic.setdefault(key,[]).append(motif_count_dic[key]/(len(sequence[i])*1.0))
#         if motif_count_dic == {}:
#             print sequence[i]
    for motif in motif_frequency_dic:
        #print len(motif_frequency_dic[motif])
        motif_frequency_dic[motif] = np.mean(motif_frequency_dic[motif])
    return motif_frequency_dic

if __name__ == '__main__':
    time_start = time.time()
    op = createHelp()
    sequence_file = open(op.fnIn_cir).readlines()
    sequence = seq_treatment(sequence_file)
    print 'Total circRNA sequence: ' + str(len(sequence))
    
    a = ['A','T','C','G']
    b = 6
    motif_list = [''.join(x) for x in itertools.product(*[a] * b)]
    motif_dic = {}
    for i in xrange(0,len(motif_list)):
        motif_dic[motif_list[i]] = 0
    print 'Total possible motifs length: ' + str(len(motif_list))
    
    motif_number_dic = motif_count(sequence, motif_dic)
    
    output = open(op.fnOut,'w')
    for motif in motif_number_dic.keys():
        output.write(motif + '\t' + str(motif_number_dic[motif]) + '\n')
    print 'Time cost: ' + str(time.time() - time_start)
    print 'Well done!'