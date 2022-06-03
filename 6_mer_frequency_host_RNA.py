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

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to shuffle sequence'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i_circRNA', '--input-file', dest='fnIn_cir', default='H:/Yangyun/IRES/transcript_ires/region_motif/IRES.mature.utr5.fa', help='input file')
    parser.add_argument('-l', '--motif-length', dest='motif_len', default=6, help='input motif length')
    parser.add_argument('-enrich', '--enrich-motif', dest='enrich', default='H:/Yangyun/IRES/6-mer_add/add_pair_end_2bp/IRES_6-mer_high_medium_vs_negative_add_enrich_motif.txt', help='input motif length')
    parser.add_argument('-deplete', '--deplete-motif', dest='deplete', default='H:/Yangyun/IRES/6-mer_add/add_pair_end_2bp/IRES_6-mer_high_medium_vs_negative_add_deplete_motif.txt', help='output file')
    parser.add_argument('-e', '--enrich-freq', dest='e', default='H:/Yangyun/IRES/6-mer_add/add_pair_end_2bp/enrich_motif_frequency/enrich_motif_frequency_CITIs_utr5_linear.txt', help='output file')
    parser.add_argument('-d', '--deplete-freq', dest='d', default='H:/Yangyun/IRES/6-mer_add/add_pair_end_2bp/enrich_motif_frequency/deplete_motif_frequency_CITIs_utr5_linear.txt', help='output file')
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
                    for j in range(i+1,len(sequence_file)):
                        if '>' in sequence_file[j]:
                            i=j-1
                            break
                        else:
                            my_seq.append(sequence_file[j].upper().strip())
                    #print my_seq[0]
                    my_seq = ''.join(my_seq[0])
                    #print my_seq
                    my_final_seq = ''.join(str(my_seq))
                    sequence.append(my_final_seq)
        i += 1
    return sequence

def motif_count(sequence,motif_dic):
    motif_frequency_dic = {}
    for i in range(0,len(sequence)):
        motif_count_dic = {}
        motif_count_dic = motif_count_dic.fromkeys(motif_dic.keys(),0)
        #print len(motif_count_dic)
        for n in range(0,len(sequence[i])-5):
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

def motif_dic_build(motif_file):
    motif_dic = {}
    for motif in motif_file:
        motif_seq = motif.strip()
        motif_dic[motif_seq] = ''
    return motif_dic

if __name__ == '__main__':
    time_start = time.time()
    op = createHelp()
    sequence_file = open(op.fnIn_cir).readlines()
    sequence = seq_treatment(sequence_file)
    print ('Total circRNA sequence: ' + str(len(sequence)))
    
    a = ['A','T','C','G']
    b = 6
    motif_list = [''.join(x) for x in itertools.product(*[a] * b)]
    motif_dic = {}
    for i in range(0,len(motif_list)):
        motif_dic[motif_list[i]] = 0
    print ('Total possible motifs length: ' + str(len(motif_list)))
    
    motif_number_dic = motif_count(sequence, motif_dic)
    
    enrich_motif = open(op.enrich).readlines()
    deplete_motif = open(op.deplete).readlines()
    enrich_dic = motif_dic_build(enrich_motif)
    deplete_dic = motif_dic_build(deplete_motif)
    
    enrich_list = []
    deplete_list = []
    
    for key in motif_number_dic:
        if key in enrich_dic:
            enrich_list.append('\t'.join([key,str(motif_number_dic[key])]))
        if key in deplete_dic:
            deplete_list.append('\t'.join([key,str(motif_number_dic[key])]))
    
    print ('Enriched list: ' + str(len(enrich_list)))
    print ('Deplete list: ' + str(len(deplete_list)))
    
    out_enrich = open(op.e,'w')
    out_deplete = open(op.d,'w')
    out_enrich.write('\n'.join(enrich_list) + '\n')
    out_deplete.write('\n'.join(deplete_list) + '\n')
    
    print ('Time cost: ' + str(time.time() - time_start))
    print ('Well done!')