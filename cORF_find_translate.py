#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2018-3-15
E-mail: fanxiaojuan@picb.ac.cn
Description: Compare circRNA ORF with host ORF
"""

"""
Input format:
#--------------------------------------------------------------------------------------------

"""

"""
Output format:
#---------------------------------------------------------------------------------------------

"""

"""
Note: remove intron or not; strand or strandless
"""

import argparse
import re
import time

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-fa', '--circRNA-seq', dest='fa', default='H:/Yangyun/IRES/CIRI2/open_pFind_full/CIRI2/merge.AS.ciri.uniq.fa', type=str,help='input file')
    parser.add_argument('-l', '--ORF-length', dest='l', default=20, type=int, help='input file')
    parser.add_argument('-s', '--start-codon', dest='s', default='TTG', type=str, help='input file')
    parser.add_argument('-o', '--output file', dest='fnOut', default='H:/Yangyun/IRES/CIRI2/open_pFind_full/CIRI2/merge.AS.ciri.uniq.TTG.orf.fa', type=str,help='output file')
    #parser.add_argument('-o_ctl', '--output-ctl', dest='fnOut_ctl', default='F:/experiment/Yangyun/IRES/CIRI2/HEK293_RNase_R_reverse_ctl.orf', type=str,help='output file')
    op=parser.parse_args()
    return op

def translate_aa(orf_seq,aa_dic):
    """
    Translate RNA to protein
    """
    aa_codon = [orf_seq[x:x+3] for x in range(0,len(orf_seq),3)]
    #print aa_codon
    aa_seq = ''.join([aa_dic[x] for x in aa_codon])
    return aa_seq

def ORF_extraction_frame(stop_index,rna_seq,op,frame):
    """
    Extract ORF from sequence
    """
    orf_list = []
    if stop_index == []:
        #print stop_index
        start_index = [m.start() for m in re.finditer(op.s, rna_seq)]
        if start_index != []:
            #print start_index
            for y in start_index:
                if y % 3 == frame:
                    #print rna_seq
                    #print y
                    orf_seq = rna_seq[y:int((len(rna_seq)-y)/3)*3+y]
                    #print orf_seq
                    #print start_index
                    orf_aa = translate_aa(orf_seq,amino_acid_dic)
                    orf_list.append(orf_aa+'rc')
                    #print orf_aa
                    break
    else:
        for x in range(len(stop_index)-1):
            start_index = [m.start() for m in re.finditer(op.s, rna_seq[(stop_index[x]+2):(stop_index[x+1])])]
            #print stop_index[x],stop_index[x+1],start_index
            #print rna_seq[(stop_index[x]+2):(stop_index[x+1])]
            for y in start_index:
                if (stop_index[x+1]-(stop_index[x]+2 + y)) % 3 == 0:
                    orf_seq = rna_seq[(stop_index[x]+2 + y):stop_index[x+1]+3]
                    #print orf_seq
                    if len(orf_seq) >= (op.l+1)*3: #20aa + stop codon
                        #print orf_seq
                        orf_aa = translate_aa(orf_seq,amino_acid_dic)
                        #print orf_seq
                        #print orf_aa
                        orf_list.append(orf_aa)
                        break
    orf_list = list(set(orf_list))
    #print orf_list
    return orf_list

def orf_detection(fa_file,flag):
    """
    Detect circRNA ORF overlapping with host ORF
    """
    temp_length = 0
    orf_list = []
    out_list = []
    for i in range(0,len(fa_file),2):
        header_line = fa_file[i].strip()
        #header_line = fa_file[i].strip().split('(')[0] + '||' + fa_file[i].strip().split('@>')[1].split('||')[0]
        #print header_line
        if flag == '0':
            seq_line = fa_file[i+1].strip().upper()
            #print seq_line
        elif flag == '1':
            seq_line = fa_file[i+1].strip()[::-1].upper()
            #print seq_line
        if 'N' not in seq_line:
            #print seq_line
            seq_line = seq_line * 4
            #print seq_line
            stop_index = [m.start() for m in re.finditer('TAA|TAG|TGA', seq_line)]
            #print stop_index
            stop_index_f1 = [x for x in stop_index if x % 3 == 0]
            stop_index_f2 = [x for x in stop_index if x % 3 == 1]
            stop_index_f3 = [x for x in stop_index if x % 3 == 2]
            orf_f1 = ORF_extraction_frame(stop_index_f1,seq_line,op,0)
            orf_f2 = ORF_extraction_frame(stop_index_f2,seq_line,op,1)
            orf_f3 = ORF_extraction_frame(stop_index_f3,seq_line,op,2)
            orf_list = list(set(orf_f1+orf_f2+orf_f3))
            #print orf_list
            for j in range(len(orf_list)):
                if orf_list[j][-2:] == 'rc':
                    #print out_list[j]
                    out_list.append(header_line+'_'+str(j)+'_rc')
                    out_list.append(orf_list[j][:-2])
                else:
                    out_list.append(header_line+'_'+str(j))
                    out_list.append(orf_list[j])
                    #print out_list
            #print out_list
    #print out_list
    return out_list

if __name__ == '__main__':
    time_start = time.time()
    op = createHelp()
    
    amino_acid_dic = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"", "TAG":"",
    "TGT":"C", "TGC":"C", "TGA":"", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    fa_file = open(op.fa).readlines()
    cir_orf = orf_detection(fa_file,'0')
#     cir_orf_ctl = orf_detection(fa_file,'1')
    
    #print cir_orf
    output = open(op.fnOut,'w')
    output.write('\n'.join(cir_orf) + '\n')
    
#     output_ctl = open(op.fnOut_ctl,'w')
#     output_ctl.write('\n'.join(cir_orf_ctl))
    
    print ('Done!')
    print ('Time used: ' + str(time.time() - time_start))