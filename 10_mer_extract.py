#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2016-3-23
E-mail: fanxiaojuan@picb.ac.cn
Description: Extract 10 mer element from fastq file
"""

"""
Input file:
#------------------------------------------------------------------------------ 
@HWI-D00599:84:C8A6GANXX:1:1101:1174:2176 1:N:0:CTTGTA
TCGGCATGGACGAGCTGTACAAGTAAAGAGCACAGAATCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACG
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFBFFFFFB/<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@HWI-D00599:84:C8A6GANXX:1:1101:1135:2184 1:N:0:CTTGTA
TCGGCATGGACGAGCTGTACAAGTAACGACTAACTGATCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACG
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#------------------------------------------------------------------------------ 
"""

"""
Output file (10-mer element):
#------------------------------------------------------------------------------ 
AGAGCACAGA
GAATATCCAG
#------------------------------------------------------------------------------ 
"""

import argparse
import re
# from compiler.ast import flatten
#from multiprocessing import Pool

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--input fastq', dest='fnIn', default='F:/experiment/Yangyun/IRES/10_mer/test.fastq', type=str,help='input fastq file')
    parser.add_argument('-f', '--input forward primer', dest='f', default='TGTACAAGTAA', type=str,help='input forward primer seq')
    parser.add_argument('-r', '--input reverse primer', dest='r', default='ATCATGGTGAG', type=str,help='input reverse primer seq')
    parser.add_argument('-o', '--output element', dest='fnOut', default='F:/experiment/Yangyun/IRES/10_mer/test_10-mer.txt', type=str,help='output element file')
    op=parser.parse_args()
    return op

def reads_extract(fastq):
    """
    Extract reads from fastq file
    """
    reads_list = []
    for i in range(0,len(fastq)-3,4):
        reads_list.append(fastq[i+1].strip())
    return reads_list

def element_extract(reads):
    """
    Extract element based forward and reverse primers
    """
    element_list = []
    count_2 = 0
    count_1 = 0
    count_0 = 0
    for i in range(0,len(reads)):
        if len(re.findall(r'(?<=GGCATGGACGAGCTGTACAAGTAA)[\S]{8,12}(?=ATCATGGTGAGCAAGGGCGAGGAGCTGTTC)', reads[i])) > 1:
            print (reads[i])
            count_2 += 1
        if len(re.findall(r'(?<=GGCATGGACGAGCTGTACAAGTAA)[\S]{8,12}(?=ATCATGGTGAGCAAGGGCGAGGAGCTGTTC)', reads[i])) == 0:
            count_0 += 1
        else:
            count_1 += 1
            element_list.append(re.findall(r'(?<=GGCATGGACGAGCTGTACAAGTAA)[\S]{8,12}(?=ATCATGGTGAGCAAGGGCGAGGAGCTGTTC)', reads[i]))
    print ('There are ' + str(count_2) + ' reads with more than 1 element.')
    print ('There are ' + str(count_1) + ' reads with 1 element.')
    print ('There are ' + str(count_0) + ' reads without element.')
    element_list = [x for x in sub for sub in element_list]
    return element_list

if __name__ == '__main__':
    op = createHelp()
    
    fastq_file = open(op.fnIn).readlines()
    reads_info = reads_extract(fastq_file)
    print ('There are ' + str(len(fastq_file)) + ' lines in fastq file.')
    print ('There are ' + str(len(reads_info)) + ' reads in fastq file.')
    
    #p = Pool(2)
    #element_info = p.map(element_extract,reads_info)
    element_info = element_extract(reads_info)
    print ('We got ' + str(len(element_info)) + ' elements')
    
    motif_length_dic = {}
    for i in range(0,len(element_info)):
        #print element_info[i]
        motif_length_dic[len(element_info[i])] = motif_length_dic.get(len(element_info[i]),0) + 1
    print (motif_length_dic)

    output_file = open(op.fnOut,'w')
    for i in range(0,len(element_info)):
        output_file.write(element_info[i][0] + '\n')
        
    print ('Done!')