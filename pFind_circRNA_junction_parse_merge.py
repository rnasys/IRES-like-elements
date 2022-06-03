#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2018-11-27
E-mail: fanxiaojuan@picb.ac.cn
Description: Merge the circRNA peptide of pFind results
"""

"""
Input format:
#--------------------------------------------------------------------------------------------
File_Name    Scan_No    Exp.MH+    Charge    Q-value    Sequence    Calc.MH+    Mass_Shift(Exp.-Calc.)    Raw_Score    Final_Score    Modification    Specificity    Proteins    Positions    Label    Target/Decoy    Miss.Clv.Sites    Avg.Frag.Mass.Shift    Others
20150410_QE3_UPLC9_DBJ_SA_46fractions_Rep1_46.29366.29366.3.0.dta    29366    3406.619953    3    0    LADQCTGLQGFLVFHSFGGGTGSGFTSLLMER    3406.619321    0.000632    50.953662    5.75208e-016    5,Carbamidomethyl[C];30,Oxidation[M];    3    sp|Q9BQE3|TBA1C_HUMAN/sp|P68363|TBA1B_HUMAN/sp|Q71U36|TBA1A_HUMAN/tr|Q6QMJ5|Q6QMJ5_HUMAN/tr|Q53GA7|Q53GA7_HUMAN/tr|F5H5D3|F5H5D3_HUMAN/tr|F8VQQ4|F8VQQ4_HUMAN/tr|Q9UQM3|Q9UQM3_HUMAN/tr|B7Z1K5|B7Z1K5_HUMAN/tr|B4DQK4|B4DQK4_HUMAN/tr|A8JZY9|A8JZY9_HUMAN/tr|B3KPS3|B3KPS3_HUMAN/    124,K,L/124,K,L/124,K,L/50,K,F/124,K,L/194,K,L/89,K,L/1,G,L/194,K,L/58,K,L/124,K,L/89,K,L/    1|0|0|    target    0    -0.309440    103934    0    0    0    2147483640    0    0    0    36
20150410_QE3_UPLC9_DBJ_SA_46fractions_Rep1_34.23374.23374.3.0.dta    23374    4051.122961    3    0    KLHEEEIQELQAQIQEQHVQIDVDVSKPDLTAALR    4051.109269    0.013692    48.783166    1.0212e-015        3    sp|P08670|VIME_HUMAN/tr|Q53HU8|Q53HU8_HUMAN/tr|B3KRK8|B3KRK8_HUMAN/tr|B0YJC5|B0YJC5_HUMAN/tr|B0YJC4|B0YJC4_HUMAN/tr|V9HWE1|V9HWE1_HUMAN/    235,K,D/235,K,D/222,K,D/53,K,D/235,K,D/235,K,D/    1|    target    2    1.985633    402653183    0    0    0    -134217857    3    0    0    68
"""

"""
Output format:
#---------------------------------------------------------------------------------------------

"""

"""
Note: remove intron or not; strand or strandless
"""

import argparse
import os

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--input-path', dest='p', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/MS_result/MSV000077546/', type=str,help='input file')
    parser.add_argument('-i_junc', '--junc', dest='fnIn_junc', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/circBase_ZFQ_junc_protein_raw.fa', type=str,help='input raw junction dataset')
    parser.add_argument('-ATG', '--ATG', dest='ATG', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/circBase_ZFQ_all_ORF_ATG.fa', type=str,help='input file')
    parser.add_argument('-CTG', '--CTG', dest='CTG', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/circBase_ZFQ_all_ORF_CTG.fa', type=str,help='input file')
    parser.add_argument('-GTG', '--GTG', dest='GTG', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/circBase_ZFQ_all_ORF_GTG.fa', type=str,help='input file')
    parser.add_argument('-TTG', '--TTG', dest='TTG', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/circBase_ZFQ_all_ORF_TTG.fa', type=str,help='input file')
    parser.add_argument('-a', '--all', dest='a', default='F:/experiment/database/circRNA_dataset/hsa_hg19_circRNA.txt', type=str,help='input circBase txt file')
    parser.add_argument('-s', '--small', dest='s', default='F:/experiment/Yangyun/IRES/CIRI2/HEK293_RNase_R.anno.bed', type=str,help='input small circRNA file')
    parser.add_argument('-o', '--outputFile', dest='fnOut', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/MS_result/MSV000077546_pFind.circRNA.junc.merge.txt', type=str,help='output file')
    op=parser.parse_args()
    return op

def dic_build(junc_file):
    """
    Build junction dictionary
    """
    junc_dic = {}
    for i in range(0,len(junc_file),2):
        header = junc_file[i].strip().split('>')[1]
        #print header
        junc_dic[header] = junc_file[i+1].strip()
        #print junc_dic
    return junc_dic

def circPeptide(ms_result,modify_list):
    """
    Extract peptides across back-splicing junction
    """
    circular_list = []
    for i in range(1,len(ms_result)):
        words = ms_result[i].strip().split('\t')
        if len(words) > 12:
            protein_name = words[12]
            q_value = float(words[4])
            #print peptide
            if q_value < 0.01 and protein_name[0:3] == 'chr' and '/REV' not in protein_name:
                modification = words[10].split(';')[:-1]
                #print modification
                modify_count = 0
                if modification != ['']:
                    for modify_types in modification:
                        #print modify_types
                        modify_type = modify_types.split(',')[1]
                        if modify_type not in modify_list:
                            modify_count += 1
                if modify_count != 0:
                    continue
                #print modification
                peptide = words[5]
                header_list = protein_name.split('/chr') #some header lines contain '(n/a)'
                header_list[-1] = header_list[-1][:-1]
                modified_header = [('chr' + x) for x in header_list[1:]]
                header_list = [header_list[0]] + modified_header
                #print header_list
                #print peptide
                for k in range(0,len(header_list)):
                    if header_list[k] not in junc_dic:
                        raise Exception('header error!')
                    else:
                        if peptide in junc_dic[header_list[k]]:
                            #print header_list[k]
                            #print junc_dic[header_list[k]]
                            #print peptide
                            if len(junc_dic[header_list[k]])/2 - 1 > junc_dic[header_list[k]].index(peptide) \
                            and len(junc_dic[header_list[k]])/2 + 1 < junc_dic[header_list[k]].index(peptide)+len(peptide):
                                #print peptide,junc_dic[header_list[j]][n].index(peptide),junc_dic[header_list[j]][n].index(peptide)+len(peptide),len(junc_dic[header_list[j]][n])/2
                                left_aa = len(junc_dic[header_list[k]])/2 - junc_dic[header_list[k]].index(peptide)
                                right_aa = junc_dic[header_list[k]].index(peptide)+len(peptide) - len(junc_dic[header_list[k]])/2
                                circular_list.append(words[0:18]+[','.join([str(left_aa),str(right_aa)])])
    return circular_list

def ORF_dic(ORF_file):
    """
    Build ORF dictionary
    """
    ORFs_dic = {}
    for i in range(0,len(ORF_file),2):
        header = ORF_file[i].strip().split('_')[0].split('>')[1]
        #print header
        ORF_seq = ORF_file[i+1].strip()
        flag = ORF_file[i].strip().split('_')[-1]
        #print flag
        if flag == 'rc':
            ORFs_dic.setdefault(header,[]).append([ORF_seq,'rc'])
        else:
            ORFs_dic.setdefault(header,[]).append([ORF_seq,'not_rc'])
    #print ORFs_dic['chr1:100349894-100366417(+)']
    return ORFs_dic

def circRNA_dic(circRNA_file,flag):
    """
    Build circRNA dic
    """
    cirRNA_dic = {}
    for i in range(0,len(circRNA_file)):
        line = circRNA_file[i].strip().split('\t')
        #print line
        if flag == 0:
            pos = line[0]+':'+line[1]+'-'+line[2]+'('+line[3]+')'
            experiment = line[12]
            gene_name = line[11]
            cell_line = line[7]
        elif flag == 1:
            pos = line[0]+':'+str(int(line[1])+1)+'-'+line[2]+'('+line[5]+')'
            experiment = 'Yuan2016'
            gene_name = line[3]
            cell_line = 'HEK293'
        cirRNA_dic[pos] = [gene_name,experiment,cell_line]
    return cirRNA_dic

def cORF_extract(cPeptide,ATG_dic,CTG_dic,GTG_dic,TTG_dic,cirBase_dic):
    """
    Extract circular ORF from total ORF dataset
    """
    circORF_list = []
    for peptide in cPeptide:
        pep_header = peptide[12].split('@')[0]
        ORF_seq = peptide[5]
        flag = 0
        #print pep_header
        #if pep_header == 'chr1:100349894-100366417(+)':
        if pep_header in ATG_dic:
            #print peptide
            for i in range(len(ATG_dic[pep_header])):
                #print ORFs_dic[pep_header][i]
                if ORF_seq in ATG_dic[pep_header][i][0]:
                    #print ORFs_dic[pep_header][i]
                    flag += 1
                    if pep_header in cirBase_dic:
                        circORF_list.append(peptide+[ATG_dic[pep_header][i][0],str(len(ATG_dic[pep_header][i][0])),\
                                             ATG_dic[pep_header][i][1],'ATG',cirBase_dic[pep_header][0],\
                                             cirBase_dic[pep_header][1],cirBase_dic[pep_header][2]])
                        break
        if flag == 0:
            if pep_header in CTG_dic:
            #print peptide
                for i in range(len(CTG_dic[pep_header])):
                    #print ORFs_dic[pep_header][i]
                    if ORF_seq in CTG_dic[pep_header][i][0]:
                        #print ORFs_dic[pep_header][i]
                        flag += 1
                        if pep_header in cirBase_dic:
                            circORF_list.append(peptide+[CTG_dic[pep_header][i][0],str(len(CTG_dic[pep_header][i][0])),\
                                                 CTG_dic[pep_header][i][1],'CTG',cirBase_dic[pep_header][0],\
                                                 cirBase_dic[pep_header][1],cirBase_dic[pep_header][2]])
                            break
        if flag == 0:
            if pep_header in GTG_dic:
                #print peptide
                for i in range(len(GTG_dic[pep_header])):
                    #print ORFs_dic[pep_header][i]
                    if ORF_seq in GTG_dic[pep_header][i][0]:
                        #print ORFs_dic[pep_header][i]
                        flag += 1
                        if pep_header in cirBase_dic:
                            circORF_list.append(peptide+[GTG_dic[pep_header][i][0],str(len(GTG_dic[pep_header][i][0])),\
                                                 GTG_dic[pep_header][i][1],'GTG',cirBase_dic[pep_header][0],\
                                                 cirBase_dic[pep_header][1],cirBase_dic[pep_header][2]])
                            break
        if flag == 0:
            if pep_header in TTG_dic:
            #print peptide
                for i in range(len(TTG_dic[pep_header])):
                    #print ORFs_dic[pep_header][i]
                    if ORF_seq in TTG_dic[pep_header][i][0]:
                        #print ORFs_dic[pep_header][i]
                        flag += 1
                        if pep_header in cirBase_dic:
                            circORF_list.append(peptide+[TTG_dic[pep_header][i][0],str(len(TTG_dic[pep_header][i][0])),\
                                                 TTG_dic[pep_header][i][1],'TTG',cirBase_dic[pep_header][0],\
                                                 cirBase_dic[pep_header][1],cirBase_dic[pep_header][2]])
                            break
#         if flag == 0:
#             circORF_list.append(peptide+['no_canonical_start_codon','no_canonical_start_codon',\
#                                  'no_canonical_start_codon','no_canonical_start_codon',cirBase_dic[pep_header][0],\
#                                  cirBase_dic[pep_header][1],cirBase_dic[pep_header][2]])
    #print len(circORF_list)
    return circORF_list

def peptide_treat(peptide_file):
    """
    Count each peptide number
    """
    count_dic = {}
    for peptide in peptide_file:
        header = peptide[12]
        count_dic[header]=count_dic.get(header,0)+1
    return count_dic

if __name__ == '__main__':
    op = createHelp()
    
    junc_raw = open(op.fnIn_junc).readlines()
    junc_dic = dic_build(junc_raw)
    print ('Junction dic done!')
     
    print ('Start to extract peptides across back-splicing junction...')
    modify_list = ['Phospho[S]','Phospho[T]','Phospho[Y]','Oxidation[M]','Acetyl[ProteinN-term]',\
                   'Carbamidomethyl[C]','Gln->pyro-Glu[AnyN-termQ]']
    cPeptide = []
    for dirNA in os.listdir(op.p):
        spectra_file = open(op.p+dirNA+'/pFind.spectra').readlines()
    #spectra_file = open(op.p+'MCF7/pFind.spectra').readlines()
        cPeptide += circPeptide(spectra_file,modify_list)
        print (dirNA + ' done.')
    print ('Total extracted peptide number: ' + str(len(cPeptide)))
    
    print ('Start to annotate each peptide...')
    circBase = open(op.a).readlines()
    zfq_circRNA = open(op.s).readlines()
    cirBase_dic = circRNA_dic(circBase,0)
    zfq_dic = circRNA_dic(zfq_circRNA,1)
    #print cirBase_dic
    #Combine cirBase and zfq database 
    for key in zfq_dic.keys():
        if key in cirBase_dic:
            #print key
            #print cirBase_dic[key][1]
            cirBase_dic[key][1] = cirBase_dic[key][1] + ', Yuan2016'
            cirBase_dic[key][2] = cirBase_dic[key][2] + ', HEK293'
        else:
            cirBase_dic[key] = zfq_dic[key]
    #print cirBase_dic['chr5:353844-376831(+)']
    
    print ('Start to add ORF information to each peptide...')
    ATG_file = open(op.ATG).readlines()
    ATG_dic = ORF_dic(ATG_file)
    
    CTG_file = open(op.CTG).readlines()
    CTG_dic = ORF_dic(CTG_file)
    
    GTG_file = open(op.GTG).readlines()
    GTG_dic = ORF_dic(GTG_file)
    
    TTG_file = open(op.TTG).readlines()
    TTG_dic = ORF_dic(TTG_file)
    circORF = cORF_extract(cPeptide,ATG_dic,CTG_dic,GTG_dic,TTG_dic,cirBase_dic)
    
    print ('Start to count each peptide number...')
    peptide_count = peptide_treat(circORF)
      
    out_list = []
    for i in range(0,len(circORF)):
        header_name = circORF[i][12]
        out_list.append('\t'.join(circORF[i] + [str(peptide_count[header_name])]))
    out_list = list(set(out_list))
        
    print ('Total peptide after filter: ' + str(len(out_list)))
    
    output = open(op.fnOut,'w')
    output.write('\n'.join(out_list) + '\n')
    print ('Done!')