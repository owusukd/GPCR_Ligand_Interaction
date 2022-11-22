#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 09:42:46 2022
This code get all the Amino Acids of a given pocket
It goes through all the ligand folders and all the proteins in the folder and all their 
predicted binding pockets
It create a tsv file at the end
@author: kwabena
"""

import os
import re
import csv
path = '/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control/3D_Comparison_Pocket_overlap_scores/Overlap_scores'
os.chdir(path)

protein = [["0HK_P08173_A_5dsg","0HK_P11229_A_5cxv"], ["7LD_P28223_A_6wgt","7LD_P41595_A_5tvn"],
           ["7MA_O43613_A_6tod","7MA_O43614_A_5wqc"], ["8NU_P14416_A_6cm4","8NU_P28223_A_6a93"],
           ["40F_Q14416_C_4xaq","40F_Q14832_C_4xar"], ["89F_P28222_A_5v54","89F_P28223_A_6wh4"],
           ["ADN_P29274_A_2ydo","ADN_P30542_A_6d9h_R"],["GGL_O00222_C_6bsz","GGL_Q14416_C_5cni","GGL_Q14832_C_5cnk"],
           ["GLU_A0A173M0G2_TR2_5x2p_A","GLU_E9P5T5_NOT_4io2_A","GLU_P41594_C_3lmk_A","GLU_P42264_NOT_3s9e_A","GLU_Q14416_C_7mtr_A"],
           ["SRO_P08908_A_7e2y_R","SRO_P37231_NOT_3adv_B"],
           ["Z99_P41594_C_7fd9_A","Z99_Q13255_C_3ks9_A","Z99_Q14831_C_3mq4","Z99_Q14832_C_7wi6_A"]]

folder = ["0HK", "7LD", "7MA", "8NU", "40F", "89F", "ADN", "GGL", "GLU", "SRO", "Z99"]

pkt_num = [[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,2,1],[1,1,1,1,2],[1,1],[1,5,1,10]]

columns = ['ligand_ID', 'Protein_PDB_ID', 'Pocket', 'AA_pkt']

results_file = 'AA_Pkt.tsv'
# create a csv file with only column heads
with open(results_file, "w") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(columns)
    
regex = r"[\sA-Z]{1,7}\d{1,}(?=\s{4})"

for indx, fold in enumerate(folder):
    path_list = [path, fold]
    paths = "/".join(path_list)
    os.chdir(paths)
    
    all_columns = {}
    
    for i in range(len(columns)):
        all_columns.setdefault(columns[i],[])

    for j in range(len(protein[indx])):
        #pkt_num = (pkt_counts[indx])[j]+1
        #for k in range(1,pkt_num):
            pkt_name = "_".join([(protein[indx])[j], 'pkt', str((pkt_num[indx])[j])])
            pkt_file_name = ".".join([pkt_name, 'txt'])
            pkt_file = open(pkt_file_name)
            pkt = pkt_file.read()
            
            AA_pkt = re.finditer(regex, pkt)
            AA_pkt_uniq = []
            AA_pkt_num = []
            for aa_pkt in AA_pkt:
                aa_pk = aa_pkt.group().strip()
                aa_pk = aa_pk.split(" ")
                aa_code = aa_pk[0]
                aa_num = int(re.sub(r"\D","", aa_pk[-1]))
                if aa_num not in AA_pkt_num:
                    AA_pkt_num.append(aa_num)
                    AA_pkt_uniq.append(aa_code)
            
            AA_pkt_uniq = ",".join(AA_pkt_uniq)
            
            all_columns['ligand_ID'].append(fold)
            all_columns['Protein_PDB_ID'].append((protein[indx])[j])
            all_columns['Pocket'].append((pkt_num[indx])[j])
            all_columns['AA_pkt'].append(AA_pkt_uniq[7:])
            
    os.chdir(path)
    with open(results_file, "a+") as outfiles:
        writer = csv.writer(outfiles, delimiter = "\t")        
        writer.writerows(zip(*[all_columns[key] for key in columns]))
outfiles.close()



