#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 20:52:19 2022
This code scores the overlap between the pockets of protein A and the part
of the protein A which is similar to protein B. Same is done for protein B.
This process is done both for flexible comaprison and rigid comparison of 
protein A and protein B.
The code loop over each of the ligand folders, and the proteins in them and 
the number of pockets for each protein.
Finally, it write all the results as one file.
@author: Kwabena Owusu Dankwah
"""

from itertools import permutations
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

columns = ['ligand_ID', 'Proteins_Compared', 'Pocket', 'P3D_similar_Coincide_pkt_flex', 'P3D_similar_Coincide_pkt_rigid']

def intersection_len(list1, list2):
    list3 = [value for value in list1 if value in list2]
    return len(list3)

def score_overlap(position_pocket, position_p3d_similar):
    len_pos_pkt = len(position_pocket)
    len_intersect = intersection_len(position_pocket, position_p3d_similar)
    score = len_intersect/len_pos_pkt
    return score
    
results_file = '3D_Similar_Overlap_Pkt.tsv'
# create a csv file with only column heads
with open(results_file, "w") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(columns)

# loop over the ligand folders
for indx, fold in enumerate(folder):
    path_list = [path, fold]
    paths = "/".join(path_list)
    os.chdir(paths)
    
    # create dict to hold the results
    all_columns = {}
    for i in range(len(columns)):
        all_columns.setdefault(columns[i],[])
    
    length = len(protein[indx])
    perm = permutations(range(0,length), 2)
    
    # loop over permutations of the pairs of protein pockets
    for j in list(perm):
        # Get the file names and read them
        proteins_comped = "_".join([(protein[indx])[j[0]], 'with', (protein[indx])[j[1]]])
        
        flex = "_".join([(protein[indx])[j[0]], 'flex_with', (protein[indx])[j[1]]])
        rigid = "_".join([(protein[indx])[j[0]], 'rigid_with', (protein[indx])[j[1]]])
        flex_file_name = ".".join([flex, "txt"])
        rigid_file_name = ".".join([rigid, "txt"])
        
        flex_file = open(flex_file_name)
        d_similar_flex = flex_file.read()
        
        rigid_file = open(rigid_file_name) 
        d_similar_rigid = rigid_file.read()
        
        # get the unique positions
        regex = r"[\sA-Z]\d{1,}(?=\s{4})"
        position_flex = re.finditer(regex, d_similar_flex)
        position_rigid = re.finditer(regex, d_similar_rigid)
        
        position_flex_uniq = []
        for pos_flex in position_flex:
            pos_fl = pos_flex.group()
            pos_fl = int(re.sub(r"\D","", pos_fl))
            if pos_fl not in position_flex_uniq:
                position_flex_uniq.append(pos_fl)
        
        position_rigid_uniq = []
        for pos_rigid in position_rigid:
            pos_ri = pos_rigid.group()
            pos_ri = int(re.sub(r"\D","", pos_ri))
            if pos_ri not in position_rigid_uniq:
                position_rigid_uniq.append(pos_ri)
        
        # compare pockets to the portions of the proteins that are 3D structurally similar
        pkt_name = "_".join([(protein[indx])[j[0]], 'pkt', str((pkt_num[indx])[j[0]])])
        pkt_file_name = ".".join([pkt_name, 'txt'])
        pkt_file = open(pkt_file_name)
        pkt = pkt_file.read()
        
        position_pkt = re.finditer(regex, pkt)
        position_pkt_uniq = []
        for pos_pkt in list(position_pkt)[2:]:
            pos_pk = pos_pkt.group()
            pos_pk = int(re.sub(r"\D","", pos_pk))
            if pos_pk not in position_pkt_uniq:
                position_pkt_uniq.append(pos_pk)
        
        # save results to dict created above
        all_columns['ligand_ID'].append(fold)
        all_columns['Proteins_Compared'].append(proteins_comped)
        all_columns['Pocket'].append(pkt_name)
        all_columns['P3D_similar_Coincide_pkt_flex'].append(score_overlap(position_pkt_uniq, position_flex_uniq))
        all_columns['P3D_similar_Coincide_pkt_rigid'].append(score_overlap(position_pkt_uniq, position_rigid_uniq))

    # change directory and write the results
    os.chdir(path)
    with open(results_file, "a+") as outfiles:
        writer = csv.writer(outfiles, delimiter = "\t")        
        writer.writerows(zip(*[all_columns[key] for key in columns]))
outfiles.close()






