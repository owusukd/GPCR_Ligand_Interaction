#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:51:15 2022
This code maps the 3D structural alignment of two proteins to their respective pdb files
and extract the portions of the 3D structures that were found to be similar.
The extracted portions are written to a pdb file.
This is done both for flexible and rigid 3D structural comparisons.
@author: kwabena
"""

import re
import csv
from itertools import combinations
from Bio import SeqIO

# set directory
import os
path = '/Users/kwabena/Research/GPCR/Entire_work_organized/3D_Comparison_Pocket_overlap_scores/3D_structural_comparison'

protein = [["0HK_P08173_A_5dsg","0HK_P11229_A_5cxv"], ["7LD_P28223_A_6wgt","7LD_P41595_A_5tvn"], 
           ["7MA_O43613_A_6tod","7MA_O43614_A_5wqc"], ["8NU_P14416_A_6cm4","8NU_P28223_A_6a93"], 
           ["40F_Q14416_C_4xaq","40F_Q14832_C_4xar"], ["89F_P28222_A_5v54","89F_P28223_A_6wh4"], 
           ["ADN_P29274_A_2ydo","ADN_P30542_A_6d9h_R"],["GGL_O00222_C_6bsz","GGL_Q14416_C_5cni","GGL_Q14832_C_5cnk"], 
           ["GLU_A0A173M0G2_TR2_5x2p_A","GLU_E9P5T5_NOT_4io2_A","GLU_P41594_C_3lmk_A","GLU_P42264_NOT_3s9e_A","GLU_Q14416_C_7mtr_A"], 
           ["SRO_P08908_A_7e2y_R","SRO_P37231_NOT_3adv_B"], 
           ["Z99_P41594_C_7fd9_A","Z99_Q13255_C_3ks9_A","Z99_Q14831_C_3mq4","Z99_Q14832_C_7wi6_A"]]

folder = ["0HK", "7LD", "7MA", "8NU", "40F", "89F", "ADN", "GGL", "GLU", "SRO", "Z99"]

columns = ['Ligand','Protein_Pair_1', 'Protein_Pair_2', '3D_RMSD_flex', '3D_RMSD_rigid', 'Aligned_Seq_Identity_flex', 
           'Aligned_Seq_Identity_rigid', 'Aligned_Seq_Similarity_flex', 'Aligned_Seq_Similarity_rigid']
results_file = '3D_Similar_RMSD.tsv'

# create a csv file with only column heads
os.chdir(path)
with open(results_file, "w") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(columns) 
#outfile.close()

for indx, fold in enumerate(folder):
    path_list = [path, fold]
    paths = "/".join(path_list)
    os.chdir(paths)
    
    # create dict to hold the results
    all_columns = {}
    for i in range(len(columns)):
        all_columns.setdefault(columns[i],[])
    
    length = len(protein[indx])
    comb = combinations(range(0,length), 2)
    
    # loop over combinations of the pairs of protein pockets
    for j in comb:
        # Get the alignment file names and read them. Also read protein files
        
        flex_aln_name = "_".join([(protein[indx])[j[0]], (protein[indx])[j[1]], 'flex.aln'])
        rigid_aln_name = "_".join([(protein[indx])[j[0]], (protein[indx])[j[1]], 'rigid.aln'])
        
        flex_aln = open(flex_aln_name, "r")
        aln_similar_flex = flex_aln.readlines()
        
        rigid_aln = open(rigid_aln_name, "r") 
        aln_similar_rigid = rigid_aln.readlines()
        
        prot_1 = [record.seq for record in SeqIO.parse(".".join([(protein[indx])[j[0]], 'pdb']), "pdb-atom")]
        prot_1_seq = str(prot_1[0]).replace("X","")
        prot_2 = [record.seq for record in SeqIO.parse(".".join([(protein[indx])[j[1]], 'pdb']), "pdb-atom")]
        prot_2_seq = str(prot_2[0]).replace("X","") 
        
        prot_11 = open(".".join([(protein[indx])[j[0]], 'pdb']), "r")
        prot_11_pdb = prot_11.read()
        prot_22 = open(".".join([(protein[indx])[j[1]], 'pdb']), "r")
        prot_22_pdb = prot_22.read()
        
        
        # get the positions of the sequences in the pdb files
        regex = r"[\sA-Z]\d{1,}(?=\s{4})"
        position_prot_1 = re.finditer(regex, prot_11_pdb)
        position_prot_2 = re.finditer(regex, prot_22_pdb)
        
        position_prot_1_uniq = []
        for pos_prot_1 in position_prot_1:
            pos_pt_1 = pos_prot_1.group()
            pos_pt_1 = int(re.sub(r"\D","", pos_pt_1))
            if pos_pt_1 not in position_prot_1_uniq:
                position_prot_1_uniq.append(pos_pt_1)
        
        position_prot_2_uniq = []
        for pos_prot_2 in position_prot_2:
            pos_pt_2 = pos_prot_2.group()
            pos_pt_2 = int(re.sub(r"\D","", pos_pt_2))
            if pos_pt_2 not in position_prot_2_uniq:
                position_prot_2_uniq.append(pos_pt_2)
        
        # map the alignment to the pdb files and subset the pdb for the blocks of the alignment
        position_chain_1_flex = []  
        position_chain_2_flex = [] 
        position_chain_1_rigid = []  
        position_chain_2_rigid = [] 
        
        prot_file_ = open(".".join([(protein[indx])[j[0]], 'pdb']), "r")
        prot_file_1 = prot_file_.readlines()
        prot_file__ = open(".".join([(protein[indx])[j[1]], 'pdb']), "r")
        prot_file_2 = prot_file__.readlines()
        
        # Flex
        temp_flex_1 = position_prot_1_uniq[:]
        temp_flex_2 = position_prot_2_uniq[:]
        for ii, line_aln_flex in enumerate(aln_similar_flex):
            
            if line_aln_flex[:7] == 'Chain 1':
                line_aln_flex = line_aln_flex.split(" ")
                ind_1 = position_prot_1_uniq.index(int(line_aln_flex[-2])) #the index of the 1st AA in aln in position_prot_1_uniq
                indx__1 = [(i-14+ind_1) for i, ij in enumerate(line_aln_flex[-1].replace("\n", "")) if ij=='-']
                if len(indx__1) > 0:  
                    for ind in indx__1:
                            temp_flex_1.insert(ind, "-")
                
                block_indx_1 = [(i-14+ind_1) for i, ij in enumerate(aln_similar_flex[(ii+1)].replace("\n", "")) if ij!=' ']
                if len(block_indx_1) > 0:  
                    for jj in block_indx_1:
                        position_chain_1_flex.append(temp_flex_1[jj])
                   
            if line_aln_flex[:7] == 'Chain 2':
                line_aln_flex = line_aln_flex.split(" ")
                ind_2 = position_prot_2_uniq.index(int(line_aln_flex[-2])) #the index of the 1st AA in aln in position_prot_2_uniq
                indx__2 = [(ik-14+ind_2) for ik, ji in enumerate(line_aln_flex[-1].replace("\n", "")) if ji=='-']
                if len(indx__2) > 0:  
                    for ind in indx__2:
                        temp_flex_2.insert(ind, "-")
                
                block_indx_2 = [(ik-14+ind_2) for ik,ji in enumerate(aln_similar_flex[(ii-1)].replace("\n", "")) if ji!=' ']
                if len(block_indx_2) > 0:
                    for ki in block_indx_2:
                        position_chain_2_flex.append(temp_flex_2[ki])
        
        #create file names
        flex_1 = "_".join([(protein[indx])[j[0]], 'flex_with', (protein[indx])[j[1]]])
        flex_file_name_1 = ".".join([flex_1, "txt"])
        
        flex_2 = "_".join([(protein[indx])[j[1]], 'flex_with', (protein[indx])[j[0]]])
        flex_file_name_2 = ".".join([flex_2, "txt"])
        
        #write 3D file of the alignment
        flex_file_1 = open(flex_file_name_1, "a+")  #protein 1 flex
        for indxx in position_chain_1_flex:
            if indxx != "-":
                for lines in prot_file_1[:-1]:
                    line = lines.split(" ")
                    lin = ' '.join(line).split()
                    len_lin = len(lin)
                    if lin[0] == "ATOM":
                        if len_lin == 12:
                            if lin[5] == str(indxx):
                                flex_file_1.write(lines)
                        if len_lin == 11:
                            lin_5 = re.sub(r"\D","", lin[4])
                            if lin_5 == str(indxx):
                                flex_file_1.write(lines)
                        
        flex_file_1.write("TER") # Important for PDB file
        flex_file_1.close()
        
        flex_file_2 = open(flex_file_name_2, "a+")    #protein 2 flex
        for indxx in position_chain_2_flex:
            if indxx != "-":
                for lines in prot_file_2[:-1]:
                    line = lines.strip().split(" ")
                    lin = ' '.join(line).split()
                    len_lin = len(lin)
                    if lin[0] == "ATOM":
                        if len_lin == 12:
                            if lin[5] == str(indxx):
                                flex_file_2.write(lines)
                        elif len_lin == 11:
                            lin_5 = re.sub(r"\D","", lin[4])
                            if lin_5 == str(indxx):
                                flex_file_2.write(lines)
                        
        flex_file_2.write("TER") # Important for PDB file
        flex_file_2.close()
        
        # Rigid 
        temp_rigid_1 = position_prot_1_uniq[:]
        temp_rigid_2 = position_prot_2_uniq[:]
        for iii, line_aln_rigid in enumerate(aln_similar_rigid):
            
            if line_aln_rigid[:7] == 'Chain 1':
                line_aln_rigid = line_aln_rigid.split(" ")
                ind_1 = position_prot_1_uniq.index(int(line_aln_rigid[-2])) #the index of the 1st AA in aln in position_prot_1_uniq
                indx__1 = [(i-14+ind_1) for i, ij in enumerate(line_aln_rigid[-1].replace("\n", "")) if ij=='-']
                if len(indx__1) > 0:
                    for ind in indx__1:
                        temp_rigid_1.insert(ind, "-")
                
                block_indx_1 = [(jk-14+ind_1) for jk ,kj in enumerate(aln_similar_rigid[(iii+1)].replace("\n", "")) if kj!=' ']
                if len(block_indx_1) > 0: 
                    for jjj in block_indx_1:
                        position_chain_1_rigid.append(temp_rigid_1[jjj])
                  
            if line_aln_rigid[:7] == 'Chain 2':
                line_aln_rigid = line_aln_rigid.split(" ")
                ind_2 = position_prot_2_uniq.index(int(line_aln_rigid[-2])) #the index of the 1st AA in aln in position_prot_2_uniq
                indx__2 = [(i-14+ind_2) for i, ij in enumerate(line_aln_rigid[-1].replace("\n", "")) if ij=='-']
                if len(indx__2) > 0:
                    for ind in indx__2:
                        temp_rigid_2.insert(ind, "-")
                
                block_indx_2 = [(kk-14+ind_2) for kk,ijk in enumerate(aln_similar_rigid[(iii-1)].replace("\n", "")) if ijk!=' ']
                if len(block_indx_2) > 0:
                    for kkk in block_indx_2:
                        position_chain_2_rigid.append(temp_rigid_2[kkk])    

        rigid_1 = "_".join([(protein[indx])[j[0]], 'rigid_with', (protein[indx])[j[1]]])
        rigid_file_name_1 = ".".join([rigid_1, "txt"])
        
        rigid_2 = "_".join([(protein[indx])[j[1]], 'rigid_with', (protein[indx])[j[0]]])
        rigid_file_name_2 = ".".join([rigid_2, "txt"])
      
        #write 3D file of the alignment
        rigid_file_1 = open(rigid_file_name_1, "a+")  #protein 1 rigid
        for indxx in position_chain_1_rigid:
            if indxx != "-":
                for lines in prot_file_1[:-1]:
                    line = lines.strip().split(" ")
                    lin = ' '.join(line).split()
                    len_lin = len(lin)
                    if lin[0] == "ATOM":
                        if len_lin == 12:
                            if lin[5] == str(indxx):
                                rigid_file_1.write(lines)
                        elif len_lin == 11:
                            lin_5 = re.sub(r"\D","", lin[4])
                            if lin_5 == str(indxx):
                                rigid_file_1.write(lines)
                        
        rigid_file_1.write("TER") # Important for PDB file
        rigid_file_1.close()
        
        rigid_file_2 = open(rigid_file_name_2, "a+")    #protein 2 rigid
        for indxx in position_chain_2_flex:
            if indxx != "-":
                for lines in prot_file_2[:-1]:
                    line = lines.strip().split(" ")
                    lin = ' '.join(line).split()
                    len_lin = len(lin)
                    if lin[0] == "ATOM":
                        if len_lin == 12:
                            if lin[5] == str(indxx):
                                rigid_file_2.write(lines)
                        elif len_lin == 11:
                            lin_5 = re.sub(r"\D","", lin[4])
                            if lin_5 == str(indxx):
                                rigid_file_2.write(lines)
                        
        rigid_file_2.write("TER") # Important for PDB file
        rigid_file_2.close()
        
        
        # save results to dict created above
        all_columns['Ligand'].append(fold)
        all_columns['Protein_Pair_1'].append(".".join([(protein[indx])[j[0]], 'pdb']))
        all_columns['Protein_Pair_2'].append(".".join([(protein[indx])[j[1]], 'pdb']))
        all_columns['3D_RMSD_flex'].append(aln_similar_flex[1].split()[9])
        all_columns['3D_RMSD_rigid'].append(aln_similar_rigid[1].split()[9])
        all_columns['Aligned_Seq_Identity_flex'].append(aln_similar_flex[2].split()[5].replace("%",""))
        all_columns['Aligned_Seq_Identity_rigid'].append(aln_similar_rigid[2].split()[5].replace("%",""))
        all_columns['Aligned_Seq_Similarity_flex'].append(aln_similar_flex[2].split()[7].replace("%","").replace("\n",""))
        all_columns['Aligned_Seq_Similarity_rigid'].append(aln_similar_rigid[2].split()[7].replace("%","").replace("\n",""))
        
    # change directory and write the RMSD results
    os.chdir(path)
    with open(results_file, "a+") as outfiles:
        writer = csv.writer(outfiles, delimiter = "\t")        
        writer.writerows(zip(*[all_columns[key] for key in columns]))

outfiles.close()
        

