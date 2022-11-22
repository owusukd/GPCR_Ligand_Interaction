###
# This python file is used to align one set of ligands with another set
# in a pairwise manner
# The RMSD from the alignment are written to a csv file
# This is done for a single ligand at a time
###

import pymol
from pymol import cmd
# set directory
import os
os.chdir('/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control/Ligand_Pose_Comformation_Docked')
import csv
from itertools import combinations

ligs = ["0HK","7LD","7MA","8NU","40F","89F","ADN","GGL","GLU","SRO","Z99"]
protein = [["P08173_A_5dsg","P11229_A_5cxv"], ["P28223_A_6wgt","P41595_A_5tvn"], ["O43613_A_6tod","O43614_A_5wqc"],["P14416_A_6cm4","P28223_A_6a93"],["Q14416_C_4xaq","Q14832_C_4xar"],["P28222_A_5v54","P28223_A_6wh4"],["P29274_A_2ydo","P30542_A_6d9h_R"], ["O00222_C_6bsz","Q14416_C_5cni","Q14832_C_5cnk"],["A0A173M0G2_TR2_5x2p_A","E9P5T5_NOT_4io2_A","P41594_C_3lmk_A","P42264_NOT_3s9e_A","Q14416_C_7mtr_A"],["P08908_A_7e2y_R","P37231_NOT_3adv_B"],["P41594_C_7fd9_A","Q13255_C_3ks9_A","Q14831_C_3mq4","Q14832_C_7wi6_A"]]

pkt_counts = [[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1,1],[1,1,1,1,1],[1,1],[1,1,1,1]]
columns = ['Ligand_ID','Compared_Pocket_1','Compared_Pocket_2','Pkt_Ligs_Aligned_RMSD']

results_file = 'Ligs_Align_Pkt_Docked.tsv'
with open(results_file, "w") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(columns)
# loop over each folder
for indx in range(0,len(pkt_counts)):
    
    # create columns for the csv file in a dictionary
    all_columns = {}
    for i in range(len(columns)):
        all_columns.setdefault(columns[i],[])
    #
    pkt_comb = combinations(enumerate(pkt_counts[indx]), 2)
    # take each combination of the number of pockets
    for j in list(pkt_comb):
        
        lig_name1 = "_".join([ligs[indx], (protein[indx])[((j)[0])[0]], "out"])
        lig_pkt_name1 = ".".join([lig_name1, 'pdbqt'])
        
        lig_name2 = "_".join([ligs[indx], (protein[indx])[((j)[1])[0]], "out"])
        lig_pkt_name2 = ".".join([lig_name2, 'pdbqt'])
        
        # PyMol begins
        cmd.load('%s'%(lig_pkt_name1))
        cmd.load('%s'%(lig_pkt_name2))
        object_list = cmd.get_names()
        
        rmsd = cmd.align('%s'%(object_list[0]), '%s'%(object_list[1]), cycles=0, object=None, target_state=0,mobile_state=0)
        
        # delete loaded ligands before start of the next comparison
        cmd.delete('all')
        # PyMol ends
        
        compared_pkt_1 = (protein[indx])[((j)[0])[0]]
        compared_pkt_2 = (protein[indx])[((j)[1])[0]]
        all_columns['Ligand_ID'].append(ligs[indx])
        all_columns['Compared_Pocket_1'].append(compared_pkt_1)
        all_columns['Compared_Pocket_2'].append(compared_pkt_2)
        all_columns['Pkt_Ligs_Aligned_RMSD'].append(rmsd[0])
        
    with open(results_file, "a+") as outfiles:
        writer = csv.writer(outfiles, delimiter = "\t")
        writer.writerows(zip(*[all_columns[key] for key in columns]))
outfiles.close()
