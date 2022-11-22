###
# This code create a PDB file of the pockets with .txt extension
# it also create configuration files for AutoDock vina
###

# set directory
import os
path = '/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control/Binding_pocket_prediction_and_comparison/Pocket_Predictions'
os.chdir(path)
import csv

protein = [["0HK_P08173_A_5dsg.pdb","0HK_P11229_A_5cxv.pdb"], ["7LD_P28223_A_6wgt.pdb","7LD_P41595_A_5tvn.pdb"], ["7MA_O43613_A_6tod.pdb","7MA_O43614_A_5wqc.pdb"],["8NU_P14416_A_6cm4.pdb","8NU_P28223_A_6a93.pdb"],["40F_Q14416_C_4xaq.pdb","40F_Q14832_C_4xar.pdb"],["89F_P28222_A_5v54.pdb","89F_P28223_A_6wh4.pdb"],["ADN_P29274_A_2ydo.pdb","ADN_P30542_A_6d9h_R.pdb"],["GGL_O00222_C_6bsz.pdb","GGL_Q14416_C_5cni.pdb","GGL_Q14832_C_5cnk.pdb"],["GLU_A0A173M0G2_TR2_5x2p_A.pdb","GLU_E9P5T5_NOT_4io2_A.pdb","GLU_P41594_C_3lmk_A.pdb","GLU_P42264_NOT_3s9e_A.pdb","GLU_Q14416_C_7mtr_A.pdb"],["SRO_P08908_A_7e2y_R.pdb","SRO_P37231_NOT_3adv_B.pdb"],["Z99_P41594_C_7fd9_A.pdb","Z99_Q13255_C_3ks9_A.pdb","Z99_Q14831_C_3mq4.pdb","Z99_Q14832_C_7wi6_A.pdb"]]

pkt_row_num = [[0,0], [0,0], [0,0], [0,0], [0,0], [0,0], [0,0], [0,0,1], [0,0,0,0,1], [0,0], [0,0,4,9]]
folder = ["0HK", "7LD", "7MA", "8NU", "40F", "89F", "ADN", "GGL", "GLU", "SRO", "Z99"]

for ii, fold in enumerate(folder):
    path_list = [path, fold]
    paths = "/".join(path_list)
    os.chdir(paths)
    for jj, prot in enumerate(protein[ii]):
        file_name = "_".join([prot, "predictions.csv"])
        file = open(file_name)
        csv_file = file.readlines()[1:]
        pkt_info = csv_file[pkt_row_num[ii][jj]].split(",")
        
        #write docking configuration file
        center_x = pkt_info[6].strip()
        center_y = pkt_info[7].strip()
        center_z = pkt_info[8].strip()
        receptor_name = "".join([prot, "qt"])
        config_name = "_".join([prot.split(".")[0], "config.txt"])
        #write docking configuration file content and close file
        config_file = open(config_name, "w")
        config_file.write("receptor = %s\n\n" %receptor_name)
        config_file.write("center_x = %s\n" %center_x)
        config_file.write("center_y = %s\n" %center_y)
        config_file.write("center_z = %s\n\n" %center_z)
        config_file.write("size_x = 20\n")
        config_file.write("size_y = 20\n")
        config_file.write("size_z = 20\n\n")
        config_file.write("num_modes = 1\n")
        config_file.write("exhaustiveness = 9")
        config_file.close()
        
        #get pockets AA and write a txt file of the atoms coordinates
        split_val = pkt_info[9].strip()[0:2]  # chain name eg A, B, C with _ eg A_, B_, C_
        AA_ind = pkt_info[9].split(split_val)[1:] # get the positions of the amino acids in pocket
        AA_ind = list(map(int, AA_ind))
        AA_ind.sort()
        AA_indx = list(map(str, AA_ind))
        prot_file_pdb = open(prot, "r")
        prot_file = prot_file_pdb.readlines()
        #create pocket file name
        pkt_num = pkt_row_num[ii][jj] + 1
        pkt_name = "_".join([prot.split(".")[0], "pkt", str(pkt_num)])
        pkt_file_name = ".".join([pkt_name, "txt"])
        #write pocket file
        pkt_file = open(pkt_file_name, "a+")
        pkt_file.write("PKT        11    101    " + pkt_name + "\n") # Important for APoc
        for indx in AA_indx:
            for lines in prot_file[:-1]:
                line = lines.split(" ")
                lin = ' '.join(line).split()
                if lin[0] == "ATOM":
                    lin_5 = lin[5].split(split_val[0])
                    lin_5 = ' '.join(lin_5).split()[0]
                    if lin_5 == indx:
                        pkt_file.write(lines)
        pkt_file.write("TER") # Important for PDB file
        pkt_file.close()
        prot_file_pdb.close()
        
