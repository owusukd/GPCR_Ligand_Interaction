###
# This python file is used to parese APoc bindiing pockets comparison results text
# files and write them to a csv file
# This is done for a single ligand at a time
###

import csv
# set directory 
import os
path = '/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control/AutoDock_ligands_Proteins'
os.chdir(path)

input_file = [["0HK_P08173_A_5dsg","0HK_P11229_A_5cxv"], ["7LD_P28223_A_6wgt","7LD_P41595_A_5tvn"], ["7MA_O43613_A_6tod","7MA_O43614_A_5wqc"],["8NU_P14416_A_6cm4","8NU_P28223_A_6a93"],["40F_Q14416_C_4xaq","40F_Q14832_C_4xar"],["89F_P28222_A_5v54","89F_P28223_A_6wh4"],["ADN_P29274_A_2ydo","ADN_P30542_A_6d9h_R"],["GGL_O00222_C_6bsz","GGL_Q14416_C_5cni","GGL_Q14832_C_5cnk"],["GLU_A0A173M0G2_TR2_5x2p_A","GLU_E9P5T5_NOT_4io2_A","GLU_P41594_C_3lmk_A","GLU_P42264_NOT_3s9e_A","GLU_Q14416_C_7mtr_A"],["SRO_P08908_A_7e2y_R","SRO_P37231_NOT_3adv_B"],["Z99_P41594_C_7fd9_A","Z99_Q13255_C_3ks9_A","Z99_Q14831_C_3mq4","Z99_Q14832_C_7wi6_A"]]

ligand = ["0HK", "7LD", "7MA", "8NU", "40F", "89F", "ADN", "GGL", "GLU", "SRO", "Z99"]

# define the APoc results text files to parse and the csv file to write to
results_file = 'AutoDock_vina_Results.tsv'
#counts = [[1,1], [1,1], [1,1], [1,1], [1,1], [1,1], [1,1], [1,1,1], [1,1,1,1,1], [1,1], [1,1,1,1]]

# function to parse results text file
columns = ['Ligand_ID', 'Protein_ID', 'Affinity_kcal_mol']


# create a csv file with only column heads
with open(results_file, "w") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(columns)
    
    all_columns = {}
    for i in range(len(columns)):
        all_columns.setdefault(columns[i],[])
        
    for ii, file in enumerate(input_file):

        paths = "/".join([path, str(ligand[ii])])
        os.chdir(paths)
        
        # read text file
        for log_file in file:
            
            file_name = '.'.join(['_'.join([log_file, 'log']),'txt'])
            
            with open(file_name, "r") as Results:
                for line in Results:
                    if not line.isspace():
                        line = line.strip()
                        column = line.split()
                        if column[0] == "1":
                            all_columns['Ligand_ID'].append(ligand[ii])
                            all_columns['Protein_ID'].append(log_file)
                            all_columns['Affinity_kcal_mol'].append(column[1])
            Results.close()
            
    # append the results to the csv file
    # zip() transpose the original data
    writer.writerows(zip(*[all_columns[key] for key in columns]))
outfile.close()
                
