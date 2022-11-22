###
# This python file is used to parese APoc bindiing pockets comparison results text
# files and write them to a csv file
# This is done for all ligands at a time
###

from itertools import combinations
import csv
# set directory 
import os
path = '/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control/Binding_pocket_prediction_and_comparison/APoc'
os.chdir(path)

protein = [["0HK_P08173_A_5dsg.pdb","0HK_P11229_A_5cxv.pdb"], ["7LD_P28223_A_6wgt.pdb","7LD_P41595_A_5tvn.pdb"], ["7MA_O43613_A_6tod.pdb","7MA_O43614_A_5wqc.pdb"],["8NU_P14416_A_6cm4.pdb","8NU_P28223_A_6a93.pdb"],["40F_Q14416_C_4xaq.pdb","40F_Q14832_C_4xar.pdb"],["89F_P28222_A_5v54.pdb","89F_P28223_A_6wh4.pdb"],["ADN_P29274_A_2ydo.pdb","ADN_P30542_A_6d9h_R.pdb"],["GGL_O00222_C_6bsz.pdb","GGL_Q14416_C_5cni.pdb","GGL_Q14832_C_5cnk.pdb"],["GLU_A0A173M0G2_TR2_5x2p_A.pdb","GLU_E9P5T5_NOT_4io2_A.pdb","GLU_P41594_C_3lmk_A.pdb","GLU_P42264_NOT_3s9e_A.pdb","GLU_Q14416_C_7mtr_A.pdb"],["SRO_P08908_A_7e2y_R.pdb","SRO_P37231_NOT_3adv_B.pdb"],["Z99_P41594_C_7fd9_A.pdb","Z99_Q13255_C_3ks9_A.pdb","Z99_Q14831_C_3mq4.pdb","Z99_Q14832_C_7wi6_A.pdb"]]

folder = ["0HK", "7LD", "7MA", "8NU", "40F", "89F", "ADN", "GGL", "GLU", "SRO", "Z99"]

columns = ['Ligand_ID','pocket_1', 'pocket_2', 'PS_score', 'P_value']
out_results_file_name = 'Combine_pocket_comp_results.tsv'

# create a csv file with only column heads
with open(out_results_file_name, "w") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    writer.writerow(columns)
    
    all_columns = {}
    for i in range(len(columns)):
        all_columns.setdefault(columns[i],[])
        
    for indx, fold in enumerate(folder):
        path_list = [path, fold]
        paths = "/".join(path_list)
        os.chdir(paths)
        
        length = len(protein[indx])
        combs = combinations(range(0,length), 2)
        for j in list(combs):
            # Get the file names and read them j = (0, 1)
            results_file_name = "_".join([(protein[indx])[j[0]], 'vs', (protein[indx])[j[1]], 'pocket_compare_results.txt'])
            lines = []
            with open(results_file_name, "r") as Results:
                for line in Results:
                    if not line.isspace():
                        line = line.strip()
                        lines.append(line)
            Results.close()
                
            # parse the text file and pick info needed
            for num, line in enumerate(lines,1):
                column = line.split()
                if column[0] == ">>>>>>>>>>>>>>>>>>>>>>>>>":
                    if column[1] == "Pocket":
                        if lines[num] == "The number of values to be sorted is not positive.":
                            continue
                        all_columns['Ligand_ID'].append(fold)
                        all_columns['pocket_1'].append(lines[num].split()[7].split(':')[1])
                        all_columns['pocket_2'].append(lines[1+num].split()[7].split(':')[1])
                        vals_1 = lines[2+num].split(',')
                        all_columns['PS_score'].append(vals_1[0].split('=')[1])
                        all_columns['P_value'].append(vals_1[1].split('=')[1])
            
    #append the results to the tsv file
    writer.writerows(zip(*[all_columns[key] for key in columns]))
outfile.close()
