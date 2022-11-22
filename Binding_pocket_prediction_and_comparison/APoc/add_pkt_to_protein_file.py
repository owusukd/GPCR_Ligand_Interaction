###
# This python file is used to add the pocket file to the GPCR 3D structure file
###

# set directory
import os
path = '/Users/kwabena/Research/GPCR/Manuscript/MDPI biomolecules/New_Data_Control/Binding_pocket_prediction_and_comparison/APoc'
os.chdir(path)

protein = [["0HK_P08173_A_5dsg.pdb","0HK_P11229_A_5cxv.pdb"], ["7LD_P28223_A_6wgt.pdb","7LD_P41595_A_5tvn.pdb"], ["7MA_O43613_A_6tod.pdb","7MA_O43614_A_5wqc.pdb"],["8NU_P14416_A_6cm4.pdb","8NU_P28223_A_6a93.pdb"],["40F_Q14416_C_4xaq.pdb","40F_Q14832_C_4xar.pdb"],["89F_P28222_A_5v54.pdb","89F_P28223_A_6wh4.pdb"],["ADN_P29274_A_2ydo.pdb","ADN_P30542_A_6d9h_R.pdb"],["GGL_O00222_C_6bsz.pdb","GGL_Q14416_C_5cni.pdb","GGL_Q14832_C_5cnk.pdb"],["GLU_A0A173M0G2_TR2_5x2p_A.pdb","GLU_E9P5T5_NOT_4io2_A.pdb","GLU_P41594_C_3lmk_A.pdb","GLU_P42264_NOT_3s9e_A.pdb","GLU_Q14416_C_7mtr_A.pdb"],["SRO_P08908_A_7e2y_R.pdb","SRO_P37231_NOT_3adv_B.pdb"],["Z99_P41594_C_7fd9_A.pdb","Z99_Q13255_C_3ks9_A.pdb","Z99_Q14831_C_3mq4.pdb","Z99_Q14832_C_7wi6_A.pdb"]]

pkt_row_num = [[0,0], [0,0], [0,0], [0,0], [0,0], [0,0], [0,0], [0,0,1], [0,0,0,0,1], [0,0], [0,0,4,9]]
folder = ["0HK", "7LD", "7MA", "8NU", "40F", "89F", "ADN", "GGL", "GLU", "SRO", "Z99"]

for ii, fold in enumerate(folder):
    path_list = [path, fold]
    paths = "/".join(path_list)
    os.chdir(paths)
    for jj, prot in enumerate(protein[ii]):
        file = open(prot, "a+")
        
        #
        pkt_num = pkt_row_num[ii][jj] + 1
        pkt_name = "_".join([prot.split(".")[0], "pkt", str(pkt_num)])
        pkt_file_name = ".".join([pkt_name, "txt"])
        pkt_file = open(pkt_file_name, "r")
        pkt = pkt_file.read()
        file.write(pkt)
        
        pkt_file.close()
        file.close()
        
