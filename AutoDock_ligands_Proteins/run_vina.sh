#! /bin/bash
# do not delete the line above

# this code performs docking between one ligand and multiple proteins
declare -a list=("0HK" "7LD" "7MA" "8NU" "40F" "89F" "ADN" "GGL" "GLU" "SRO" "Z99")

declare -a num_prots
num_prots[0]="0HK_P08173_A_5dsg.pdbqt;0HK_P11229_A_5cxv.pdbqt"
num_prots[1]="7LD_P28223_A_6wgt.pdbqt;7LD_P41595_A_5tvn.pdbqt"
num_prots[2]="7MA_O43613_A_6tod.pdbqt;7MA_O43614_A_5wqc.pdbqt"
num_prots[3]="8NU_P14416_A_6cm4.pdbqt;8NU_P28223_A_6a93.pdbqt"
num_prots[4]="40F_Q14416_C_4xaq.pdbqt;40F_Q14832_C_4xar.pdbqt"
num_prots[5]="89F_P28222_A_5v54.pdbqt;89F_P28223_A_6wh4.pdbqt"
num_prots[6]="ADN_P29274_A_2ydo.pdbqt;ADN_P30542_A_6d9h_R.pdbqt"
num_prots[7]="GGL_O00222_C_6bsz.pdbqt;GGL_Q14416_C_5cni.pdbqt;GGL_Q14832_C_5cnk.pdbqt"
num_prots[8]="GLU_A0A173M0G2_TR2_5x2p_A.pdbqt;GLU_E9P5T5_NOT_4io2_A.pdbqt;GLU_P41594_C_3lmk_A.pdbqt;GLU_P42264_NOT_3s9e_A.pdbqt;GLU_Q14416_C_7mtr_A.pdbqt"
num_prots[9]="SRO_P08908_A_7e2y_R.pdbqt;SRO_P37231_NOT_3adv_B.pdbqt"
num_prots[10]="Z99_P41594_C_7fd9_A.pdbqt;Z99_Q13255_C_3ks9_A.pdbqt;Z99_Q14831_C_3mq4.pdbqt;Z99_Q14832_C_7wi6_A.pdbqt"

i=0
for dir in "${list[@]}"; do
    cd ./$dir
    IFS=";" read -r -a arr <<< "${num_prots[$i]}"
    for prot in "${arr[@]}"; do  # in the ligand folder of ligands that bind to the protein
        b=`basename $prot .pdbqt`      # separate the name from the file extension
        echo Processing Protein $b and $dir  # print which protein is being worked on

        vina --config ${b}_config.txt --ligand ../Control_Data_Ligands/${dir}.pdbqt --out ${b}_out.pdbqt --log ${b}_log.txt

    done
    ((i=$i+1))
    cd ../
done
