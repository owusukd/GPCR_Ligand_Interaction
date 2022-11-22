#! /bin/bash
# Run APoc

declare -a list=("0HK" "7LD" "7MA" "8NU" "40F" "89F" "ADN" "GGL" "GLU" "SRO" "Z99")
declare -a indx
declare -a num_prots
num_prots[0]="0HK_P08173_A_5dsg.pdb;0HK_P11229_A_5cxv.pdb"
num_prots[1]="7LD_P28223_A_6wgt.pdb;7LD_P41595_A_5tvn.pdb"
num_prots[2]="7MA_O43613_A_6tod.pdb;7MA_O43614_A_5wqc.pdb"
num_prots[3]="8NU_P14416_A_6cm4.pdb;8NU_P28223_A_6a93.pdb"
num_prots[4]="40F_Q14416_C_4xaq.pdb;40F_Q14832_C_4xar.pdb"
num_prots[5]="89F_P28222_A_5v54.pdb;89F_P28223_A_6wh4.pdb"
num_prots[6]="ADN_P29274_A_2ydo.pdb;ADN_P30542_A_6d9h_R.pdb"
num_prots[7]="GGL_O00222_C_6bsz.pdb;GGL_Q14416_C_5cni.pdb;GGL_Q14832_C_5cnk.pdb"
num_prots[8]="GLU_A0A173M0G2_TR2_5x2p_A.pdb;GLU_E9P5T5_NOT_4io2_A.pdb;GLU_P41594_C_3lmk_A.pdb;GLU_P42264_NOT_3s9e_A.pdb;GLU_Q14416_C_7mtr_A.pdb"
num_prots[9]="SRO_P08908_A_7e2y_R.pdb;SRO_P37231_NOT_3adv_B.pdb"
num_prots[10]="Z99_P41594_C_7fd9_A.pdb;Z99_Q13255_C_3ks9_A.pdb;Z99_Q14831_C_3mq4.pdb;Z99_Q14832_C_7wi6_A.pdb"

declare -i i=0
for dir in "${list[@]}"; do
    cd ./$dir
    
        IFS=";" read -r -a arr <<< "${num_prots[$i]}"
        
        #get the length of the array arr
        len=`expr ${#arr[@]} - 1`
        
        #create all possible pairs of the indices of the protein and perform APoc
        set -- `seq 0 $len`
        for a; do
                shift
                for b; do
                        /Users/kwabena/apoc/bin/apoc ${arr[$a]} ${arr[$b]} -pvol 50 -plen 5 > ${arr[$a]}_vs_${arr[$b]}_pocket_compare_results.txt
                done
        done
    i=`expr $i + 1`
    cd ../
done
