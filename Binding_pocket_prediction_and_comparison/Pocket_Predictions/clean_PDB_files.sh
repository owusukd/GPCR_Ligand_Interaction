#! /bin/bash
## this script removes unwanted lines in the pdb files to put them in the right
## input format for P2Rank

declare -a list=("0HK" "7LD" "7MA" "8NU" "40F" "89F" "97V" "ADN" "GGL" "GLU" "SRO" "Z99")

for dir in "${list[@]}"; do
    cd ./$dir
    
    ## convert pdb files to txt files
    for f in ./*.pdb; do
        mv "$f" "${f%.*}.txt"
    done

    for file in ./*.txt; do
        echo $file
        b=`basename $file .txt`
        grep ATOM $file > ${b}.pdb
        obabel -i pdb ${b}.pdb -o pdbqt -O ${b}.pdbqt -xh --partialcharge gasteiger
    done
    rm *.txt
    
    cd ../
done

