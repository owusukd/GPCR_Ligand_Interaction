#! /bin/bash
# do not delete the line above

declare -a list=("0HK" "7LD" "7MA" "8NU" "40F" "89F" "ADN" "GGL" "GLU" "SRO" "Z99")

i=0
for dir in "${list[@]}"; do
    cd ./$dir
    chmod a+w *
    chmod a+r *
    cd ../
done
