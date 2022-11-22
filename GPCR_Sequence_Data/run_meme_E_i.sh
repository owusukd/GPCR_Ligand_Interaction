#!/bin/bash

for dir in ./fasta_MEME_E_i/*/
do
	for file in $dir*.txt
	do
		meme $file -protein -oc ${file%.*} -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 3 -maxw 20
	done 
done
