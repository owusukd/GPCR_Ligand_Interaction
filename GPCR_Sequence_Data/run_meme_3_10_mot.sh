#!/bin/bash

for dir in ./fasta_MEME_3_10_mot/*/
do
	for file in $dir*.txt
	do
		meme $file -protein -oc ${file%.*} -nostatus -time 18000 -mod zoops -nmotifs 10 -minw 3 -maxw 20
	done 
done
