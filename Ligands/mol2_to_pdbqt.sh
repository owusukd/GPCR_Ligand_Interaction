#! /bin/bash
# convert files from .mol2 to .pdbqt
for f in ./*.mol2; do
	b=`basename $f .mol2`
	echo Processing ligand $b
	obabel -i mol2 $f -o pdbqt -O $b.pdbqt -xh --partialcharge gasteiger
done
