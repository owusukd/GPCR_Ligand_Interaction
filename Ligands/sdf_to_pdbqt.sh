#! /bin/bash
# convert files from .mol2 to .pdbqt
for f in ./*.sdf; do
	b=`basename $f .sdf`
	echo Processing ligand $b
	obabel -i sdf $f -o pdbqt -O $b.pdbqt -xh --partialcharge gasteiger
done
