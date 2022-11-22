#! /bin/bash
# convert ligand SMILES into 3D strucutres with .pdbqt extension
for f in ./*.smi; do
	b=`basename $f .smi`
	echo Processing ligand $b
	obabel -i smi $f -o pdbqt -O $b.pdbqt -m -xh --gen3d --partialcharge gasteiger
done
