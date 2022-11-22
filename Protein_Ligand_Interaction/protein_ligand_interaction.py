#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 13:17:06 2022
work in progress
@author: kwabena
"""


import os

from openbabel import OBMol, OBConversion, OBMolAtomIter, OBResidueIter, OBResidueAtomIter, OBAtomAtomIter

path = '/Users/kwabena/Research/GPCR/Entire_work_organized/Pipeline/Protein_Ligand_Interaction_Actual/Test'
os.chdir(path)

def isnonpolar(atom):
    if atom.MatchesSMARTS('[#6,#16,F,Cl,Br,I]'):
        return True

    return False


def res_protein_ligand_interaction(protein_file_name, ligand_file_name):
   
    # opening the molecules files
    file_format = protein_file_name.split(".")[1]
    protein = OBMol()
    ligand = OBMol()
    conv = OBConversion()
    conv.SetInFormat(file_format)
    conv.ReadFile(protein, protein_file_name)
    conv.ReadFile(ligand, ligand_file_name)
    
    # get the residue of the protein interacting with the atom(s) of the ligand
    reslist = []
    for residue in OBResidueIter(protein):
        residuename = residue.GetName()
        for atom in OBResidueAtomIter(residue):
            if not atom.GetAtomicNum() == 1: # check if not hydrogen atom
                for atomlig in OBMolAtomIter(ligand):
                    if not atomlig.GetAtomicNum() == 1: # check if not hydrogen atom
                        distance = atom.GetDistance(atomlig)
                        if distance <= 5:
                            if isnonpolar(atomlig) & isnonpolar(atom):
                                if atom.IsHbondDonor() & atomlig.IsHbondAcceptor():
                                    for neighborDon in OBAtomAtomIter(atom):
                                        if neighborDon.IsHydrogen():
                                            angle = atom.GetAngle(neighborDon, atomlig)
                                            if angle>135.0:
                                                reslist.append(residuename)
                                if atom.IsHbondAcceptor() & atomlig.IsHbondDonor():
                                    for neighborDon in OBAtomAtomIter(atomlig):
                                        if neighborDon.IsHydrogen():
                                            angle = atomlig.GetAngle(neighborDon, atom)
                                            if angle>135.0:
                                                reslist.append(residuename)
    return reslist

# test function
res_5dsg_0HK = res_protein_ligand_interaction("0HK_P08173_A_5dsg.pdbqt", "0HK_P08173_A_5dsg_out.pdbqt")
res_5cxv_0HK = res_protein_ligand_interaction("0HK_P11229_A_5cxv.pdbqt", "0HK_P11229_A_5cxv_out.pdbqt")




