#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 09 19:53:42 2020

@authors: Yann Vander Meersche - Debbah Nagi

This script calculate intra-proteins interactions in a pdb file given

Inputs:
=======
pdb file


Output:
========
interactions table

Usage:
======
    $ python PIC.py -p pdb_file
"""

################################################################################
import os
import sys
import argparse
import re

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from tabulate import tabulate
################################################################################


# Variable declarations
hdc_aa = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"]
ani_aa = ["ASP", "GLU"]
cat_aa = ["ARG", "HIS", "LYS"]
ion_aa = ["ASP", "GLU", "ARG", "HIS", "LYS"]

def args():
    """Parse the command-line arguments.

    Arguments
    ---------
    Python command-line

    Returns
    -------
    pdb_file: string
    """

    # Declaration of expexted arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", help="Path to the PDB file.", type=str, required=True)
    args = parser.parse_args()

    # Check if the PDB files directory is valid
    pdb_file = args.pdb
    if not os.path.isfile(pdb_file):
        sys.exit(f"{pdb_file} does not exist.\n"
                "Please enter a valid PDB file.")

    return pdb_file

def parse_pdb(pdb_file):
    """Reads a PDB file and returns a pandas data frame.

    Arguments
    ---------
    pdb_file (string): Path to the PDB file.

    Returns
    -------
    arr_coors (Numpy array): Coordinates of each atom of the PDB file. 
    """
    # List of list containing information about atoms from the pdb file
    rows = []
    with open(pdb_file, "r") as f_in:
        # Parsing through the file 
        for line in f_in:
            #Take the first NMR structure
            if line.startswith("TER"):
                break
            # Sorting by ATOM
            if line.startswith("ATOM"):
                # Extracting informations from the pdb
                atom_num = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22]
                res_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                # Appending those informations into rows list
                rows.append([atom_num, atom_name, res_name, chain_id, res_num, x, y, z])
    # Create a numpy_array containing informations from rows
    arr_coors = pd.DataFrame(rows, columns=["atom_num", "atom_name", "res_name", "chain_id", "res_num", "x", "y", "z"])[["x", "y", "z"]].to_numpy()
    # Return the list and dataframe
    return arr_coors, rows


# Main program
if __name__ == "__main__":
    
    # Taking the pdb_file
    pdb_file = args()

    # Calling parse_pdb
    arr_coors, rows = parse_pdb(pdb_file)

    # Creating a distance_matrix with the numpy_array
    dist_mat = distance_matrix(arr_coors, arr_coors)

    rows_all = []
    for i in range(arr_coors.shape[0]):
        for j in range(i+1, arr_coors.shape[0]):
            if dist_mat[i,j] < 10:
                rows_all.append(rows[i]+rows[j]+[dist_mat[i,j]])

    df_all = pd.DataFrame(rows_all)
    


    #Intraprotein Hydrophobic Interactions
    print("Intraprotein Hydrophobic Interactions\n".center(74))
    print("Hydrophobic Interactions within 5 Angstroms")

    df_hydrophobic = df_all[[4,2,3,12,10,11]][(df_all[2].isin(hdc_aa)) & (df_all[10].isin(hdc_aa)) & 
                                              (df_all[1].str.contains("C[BGDE]")) & (df_all[9].str.contains("C[BGDE]")) &
                                              (df_all[16] <= 5.0) & (df_all[4] != df_all[12])].drop_duplicates()
    if df_hydrophobic.empty:
        print("")
        print("NO INTRAPROTEIN HYDROPHOBIC INTERACTIONS FOUND\n\n".center(74))
    else:
        header_hydrophobic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain"]
        table_hydrophobic = tabulate(df_hydrophobic, headers = header_hydrophobic, showindex=False, numalign="left", tablefmt="rst")
        print(table_hydrophobic, "\n\n\n")



    #Intraprotein Disulphide Bridges
    print("Intraprotein Disulphide Bridges\n".center(74))
    print("Disulphide bridges: Between sulphur atoms of cysteines within 2.2 Angstroms")

    df_disulphide = df_all[[4,2,3,12,10,11,16]][(df_all[2] == "CYS") & (df_all[10] == "CYS") &
                                                (df_all[1] == "SG") & (df_all[9] == "SG") &
                                                (df_all[16] <= 2.2)]
    if df_disulphide.empty:
        print("")
        print("NO INTRAPROTEIN DISULPHIDE BRIDGES FOUND\n\n\n".center(74))
    else:
        header_disulphide = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "Distance"]
        table_disulphide = tabulate(df_disulphide, headers = header_disulphide, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
        print(table_disulphide, "\n\n\n")



    #Intraprotein Ionic Interactions
    print("Intraprotein Ionic Interactions\n".center(74))
    print("Ionic Interactions within 6 Angstroms")


    df_ionic = df_all[[4,2,3,12,10,11]][(df_all[2].isin(ion_aa)) & (df_all[10].isin(ion_aa)) & 
                                        (df_all[1].str.match("[NO][HEZ][^D]*")) & (df_all[9].str.match("[NO][HEZ][D^]*")) &
                                        (df_all[16] < 6) & (df_all[4] != df_all[12])].drop_duplicates()
    if df_disulphide.empty:
        print("")
        print("NO INTRAPROTEIN IONIC INTERACTIONS FOUND\n\n\n".center(74))
    else:
        header_ionic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain"]
        table_ionic = tabulate(df_ionic, headers = header_ionic, showindex=False, numalign="left", tablefmt="rst")
        print(table_ionic, "\n\n\n")





    #Intraprotein Aromatic-Aromatic Interactions
    print("Intraprotein Aromatic-Aromatic Interactions\n".center(74))
    print("Aromatic-Aromatic Interactions within 4.5 and 7 Angstroms")






    centroids_PHE = df_all[[4,5,6,7]][(df_all[2] == "PHE") &
                             (df_all[1].str.contains("C[GDEZ]"))  &
                             (df_all[4] != df_all[12])].drop_duplicates().groupby([4]).mean()
    centroids_TYR = df_all[[1,2,4,5,6,7]][(df_all[2] == "TYR") &
                             (df_all[1].str.contains("C[GDEZ]"))  &
                             (df_all[4] != df_all[12])].drop_duplicates().groupby([4]).mean()
    centroids_TRP = df_all[[1,2,4,5,6,7]][(df_all[2] == "TRP") &
                             (df_all[1].str.contains("C[DEZH][23]"))  &
                             (df_all[4] != df_all[12])].drop_duplicates().groupby([4]).mean()

    centroids_arro = pd.concat([centroids_PHE, centroids_TYR, centroids_TRP])
    print(centroids_arro)


    liste = list(centroids_arro.index)
    print(liste)


    arr_centro = centroids_arro.to_numpy()
    dist_mat_centro = distance_matrix(arr_centro, arr_centro)
    print(dist_mat)

    for i in  range(len(liste)):
        for j in range(i+1, len(liste)):
            if (dist_mat_centro[i,j] < 7) & (dist_mat_centro[i,j] > 4.5):
                print(liste[i], liste[j], dist_mat_centro[i,j])


    """
    centroids_coors = df_all[[1,2,4,5,6,7,9,10,12,13,14,15]][(df_all[2].isin(arom_aa)) & (df_all[10].isin(arom_aa)) &
                             (df_all[1].str.contains("C[GDEZ]")) & (df_all[9].str.contains("C[GDEZ]")) &
                             (df_all[4] != df_all[12])]
    print(centroids_coors)
    df_aromatic = df_all[[4,2,3,12,10,11]][(df_all[2].isin(arom_aa)) & (df_all[10].isin(arom_aa)) & 
                                           (df_all[1].str.contains("C[GDE][12]")) & (df_all[9].str.contains("C[BGDE]")) &
                                           (df_all[16] < 6) & (df_all[4] != df_all[12])].drop_duplicates()
    if df_aromatic.empty:
        print("")
        print("NO INTRAPROTEIN AROMATIC-AROMATIC INTERACTIONS FOUND\n\n\n".center(74))
    else:
        header_ionic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "D(centroid-centroid)", "Dihedral Angle"]
        print(df_aromatic)
    """