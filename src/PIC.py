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
import math

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
            #if line.startswith("TER"):
            #   break
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



def get_dihedral(coor1, coor2, coor3, coor4):
    """
    Calcul d'un angle dihèdre
    """
    vec12 = coor1[0] - coor2[0], coor1[1] - coor2[1], coor1[2] - coor2[2]
    vec23 = coor2[0] - coor3[0], coor2[1] - coor3[1], coor2[2] - coor3[2]
    vec34 = coor3[0] - coor4[0], coor3[1] - coor4[1], coor3[2] - coor4[2]
    
    vec_norm1 = np.cross(vec12, vec23)
    vec_norm2 = np.cross(vec23, vec34)
     
    angle_rad = math.acos(np.dot(vec_norm1, vec_norm2) / (np.linalg.norm(vec_norm1) * np.linalg.norm(vec_norm2)))
    angle_deg = angle_rad * 180 / math.pi

    if np.dot(vec12, vec_norm2) < 0:    #Calcul du produit mixte, pour assigner le bon signe à l'angle.
        return angle_deg
    else:
        return -angle_deg



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


    centroids_PHE_TYR = df_all[[4,5,6,7]][(df_all[2].isin(["PHE", "TYR"])) &
                             (df_all[1].str.contains("C[GDEZ]"))  &
                             (df_all[4] != df_all[12])].drop_duplicates().groupby([4]).mean()
    centroids_TRP = df_all[[4,5,6,7]][(df_all[2] == "TRP") &
                             (df_all[1].str.contains("C[DEZH][23]"))  &
                             (df_all[4] != df_all[12])].drop_duplicates().groupby([4]).mean()

    centroids_arro = pd.concat([centroids_PHE_TYR, centroids_TRP])


    liste = list(centroids_arro.index)
    arr_centro = centroids_arro.to_numpy()
    dist_mat_centro = distance_matrix(arr_centro, arr_centro)


    res_1 = []
    res_2 = []
    list_dist = []

    for i in range(len(liste)):
        for j in range(i+1, len(liste)):
            if (dist_mat_centro[i,j] > 4.5) & (dist_mat_centro[i,j] < 7):
                res_1.append(liste[i])
                res_2.append(liste[j])
                list_dist.append(dist_mat_centro[i,j])

    df_aromatic = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[[4,2,3,12,10,11]][((df_all[4] == res_1[i]) & (df_all[12] == res_2[i])) | ((df_all[4] == res_2[i]) & (df_all[12] == res_1[i]))].drop_duplicates()
        df_aromatic = df_aromatic.append(row)


    #Insert the dist column
    df_aromatic["D(centroid-centroid)"] = list_dist



    #Calculate the dihedral angles NOT WORKING
    res_1 = df_aromatic.iloc[:,0].to_numpy()
    res_2 = df_aromatic.iloc[:,3].to_numpy()

    dihedral_angle = []

    for i in range(len(res_1)):
        if df_aromatic.iloc[i,1] == "PHE" or df_aromatic.iloc[i,1] == "TYR":
            C1 = df_all[[5,6,7]][(df_all[4] == res_1[i]) & (df_all[1].str.match("CG"))].drop_duplicates().to_numpy().flatten()
            C2 = df_all[[5,6,7]][(df_all[4] == res_1[i]) & (df_all[1].str.match("CZ"))].drop_duplicates().to_numpy().flatten()
        if df_aromatic.iloc[i,4] == "PHE" or df_aromatic.iloc[i,4] == "TYR":
            C3 = df_all[[5,6,7]][(df_all[4] == res_2[i]) & (df_all[1].str.match("CZ"))].drop_duplicates().to_numpy().flatten()
            C4 = df_all[[5,6,7]][(df_all[4] == res_2[i]) & (df_all[1].str.match("CG"))].drop_duplicates().to_numpy().flatten()
        if df_aromatic.iloc[i,1] == "TRP":
            C1 = df_all[[5,6,7]][(df_all[4] == res_1[i]) & (df_all[1].str.match("CE3"))].drop_duplicates().to_numpy().flatten()
            C2 = df_all[[5,6,7]][(df_all[4] == res_1[i]) & (df_all[1].str.match("CZ2"))].drop_duplicates().to_numpy().flatten()
        if df_aromatic.iloc[i,4] == "TRP":
            C3 = df_all[[5,6,7]][(df_all[4] == res_2[i]) & (df_all[1].str.match("CZ2"))].drop_duplicates().to_numpy().flatten()
            C4 = df_all[[5,6,7]][(df_all[4] == res_2[i]) & (df_all[1].str.match("CE3"))].drop_duplicates().to_numpy().flatten()

        angle = get_dihedral(C1, C2, C3, C4)
        dihedral_angle.append(angle)

    #Insert the dihedral angle column
    df_aromatic["Dihedral Angle"] = dihedral_angle


    print(df_aromatic)





    #Intraprotein Aromatic-Sulphur Interactions

    centroids_PHE_TYR = df_all[[4,5,6,7]][(df_all[2].isin(["PHE", "TYR"])) &
                             (df_all[1].str.contains("C[GDEZ]"))  &
                             (df_all[4] != df_all[12])].drop_duplicates().groupby([4]).mean()
    centroids_TRP = df_all[[4,5,6,7]][(df_all[2] == "TRP") &
                             (df_all[1].str.contains("C[DEZH][23]"))  &
                             (df_all[4] != df_all[12])].drop_duplicates().groupby([4]).mean()

    centroids_arro = pd.concat([centroids_PHE_TYR, centroids_TRP])


    liste = list(centroids_arro.index)
    arr_centro = centroids_arro.to_numpy()


    df_coors_s = df_all[[4,5,6,7]][(df_all[2] == "CYS") & (df_all[1] == "SG")].drop_duplicates().groupby([4]).mean()

    index_s = list(df_coors_s.index)
    arr_coors_s = df_coors_s.to_numpy()


    dist_mat_centro_s = distance_matrix(arr_centro, arr_coors_s)


    res_1 = []
    res_2 = []
    list_dist = []

    for i in range(len(liste)):
        for j in range(len(index_s)):
            if dist_mat_centro_s[i,j] < 5.3:
                res_1.append(liste[i])
                res_2.append(index_s[j])
                list_dist.append(dist_mat_centro_s[i,j])
                print(liste[i], index_s[j], dist_mat_centro_s[i,j])

    df_aromatic = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[[4,2,3,12,10,11]][((df_all[4] == res_1[i]) & (df_all[12] == res_2[i])) | ((df_all[4] == res_2[i]) & (df_all[12] == res_1[i]))].drop_duplicates()
        df_aromatic = df_aromatic.append(row)