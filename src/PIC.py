"""
This script ...

Inputs:
=======


Output:
========


Usage:
======
    $ python PIC.py ...
"""

################################################################################
import os
import sys
import argparse

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from tabulate import tabulate
################################################################################



def args():
    """Parse the command-line arguments.

    Arguments
    ---------
    Python command-line

    Returns
    -------
    pdb_file: string
    """

    #Declaration of expexted arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", help="Path to the PDB file.", type=str, required=True)
    args = parser.parse_args()

    #Check if the PDB files directory is valid
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
    rows = []
    with open(pdb_file, "r") as f_in:
        for line in f_in:
            if line.startswith("ATOM"):
                atom_num = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22]
                res_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                rows.append([atom_num, atom_name, res_name, chain_id, res_num, x, y, z])
    arr_coors = pd.DataFrame(rows, columns=["atom_num", "atom_name", "res_name", "chain_id", "res_num", "x", "y", "z"])[["x", "y", "z"]].to_numpy()
    return arr_coors, rows




if __name__ == "__main__":
    pdb_file = args()


    arr_coors, rows = parse_pdb(pdb_file)

    dist_mat = distance_matrix(arr_coors, arr_coors)

    print("Matrix shape", dist_mat.shape)

    rows_all = []

    for i in range(arr_coors.shape[0]):
        for j in range(i+1, arr_coors.shape[0]):
            if dist_mat[i,j] < 10 and (rows[i][4] != rows[j][4]):
                rows_all.append(rows[i]+rows[j]+[dist_mat[i,j]])

    df_all = pd.DataFrame(rows_all)

    print(df_all)


    df_disulphide = df_all[[4,2,3,12,10,11,16]][(df_all[2] == "CYS") & (df_all[10] == "CYS") & (df_all[1] == "SG") & (df_all[9] == "SG") & (df_all[16] <= 2.2)]
    header_disulphide = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "Distance"]
    table_disulphide = tabulate(df_disulphide, headers = header_disulphide, showindex=False, floatfmt=".2f", tablefmt="rst")
    print(table_disulphide)