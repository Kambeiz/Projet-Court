"""
This script calculates intra-proteins interactions in a given PDB file.

This script, given the coordinates of protein 3D structure, computes  various 
interactions such as disulphide bonds, interactions between hydrophobic 
residues, ionic interactions, hydrogen bonds, aromatic-aromatic interactions, 
aromatic-sulphur interactions and cation-pi interactions within a protein or 
between proteins in a complex.


Input:
=======
PDB file

Output:
========
Protein interaction tables

Usage:
======
    $ python PIC.py -p pdb_file
"""

__author__ = "Debbah Nagi, Vander Meersche Yann"

__license__ = "M2BI"
__version__ = "1.0.0"
__date__ = "2020-09-09"
__email__ = "debbah.nag@gmail.com, yann-vm@hotmail.fr"
__copyright__ = "Copyright 2020, The Short Project Inc"


# MODULES ######################################################################
import os
import sys
import argparse
import re
import math

import numpy as np
import pandas as pd
from tabulate import tabulate
from scipy.spatial import distance_matrix
from selenium import webdriver
################################################################################


# CONSTANTS ####################################################################
hdc_aa = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"]
ani_aa = ["ASP", "GLU"]
cat_aa = ["ARG", "LYS"]
ion_aa = ["ASP", "GLU", "ARG", "HIS", "LYS"]
sulph_aa = ["CYS", "MET"]
################################################################################



def args():
    """
    Parse the command-line arguments.

    Arguments
    ---------
    Python command-line

    Returns
    -------
    pdb_file (string): Path to the PDB file
    """
    # Declaration of expexted arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", help="Path to the PDB file.", type=str, required=True)
    args = parser.parse_args()

    # Check if the PDB file directory is valid
    pdb_file = args.pdb
    if not os.path.isfile(pdb_file):
        sys.exit(f"{pdb_file} does not exist.\n"
                "Please enter a valid PDB file.")

    return pdb_file


def parse_pdb(pdb_file):
    """
    Reads a PDB file and returns a Pandas dataframe.

    Arguments
    ---------
    pdb_file (string): Path to the PDB file

    Returns
    -------
    arr_coors (NumPy array): Coordinates of each atom of the PDB file
    """
    # List of list containing information about atoms from the PDB file
    rows = []

    with open(pdb_file, "r") as f_in:
        # Parse through the file 
        for line in f_in:
            # If requiered take the first NMR structure
            if line.startswith("ENDMDL"):
               break
            # Extractes informations from the PDB
            if line.startswith("ATOM"):
                atom_num = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22]
                res_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                # Appends these informations into a list
                rows.append([atom_num, atom_name, res_name, chain_id, res_num, x, y, z])

    # Creates a NumPy array containing atoms coordinates
    arr_coors = pd.DataFrame(rows, columns=["atom_num", "atom_name", "res_name", "chain_id", "res_num", "x", "y", "z"])[["x", "y", "z"]].to_numpy()

    return arr_coors, rows


def get_dihedral(coor1, coor2, coor3, coor4):
    """
    Dihedral angle calculation.

    Arguments
    ---------
    coors* (NumPy array): 3D Coordinate of an atom

    Returns
    -------
    angle_deg (float): Dihedral angle of the plans
    """
    vec12 = coor1[0] - coor2[0], coor1[1] - coor2[1], coor1[2] - coor2[2]
    vec23 = coor2[0] - coor3[0], coor2[1] - coor3[1], coor2[2] - coor3[2]
    vec34 = coor3[0] - coor4[0], coor3[1] - coor4[1], coor3[2] - coor4[2]
    
    vec_norm1 = np.cross(vec12, vec23)
    vec_norm2 = np.cross(vec23, vec34)
     
    angle_rad = math.acos(np.dot(vec_norm1, vec_norm2) / (np.linalg.norm(vec_norm1) * np.linalg.norm(vec_norm2)))
    angle_deg = angle_rad * 180 / math.pi

    #Calculate the dot product to determine the right sign
    if np.dot(vec12, vec_norm2) < 0:
        return angle_deg
    else:
        return -angle_deg


def launching_HBONDS(pdbfile):
    """
    Open PIC in a web-browser (Firefox) and get the Hbond output text.

    Arguments
    ---------
    pdbfile (string): Path to the PDB file

    Returns
    -------
    body (string): Body of the HTML page
    """
    url = "http://pic.mbu.iisc.ernet.in/job.html"
    list_elements = ["hbond3", "hbond4", "hbond5", "Submit"]

    firefox_options = webdriver.FirefoxOptions()  
    firefox_options.add_argument("--headless")
    # Chose between openind Firefox and not opening it
    driver = webdriver.Firefox(options=firefox_options)    #Graphical
    #driver = webdriver.Firefox()                          #Hidden
    
    driver.get(url)
    upload = driver.find_element_by_name("pdbname")

    absolute_path = os.path.abspath(pdbfile)
    upload.send_keys(absolute_path)

    for elem in list_elements:
        click_elem = driver.find_element_by_name(elem)
        click_elem.click()

    body = (driver.find_element_by_xpath("//body").text).split("\n")
    driver.close()

    return body

'''
    def launching_HBONDS(pdbfile):
    """
    Open PIC in a browser and get the Hbond output.
    Arguments
    ---------
    pdbfile (string): Path to the PDB file
    Returns
    -------
    body (string): Body of the HTML page
    """
    url = "http://pic.mbu.iisc.ernet.in/job.html"
    url_hbd = "http://pic.mbu.iisc.ernet.in/TEMP/" + os.path.basename(pdbfile).split(".")[0] + ".hbd" 
    list_elements = ["hbond3", "hbond4", "hbond5", "Submit"]
    firefox_options = webdriver.FirefoxOptions()  
    firefox_options.add_argument("--headless")
    # Chose between openind Firefox and not opening it
    #driver = webdriver.Firefox(options=firefox_options)    #Hidden 
    driver = webdriver.Firefox()                          #Graphical
    
    driver.get(url)
    upload = driver.find_element_by_name("pdbname")
    absolute_path = os.path.abspath(pdbfile)
    upload.send_keys(absolute_path)
    for elem in list_elements:
        click_elem = driver.find_element_by_name(elem)
        click_elem.click()
    driver.get(url_hbd)
    body = (driver.find_element_by_xpath("//body").text).split("\n")
    driver.close()
    return body
'''

def body_to_list(body):
    """
    Parse the Hbond output into 3 Pandas dataframes.

    Arguments
    ---------
    body (string): Body of the HTML page

    Returns
    -------
    df_ch_main, df_ch_side, df_ch_side (dataframes): PIC Hydrogen bond interactions
    """
    count = 0
    main_ch_main = []
    main_ch_side = []
    side_ch_side = []
    

    ####################
    for line in body:
        if "Main Chain-Main" in line:
            main_main = True
        elif "Main Chain-Side" in line:
            main_side = True
        elif "Side Chain-Side" in line:
            side_side = True
        if len(line) > 9:
            if line[:1].isdigit() and main_main == True:
                main_ch_main.append(line)
            elif line.startswith("Dd-a") and count == 0:
                main_main = False
                count += 1
            elif line[:1].isdigit() and main_side == True:
                main_ch_side.append(line)
            elif line.startswith("Dd-a") and count == 1:
                main_side = False
                count += 1
            elif line[:1].isdigit() and side_side == True:
                side_ch_side.append(line)
            elif line.startswith("Dd-a") and count == 2:
                side_side=False

    return list_to_df(main_ch_main), list_to_df(main_ch_side), list_to_df(side_ch_side)    


def list_to_df(list_interactions):
    """
    Convert a list of list
    """
    df = pd.DataFrame(list_interactions)
    df = df[0].str.split(pat=" ", expand=True)

    return df



# Main program
if __name__ == "__main__":
    
    # Taking the pdb_file and commandline options / parameters
    pdb_file = args()

    # Calling parse_pdb
    arr_coors, rows_list = parse_pdb(pdb_file)

    # Creating a distance_matrix with the NumPy array
    dist_mat = distance_matrix(arr_coors, arr_coors)

    # Create list with the possible pairwise atom interactions in a given threshold
    rows_all = []
    for i in range(arr_coors.shape[0]):
        for j in range(i+1, arr_coors.shape[0]):
            # Add a row to the list only if the distance between them is below 
            # the threshold, and if atoms are form a different residue
            if (dist_mat[i,j] < 10):
                # [Atom 1] [Atom 2] [Distance between Atom 1 and Atom 2]
                rows_all.append(rows_list[i]+rows_list[j]+[dist_mat[i,j]])
    # Convert the list of list into a dataframe
    df_all = pd.DataFrame(rows_all, columns = ["atom_num1", "atom_name1", "res_name1", "chain_id1", "res_num1", "x1", "y1", "z1",
                                               "atom_num2", "atom_name2", "res_name2", "chain_id2", "res_num2", "x2", "y2", "z2",
                                               "atom_dist"])
    


    # Intraprotein Hydrophobic Interactions ####################################
    print("")
    print("Intraprotein Hydrophobic Interactions\n".center(106))
    print("Hydrophobic Interactions within 5 Angstroms")


    df_hydrophobic = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(df_all["res_name1"].isin(hdc_aa)) & (df_all["res_name2"].isin(hdc_aa)) & 
                                                                                                     (df_all["atom_name1"].str.contains("C[BGDE]")) & (df_all["atom_name2"].str.contains("C[BGDE]")) &
                                                                                                     (df_all["atom_dist"] <= 5.0) & (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates()
    if df_hydrophobic.empty:
        print("")
        print("NO INTRAPROTEIN HYDROPHOBIC INTERACTIONS FOUND\n\n".center(106))
    else:
        header_hydrophobic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain"]
        table_hydrophobic = tabulate(df_hydrophobic, headers = header_hydrophobic, showindex=False, numalign="left", tablefmt="rst")
        print(table_hydrophobic, "\n\n\n")



    # Intraprotein Disulphide Bridges ##########################################
    print("Intraprotein Disulphide Bridges\n".center(106))
    print("Disulphide bridges: Between sulphur atoms of cysteines within 2.2 Angstroms")

    df_disulphide = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2","atom_dist"]][(df_all["res_name1"] == "CYS") & (df_all["res_name2"] == "CYS") &
                                                                                                                (df_all["atom_name1"] == "SG") & (df_all["atom_name2"] == "SG") &
                                                                                                                (df_all["atom_dist"] <= 2.2)]
    if df_disulphide.empty:
        print("")
        print("NO INTRAPROTEIN DISULPHIDE BRIDGES FOUND\n\n\n".center(106))
    else:
        header_disulphide = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "Distance"]
        table_disulphide = tabulate(df_disulphide, headers = header_disulphide, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
        print(table_disulphide, "\n\n\n")



    #============================= HBOND RESULTS ==============================#
    page_results = launching_HBONDS(pdb_file)
    df_main_main, df_main_side, df_side_side = body_to_list(page_results)
    legend_hbond = """
Dd-a     =   Distance Between Donor and Acceptor
Dh-a     =   Distance Between Hydrogen and Acceptor
A(d-H-N) =   Angle Between Donor-H-N
A(a-O=C) =   Angle Between Acceptor-O=C
MO       =   Multiple Occupancy
Note that angles that are undefined are written as 999.99


"""

    #Intraprotein Main Chain-Main Chain Hydrogen Bonds #########################
    print("Intraprotein Main Chain-Main Chain Hydrogen Bonds\n".center(106))

    if df_main_main.empty:
        print("")
        print("NO INTRAPROTEIN MAIN CHAIN-MAIN CHAIN HYDROGEN BONDS FOUND\n\n\n".center(106))
    else:
        header_main_main = ["POS", "CHAIN", "RES", "ATOM", "POS", "CHAIN", "RES", "ATOM", "Dd-a", "Dh-h", "A(d-H-N)", "A(a-O=C)"]
        table_main_main = tabulate(df_main_main, headers = header_main_main, showindex=False, numalign="left", tablefmt="rst")
        print("            DONOR                         ACCEPTOR                          PARAMETERS              ")
        print(table_main_main)
        print(legend_hbond)


    # Intraprotein Main Chain-Side Chain Hydrogen Bonds ########################
    print("Intraprotein Main Chain-Side Chain Hydrogen Bonds\n".center(106))

    if df_main_side.empty:
        print("")
        print("NO INTRAPROTEIN MAIN CHAIN-SIDE CHAIN HYDROGEN BONDS FOUND\n\n\n".center(106))
    else:
        header_main_side = ["POS", "CHAIN", "RES", "ATOM", "POS", "CHAIN", "RES", "ATOM", "MO", "Dd-a", "Dh-h", "A(d-H-N)", "A(a-O=C)"]
        table_main_side = tabulate(df_main_side, headers = header_main_side, showindex=False, numalign="left", tablefmt="rst")
        print("             DONOR                         ACCEPTOR                             PARAMETERS                ")
        print(table_main_side)
        print(legend_hbond)


    # Intraprotein Side Chain-Side Chain Hydrogen Bonds ########################
    print("Intraprotein Side Chain-Side Chain Hydrogen Bonds\n".center(106))

    if df_side_side.empty:
        print("")
        print("NO INTRAPROTEIN SIDE CHAIN-SIDE CHAIN HYDROGEN BONDS FOUND\n\n\n".center(106))
    else:
        header_side_side = ["POS", "CHAIN", "RES", "ATOM", "POS", "CHAIN", "RES", "ATOM", "MO", "Dd-a", "Dh-h", "A(d-H-N)", "A(a-O=C)"]
        table_side_side = tabulate(df_side_side, headers = header_side_side, showindex=False, numalign="left", tablefmt="rst")
        print("             DONOR                         ACCEPTOR                             PARAMETERS                ")
        print(table_side_side)
        print(legend_hbond)



    # Intraprotein Ionic Interactions ##########################################
    print("Intraprotein Ionic Interactions\n".center(106))
    print("Ionic Interactions within 6 Angstroms")

    df_ionic = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(df_all["res_name1"].isin(ion_aa)) & (df_all["res_name2"].isin(ion_aa)) & 
                                                                                               (df_all["atom_name1"].str.match("[NO][HEZ][^D]*")) & (df_all["atom_name2"].str.match("[NO][HEZ][D^]*")) &
                                                                                               (df_all["atom_dist"] < 6) & (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates()
    if df_ionic.empty:
        print("")
        print("NO IONIC INTERACTIONS FOUND\n\n\n".center(106))
    else:
        header_ionic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain"]
        table_ionic = tabulate(df_ionic, headers = header_ionic, showindex=False, numalign="left", tablefmt="rst")
        print(table_ionic, "\n\n\n")



    # Intraprotein Aromatic-Aromatic Interactions ##############################
    print("Intraprotein Aromatic-Aromatic Interactions\n".center(106))
    print("Aromatic-Aromatic Interactions within 4.5 and 7 Angstroms")


    # Calculate the aromatics centroids coordinates
    centroids_PHE_TYR = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"].isin(["PHE", "TYR"])) &
                                                                        (df_all["atom_name1"].str.contains("C[GDEZ]")) & (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()
    centroids_TRP = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"] == "TRP") &
                                                                    (df_all["atom_name1"].str.contains("C[DEZH][23]")) & (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()


    # Creates a dataframe with every aromatics centroids
    centroids_arom = pd.concat([centroids_PHE_TYR, centroids_TRP])
    # Grabs the index (chain_id, res_num)
    index_a = list(centroids_arom.index)
    # Transforms the dataframe into a coordinate array
    arr_centro = centroids_arom.to_numpy()


    # Calculates the distance matrix between centroids
    distmat_centro = distance_matrix(arr_centro, arr_centro)

    res_1, res_2 = [], []
    chain_1, chain_2 = [], []
    list_dist = []
    #Selects the residues names and chain ID of every pairs of centroids in the correct threshold
    for i in range(len(index_a)):
        for j in range(i+1, len(index_a)):
            if (distmat_centro[i,j] > 4.5) & (distmat_centro[i,j] < 7):
                res_1.append(index_a[i][1])
                res_2.append(index_a[j][1])
                chain_1.append(index_a[i][0])
                chain_2.append(index_a[j][0])
                list_dist.append(distmat_centro[i,j])

    #Selects the lines of the dataframe corresponding to these pairs
    df_aromatic = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & 
                                                                                               ((df_all["res_num2"] == res_2[i]) & (df_all["chain_id2"] == chain_2[i]))) | 
                                                                                              (((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) &
                                                                                               ((df_all["res_num2"] == res_1[i]) & (df_all["chain_id2"] == chain_1[i])))].drop_duplicates()
        df_aromatic = df_aromatic.append(row)


    if df_aromatic.empty:
        print("")
        print("NO INTRAPROTEIN AROMATIC-AROMATIC INTERACTIONS FOUND\n\n\n".center(106))
    else:
        #Insert the dist column
        df_aromatic["D(centroid-centroid)"] = list_dist

        """
        #Calculate the dihedral angles NOT WORKING
        res_1 = df_aromatic.iloc[:,0].to_numpy()
        res_2 = df_aromatic.iloc[:,3].to_numpy()

        # Calculate the dihedral angles between 2 points of each aromatic cycle
        dihedral_angle = []
        for i in range(len(res_1)):
            if df_aromatic.iloc[i,1] == "PHE" or df_aromatic.iloc[i,1] == "TYR":
                C1 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & (df_all["atom_name1"].str.match("CG"))].drop_duplicates().to_numpy().flatten()
                C2 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & (df_all["atom_name1"].str.match("CZ"))].drop_duplicates().to_numpy().flatten()
            if df_aromatic.iloc[i,1] == "PHE" or df_aromatic.iloc[i,4] == "TYR":
                C3 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & (df_all["atom_name1"].str.match("CZ"))].drop_duplicates().to_numpy().flatten()
                C4 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & (df_all["atom_name1"].str.match("CG"))].drop_duplicates().to_numpy().flatten()
            if df_aromatic.iloc[i,1] == "TRP":
                C1 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & (df_all["atom_name1"].str.match("CE3"))].drop_duplicates().to_numpy().flatten()
                C2 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & (df_all["atom_name1"].str.match("CZ2"))].drop_duplicates().to_numpy().flatten()
            if df_aromatic.iloc[i,4] == "TRP":
                C3 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & (df_all["atom_name1"].str.match("CZ2"))].drop_duplicates().to_numpy().flatten()
                C4 = df_all[["x1","y1","z1"]][((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & (df_all["atom_name1"].str.match("CE3"))].drop_duplicates().to_numpy().flatten()

            angle = get_dihedral(C1, C2, C3, C4)
            dihedral_angle.append(angle)

        #Insert the dihedral angle column
        df_aromatic["Dihedral Angle"] = dihedral_angle
        """
  
        #Creates a beautiful table with tabulate :)
        header_aromatic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "D(centroid-centroid)", "Dihedral Angle"]
        table_aromatic = tabulate(df_aromatic, headers = header_aromatic, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
        print(table_aromatic, "\n\n\n")



    # Intraprotein Aromatic-Sulphur Interactions ###############################
    print("Intraprotein Aromatic-Sulphur Interactions\n".center(106))
    print("Aromatic-Sulphur Interactions within 5.3 Angstroms")

    # Calculates the distance matrix between aromatics centroids and slphur atoms
    df_coors_s = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"].isin(sulph_aa)) & (df_all["atom_name1"].str.contains("S"))].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()
    # Grabs the index (chain_id, res_num)
    index_s = list(df_coors_s.index)
    # Transforms the dataframe into a coordinate array
    arr_coors_s = df_coors_s.to_numpy()

    # Calculates the distance matrix between aromatic centroids and sulphur atoms
    distmat_centro_sulphur = distance_matrix(arr_centro, arr_coors_s)

    res_1, res_2 = [], []
    chain_1, chain_2 = [], []
    list_dist = []
    #Selects the residues names and chain ID of every pairs of centroids/sulphurs in the correct threshold
    for i in range(len(index_a)):
        for j in range(len(index_s)):
            if distmat_centro_sulphur[i,j] < 5.3:
                res_1.append(index_a[i][1])
                res_2.append(index_s[j][1])
                chain_1.append(index_a[i][0])
                chain_2.append(index_s[j][0])
                list_dist.append(distmat_centro_sulphur[i,j])

    #Selects the lines of the dataframe corresponding to these pairs
    df_aromatic_sulphur = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & 
                                                                                               ((df_all["res_num2"] == res_2[i]) & (df_all["chain_id2"] == chain_2[i]))) | 
                                                                                              (((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & 
                                                                                               ((df_all["res_num2"] == res_1[i]) & (df_all["chain_id2"] == chain_1[i])))].drop_duplicates()
        df_aromatic_sulphur = df_aromatic_sulphur.append(row)



    if df_aromatic_sulphur.empty:
        print("")
        print("NO INTRAPROTEIN AROMATIC-SULPHUR INTERACTIONS FOUND\n\n\n".center(106))
    else:
        #Insert the dist column
        df_aromatic_sulphur["D(centroid-centroid)"] = list_dist

        header_aromatic_s = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "D(Centroid-Sulphur)", "Angle"]
        table_aromatic_s = tabulate(df_aromatic_sulphur, headers = header_aromatic_s, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
        print(table_aromatic_s, "\n\n\n")



    #Intraprotein Cation-Pi Interactions #######################################
    print("Intraprotein Cation-Pi Interactions\n".center(106))
    print("Cation-Pi Interactions within 6 Angstroms")


    df_coors_i = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"].isin(cat_aa)) & (df_all["atom_name1"].isin(["NH1", "NH2", "NZ"]))].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()


    index_i = list(df_coors_i.index)
    arr_coors_i = df_coors_i.to_numpy()

    dist_mat_centro_i = distance_matrix(arr_centro, arr_coors_i)

    res_1, res_2 = [], []
    chain_1, chain_2 = [], []
    list_dist = []
    #Selects the residues names and chain ID of every pairs of centroids/cations in the correct threshold
    for i in range(len(index_a)):
        for j in range(len(index_i)):
            if dist_mat_centro_i[i,j] < 6:
                res_1.append(index_a[i][1])
                res_2.append(index_i[j][1])
                chain_1.append(index_a[i][0])
                chain_2.append(index_i[j][0])
                list_dist.append(dist_mat_centro_i[i,j])

    df_arom_i = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & 
                                                                                               ((df_all["res_num2"] == res_2[i]) & (df_all["chain_id2"] == chain_2[i]))) | 
                                                                                              (((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & 
                                                                                               ((df_all["res_num2"] == res_1[i]) & (df_all["chain_id2"] == chain_1[i])))].drop_duplicates()
        df_arom_i = df_arom_i.append(row)



    if df_arom_i.empty:
        print("")
        print("NO INTRAPROTEIN CATION-PI INTERACTIONS FOUND\n\n\n".center(106))
    else:
        #Insert the dist column
        df_arom_i["D(cation-Pi)"] = list_dist

        header_arom_i = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "D(cation-Pi)", "Angle"]
        table_arom_i = tabulate(df_arom_i, headers = header_arom_i, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
        print(table_arom_i, "\n")
