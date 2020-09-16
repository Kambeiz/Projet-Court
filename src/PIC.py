"""
This script calculates intra- and inter-proteins interactions in a given PDB file.

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
    $ python PIC.py [-h] -pdb PDBFILE [-intra | -inter] [-hydro HYDROPHOBIC] [-ion IONIC] [-AA A A] [-AS AROMSULPH] [-AC AROMCATION]

"""

__author__ = "Debbah Nagi, Vander Meersche Yann"

__license__ = "M2-BI"
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

import selenium
import numpy as np
import pandas as pd
from tabulate import tabulate
from scipy.spatial import distance_matrix
from selenium import webdriver
################################################################################


# CONSTANTS ####################################################################
hdc_aa = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"]
ani_aa = ["ASP", "GLU"]
cat_aa = ["ARG", "LYS", "HIS"]
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
    parser.add_argument("-pdb", "--pdbfile", help="Path to the PDB file.", type=str, required=True)
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-intra", "--intrachain", help="Use this option to only print the intra-chain interactions.", action="store_true")
    group.add_argument("-inter", "--interchain", help="Use this option to only print the inter-chain interactions.", action="store_true")
    parser.add_argument("-hydro", "--hydrophobic", help="Enter the interaction cut-off value (Default 5A).", type=float, required=False, default=5)
    parser.add_argument("-ion", "--ionic", help="Enter the interaction cut-off value (Default 6A).", type=float, required=False, default=6)
    parser.add_argument("-AA", "--aromarom", metavar="A", type=float, nargs=2, help="Enter the interaction cut-off value (Default 4.5A to 7A).", required=False, default=(4.5, 7))
    parser.add_argument("-AS", "--aromsulph", help="Enter the interaction cut-off value (Default 5.3A).", type=float, required=False, default=5.3)
    parser.add_argument("-AC", "--aromcation", help="Enter the interaction cut-off value (Default 6A).", type=float, required=False, default=6)
    args = parser.parse_args()

    # Check if the PDB file directory is valid
    pdb_file = args.pdbfile
    if not os.path.isfile(pdb_file):
        sys.exit(f"{pdb_file} does not exist.\n"
                "Please enter a valid PDB file.")


    # Arguments cutoff declaration
    hydro_cutoff = args.hydrophobic
    ionic_cutoff = args.ionic
    aromarom_cutoff_min, aromarom_cutoff_max = args.aromarom
    aromsulph_cutoff = args.aromsulph
    aromcation_cutoff = args.aromcation

    # Interaction type declaration
    intrachain = args.intrachain
    interchain = args.interchain

    return pdb_file, hydro_cutoff, ionic_cutoff, aromarom_cutoff_min, aromarom_cutoff_max, aromsulph_cutoff, aromcation_cutoff, intrachain, interchain


def parse_pdb(pdb_file):
    """
    Reads a PDB file and returns a list of list (PDB file lines) and a NumPy array (PDB file coordinates).

    Arguments
    ---------
    pdb_file (string): Path to the PDB file

    Returns
    -------
    arr_coors (NumPy array): Coordinates of each atom of the PDB file
    rows (list): All sequence information
    """
    # List of list containing information about atoms from the PDB file
    rows = []

    with open(pdb_file, "r") as f_in:
        # Go through the file 
        for line in f_in:
            # If requiered take the first NMR structure
            if line.startswith("ENDMDL"):
               break
            # Extracts informations from the PDB
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

    # Create a NumPy array containing atoms coordinates
    arr_coors = pd.DataFrame(rows, columns=["atom_num", "atom_name", "res_name", "chain_id", "res_num", "x", "y", "z"])[["x", "y", "z"]].to_numpy()

    return arr_coors, rows


def launching_HBONDS(pdbfile):
    """
    Open PIC in a browser (headless) and get the HBOND file (.hbd).

    Arguments
    ---------
    pdbfile (string): Path to the PDB file
    
    Returns
    -------
    body (string): HBOND output file (.hbd)
    """
    
    # URL from where we are launching HBOND
    url = "http://pic.mbu.iisc.ernet.in/job.html"
    # URL where we will get our hydrogen interaction results
    url_hbd = "http://pic.mbu.iisc.ernet.in/TEMP/" + os.path.basename(pdbfile).split(".")[0] + ".hbd" 
    
    # Elements name that we will use in the page form for choosing our interactions
    list_elements = ["hbond3", "hbond4", "hbond5", "Submit"]
    
    # Options for the browser 
    firefox_options = webdriver.FirefoxOptions()  
    firefox_options.add_argument("--headless")
    # Chose between opening a window of Firefox and not opening it 
    driver = webdriver.Firefox(options=firefox_options)    #Hidden 
    #driver = webdriver.Firefox()                          #Graphical
    
    
    # Going to the chosen URL
    # Handling exception if the page is not loaded correctly
    try:
        driver.get(url)
    # Return a empty file, and close the browser
    except:
        print("Unable to extract HBOND output file.\n(PIC might be on heavy load, please try later. Are you connected to internet?)\n\n".center(106))
        driver.close()
        body = ""

        return body

    
    # Search the button for sendig our PDB file
    upload = driver.find_element_by_name("pdbname")
    
    # Get absolute the path of our PDB file
    absolute_path = os.path.abspath(pdbfile)
    
    # Send the path of the file into the right form
    upload.send_keys(absolute_path)
    
    # Choose our interactions
    for elem in list_elements:
        click_elem = driver.find_element_by_name(elem)
        click_elem.click()
    
    # Move to the result page
    driver.get(url_hbd)
    
    # Extract the body text and split it in a list of line
    body = (driver.find_element_by_xpath("//body").text).split("\n")
    
    # Close the browser
    driver.close()

    return body


def body_to_list(body):
    """
    Parse the HBOND output into a Pandas dataframe.

    Arguments
    ---------
    body (string): HBOND output file (.hbd)

    Returns
    -------
    df_hbond (dataframe): PIC Hydrogen bond interactions
    """
    
    if body == "":
        df_hbond = pd.DataFrame()

        return df_hbond

    else:
        # Creating a list that will contain our information sorted
        list_hbond = []
                       
        # Ignore header and loop over the body text
        for line in body:
            if not line.startswith("#"):
                # Extract information per position in line
                res_num_d = int(line[4:9]) 
                chain_id_d = line[11]
                res_name_d = line[13]
                atom_name_d = line[15:18].strip()
                res_num_a = int(line[23:27])
                chain_id_a = line[29]
                res_name_a = line[31]
                atom_name_a = line[33:36].strip()
                typ = line[37:39]
                MO = line[40]
                Dd_a = float(line[46:50])
                Dh_a = float(line[51:55])
                dHN = float(line[56:62])
                aOC = float(line[63:69])
    
                # Store those informations into a list of list
                list_hbond.append([res_num_d, chain_id_d, res_name_d, atom_name_d, res_num_a, chain_id_a, res_name_a, atom_name_a, MO, typ, Dd_a, Dh_a, dHN, aOC])
    
        aa_dict = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
                   "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
        # Cast the list of lists into a Pandas dataframe and replace the one letter amino acids column by it's 3 letter code
        df_hbond = pd.DataFrame(list_hbond).replace({2: aa_dict, 6: aa_dict})

        return df_hbond




################################# MAIN PROGRAM #################################

if __name__ == "__main__":

    # Take the pdb_file and command line options/parameters
    pdb_file, hydro_cutoff, ionic_cutoff, aromarom_cutoff_min, aromarom_cutoff_max, aromsulph_cutoff, aromcation_cutoff, intrachain, interchain = args()

    # Setup variable to print the correct type on interaction calculation:
    if intrachain:
        intrainter = "Intraprotein"
        INTRAINTER = "INTRAPROTEIN"
    elif interchain:
        intrainter = "Interprotein"
        INTRAINTER = "INTERPROTEIN"
    else:
        intrainter = "Intraprotein & Interprotein"
        INTRAINTER = "INTRAPROTEIN & INTERPROTEIN"

    # Parse the PDB file
    arr_coors, rows_list = parse_pdb(pdb_file)

    # Create a distance_matrix with the coordinates NumPy arrays
    dist_mat = distance_matrix(arr_coors, arr_coors)

    # Create list with the possible pairwise atom interactions within a given threshold
    rows_all = []

    for i in range(arr_coors.shape[0]):
        for j in range(i+1, arr_coors.shape[0]):
            # Add a row to the list only if the distance between them is below the threshold (speed-up calculation), and if atoms are form a different residue
            if (dist_mat[i,j] < max(hydro_cutoff, ionic_cutoff, aromarom_cutoff_min, aromarom_cutoff_max, aromsulph_cutoff, aromcation_cutoff) + 3): # +3 to  be sure to select all residues in arom_arom interactions
                # [Atom 1] [Atom 2] [Distance between Atom 1 and Atom 2]
                rows_all.append(rows_list[i]+rows_list[j]+[dist_mat[i,j]])
    # Convert the list of list into a dataframe
    df_all = pd.DataFrame(rows_all, columns = ["atom_num1", "atom_name1", "res_name1", "chain_id1", "res_num1", "x1", "y1", "z1",
                                               "atom_num2", "atom_name2", "res_name2", "chain_id2", "res_num2", "x2", "y2", "z2",
                                               "atom_dist"])
    


    # Hydrophobic Interactions #################################################
    print("")
    print(f"{intrainter} Hydrophobic Interactions\n".center(106))
    print(f"Hydrophobic Interactions within {hydro_cutoff} Angstroms")

    # Create a dataframe with hydrophobic interactions
    df_hydrophobic = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(df_all["res_name1"].isin(hdc_aa)) & (df_all["res_name2"].isin(hdc_aa)) &    # Select the correct amino acids pairs
                                                                                                     (df_all["atom_name1"].str.match("C[BGDEHZ]+")) & (df_all["atom_name2"].str.match("C[BGDEHZ]+")) &    # Select the correct atoms
                                                                                                     (df_all["atom_dist"] <= hydro_cutoff) & (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates()    # Select pairs within the cutoff with are not in the same residue
    
    # Use an option to select the correct subset of data   
    if intrachain:
        df_hydrophobic = df_hydrophobic[(df_hydrophobic["chain_id1"] == df_hydrophobic["chain_id2"])]
    elif interchain:
        df_hydrophobic = df_hydrophobic[(df_hydrophobic["chain_id1"] != df_hydrophobic["chain_id2"])]

    # Check if the dataframe is empty and print accordingly
    if df_hydrophobic.empty:
        print("")
        print(f"NO {INTRAINTER} HYDROPHOBIC INTERACTIONS FOUND\n\n".center(106))
    else:
        header_hydrophobic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain"]
        table_hydrophobic = tabulate(df_hydrophobic, headers = header_hydrophobic, showindex=False, numalign="left", tablefmt="rst")
        print(table_hydrophobic, "\n\n\n")



    # Disulphide Bridges #######################################################
    print(f"{intrainter} Disulphide Bridges\n".center(106))
    print("Disulphide bridges: Between sulphur atoms of cysteines within 2.2 Angstroms")

    # Create a dataframe with disulphide bridges
    df_disulphide = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2","atom_dist"]][(df_all["res_name1"] == "CYS") & (df_all["res_name2"] == "CYS") &    # Select the cystein pairs
                                                                                                                (df_all["atom_name1"] == "SG") & (df_all["atom_name2"] == "SG") &    # Select the sulphurs atoms pairs
                                                                                                                (df_all["atom_dist"] <= 2.2)]    # Select pairs within the cutoff
    
    # Use an option to select the correct subset of data 
    if intrachain:
        df_disulphide = df_disulphide[(df_disulphide["chain_id1"] == df_disulphide["chain_id2"])]
    elif interchain:
        df_disulphide = df_disulphide[(df_disulphide["chain_id1"] != df_disulphide["chain_id2"])]


    # Check if the dataframe is empty and print accordingly
    if df_disulphide.empty:
        print("")
        print(f"NO {INTRAINTER} DISULPHIDE BRIDGES FOUND\n\n\n".center(106))
    else:
        header_disulphide = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "Distance"]
        table_disulphide = tabulate(df_disulphide, headers = header_disulphide, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
        print(table_disulphide, "\n\n\n")



    #============================= HBOND RESULTS ==============================#
    
    # Get the HBOND output file
    page_results = launching_HBONDS(pdb_file)
    # Parse the HBOND file
    df_hbond = body_to_list(page_results)

    
    legend_hbond = """
Dd-a     =   Distance Between Donor and Acceptor
Dh-a     =   Distance Between Hydrogen and Acceptor
A(d-H-N) =   Angle Between Donor-H-N
A(a-O=C) =   Angle Between Acceptor-O=C
MO       =   Multiple Occupancy
Note that angles that are undefined are written as 999.99


"""

    if df_hbond.empty:
        print("")
        print(f"NO {INTRAINTER} HYDROGEN BONDS FOUND\n\n\n".center(106))
    else:

        # Main Chain-Main Chain Hydrogen Bonds #################################

        # Use an option to select the correct subset of data 
        if intrachain:
            df_hbond = df_hbond[(df_hbond[1] == df_hbond[5])]
        elif interchain:
            df_hbond = df_hbond[(df_hbond[1] != df_hbond[5])]

        print(f"{intrainter} Main Chain-Main Chain Hydrogen Bonds\n".center(106))
        
        # Create a dataframe with main chain-main chain hydrogen bonds
        df_main_main = df_hbond[[0,1,2,3,4,5,6,7,10,11,12,13]][df_hbond[9] == "MM"]
    
        # Check if the dataframe is empty and print accordingly
        if df_main_main.empty:
            print("")
            print(f"NO {INTRAINTER} MAIN CHAIN-MAIN CHAIN HYDROGEN BONDS FOUND\n\n\n".center(106))
        else:
            header_main_main = ["POS", "CHAIN", "RES", "ATOM", "POS", "CHAIN", "RES", "ATOM", "Dd-a", "Dh-h", "A(d-H-N)", "A(a-O=C)"]
            table_main_main = tabulate(df_main_main, headers = header_main_main, showindex=False, numalign="left", tablefmt="rst")
            print("            DONOR                         ACCEPTOR                          PARAMETERS              ")
            print(table_main_main)
            print(legend_hbond)
    

        # Main Chain-Side Chain Hydrogen Bonds #################################
        print("{intrainter} Main Chain-Side Chain Hydrogen Bonds\n".center(106))
    
        # Create a dataframe with main chain-side chain hydrogen bonds
        df_main_side = df_hbond[[0,1,2,3,4,5,6,7,8,10,11,12,13]][(df_hbond[9] == "SO") | (df_hbond[9] == "SN")]
        
        # Check if the dataframe is empty and print accordingly
        if df_main_side.empty:
            print("")
            print(f"NO {INTRAINTER} MAIN CHAIN-SIDE CHAIN HYDROGEN BONDS FOUND\n\n\n".center(106))
        else:
            header_main_side = ["POS", "CHAIN", "RES", "ATOM", "POS", "CHAIN", "RES", "ATOM", "MO", "Dd-a", "Dh-h", "A(d-H-N)", "A(a-O=C)"]
            table_main_side = tabulate(df_main_side, headers = header_main_side, showindex=False, numalign="left", tablefmt="rst")
            print("             DONOR                         ACCEPTOR                             PARAMETERS                ")
            print(table_main_side)
            print(legend_hbond)
    
    
        # Side Chain-Side Chain Hydrogen Bonds #################################
        print(f"{intrainter} Side Chain-Side Chain Hydrogen Bonds\n".center(106))
    
        # Create a dataframe with side chain-side chain hydrogen bonds
        df_side_side = df_hbond[[0,1,2,3,4,5,6,7,8,10,11,12,13]][df_hbond[9] == "SS"]
        
        # Check if the dataframe is empty and print accordingly
        if df_side_side.empty:
            print("")
            print(f"NO {INTRAINTER} SIDE CHAIN-SIDE CHAIN HYDROGEN BONDS FOUND\n\n\n".center(106))
        else:
            header_side_side = ["POS", "CHAIN", "RES", "ATOM", "POS", "CHAIN", "RES", "ATOM", "MO", "Dd-a", "Dh-h", "A(d-H-N)", "A(a-O=C)"]
            table_side_side = tabulate(df_side_side, headers = header_side_side, showindex=False, numalign="left", tablefmt="rst")
            print("             DONOR                         ACCEPTOR                             PARAMETERS                ")
            print(table_side_side)
            print(legend_hbond)



    # Ionic Interactions #######################################################
    print(f"{intrainter} Ionic Interactions\n".center(106))
    print(f"Ionic Interactions within {ionic_cutoff} Angstroms")

    # Create a dataframe with ionic interactions
    df_ionic = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][((((df_all["res_name1"].isin(ani_aa)) & (df_all["res_name2"].isin(cat_aa))) &   # Select the correct residues pairs
                                                                                                 ((df_all["atom_name1"].isin(["OD1", "OD2", "OE1", "OE2"])) &                  # Select correct atoms
                                                                                                  (df_all["atom_name2"].isin(["NH1", "NH2", "NZ", "ND1", "NE2"])))) |
                                                                                                (((df_all["res_name2"].isin(ani_aa)) & (df_all["res_name1"].isin(cat_aa))) &   # Select the correct residues pairs
                                                                                                 ((df_all["atom_name2"].isin(["OD1", "OD2", "OE1", "OE2"])) &                  # Select correct atoms
                                                                                                  (df_all["atom_name1"].isin(["NH1", "NH2", "NZ", "ND1", "NE2"]))))) &
                                                                                                  (df_all["atom_dist"] < ionic_cutoff) & (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates()   # Select atoms in the correct threshold

    # Use an option to select the correct subset of data 
    if intrachain:
        df_ionic = df_ionic[(df_ionic["chain_id1"] == df_ionic["chain_id2"])]
    elif interchain:
        df_ionic = df_ionic[(df_ionic["chain_id1"] != df_ionic["chain_id2"])]

    # Check if the dataframe is empty and print accordingly
    if df_ionic.empty:
        print("")
        print("NO IONIC INTERACTIONS FOUND\n\n\n".center(106))
    else:
        header_ionic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain"]
        table_ionic = tabulate(df_ionic, headers = header_ionic, showindex=False, numalign="left", tablefmt="rst")
        print(table_ionic, "\n\n\n")



    # Aromatic-Aromatic Interactions ###########################################
    print(f"{intrainter} Aromatic-Aromatic Interactions\n".center(106))
    print(f"Aromatic-Aromatic Interactions within {aromarom_cutoff_min} and {aromarom_cutoff_max} Angstroms")

    # Calculate the aromatics centroids coordinates
    res_1, res_2 = [], []
    chain_1, chain_2 = [], []
    list_dist = []
    
    # Check if the Dataframe is empty to not raise an exception
    try:
        # Create 2 lists of residus, 2 lists of chains and one list of distance
        centroids_PHE_TYR = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"].isin(["PHE", "TYR"])) &    # Select the correct amino acids
                                                                            (df_all["atom_name1"].isin(["CG", "CD1", "CD2", "CE1", "CE2", "CZ"])) &    # Select the correct atoms
                                                                            (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()    # Group them by chain_id and res num and calculate the mean coordinate value (centroid)
        centroids_TRP = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"] == "TRP") &    # Select the correct amino acids
                                                                        (df_all["atom_name1"].isin(["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"])) &    # Select the correct atoms
                                                                        (df_all["res_num1"] != df_all["res_num2"])].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()    # Group them by chain_id and res num and calculate the mean coordinate value (centroid)
    
        # Create a dataframe with every aromatics centroids
        centroids_arom = pd.concat([centroids_PHE_TYR, centroids_TRP])
        # Grab the index (chain_id, res_num)
        index_a = list(centroids_arom.index)
        # Transform the dataframe into a coordinate array
        arr_centro = centroids_arom.to_numpy()
        # Calculate the distance matrix between centroids
        distmat_centro = distance_matrix(arr_centro, arr_centro)

        #Select the residues names and chain ID of every pairs of centroids in the correct threshold
        for i in range(len(index_a)):
            for j in range(i+1, len(index_a)):
                if (distmat_centro[i,j] > aromarom_cutoff_min) & (distmat_centro[i,j] < aromarom_cutoff_max):
                    res_1.append(index_a[i][1])
                    res_2.append(index_a[j][1])
                    chain_1.append(index_a[i][0])
                    chain_2.append(index_a[j][0])
                    list_dist.append(distmat_centro[i,j])

    except pd.core.base.DataError:
        pass

    # Create a dataframe with aromatic-aromatic interactions according to the residues corresponding to the residues prviously identified (centroids)
    df_aromatic = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & 
                                                                                               ((df_all["res_num2"] == res_2[i]) & (df_all["chain_id2"] == chain_2[i]))) | 
                                                                                              (((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) &
                                                                                               ((df_all["res_num2"] == res_1[i]) & (df_all["chain_id2"] == chain_1[i])))].drop_duplicates()
        df_aromatic = df_aromatic.append(row)


    if not df_aromatic.empty:
        #Insert the dist column
        df_aromatic["D(centroid-centroid)"] = list_dist

        # Use an option to select the correct subset of data 
        if intrachain:
            df_aromatic = df_aromatic[(df_aromatic["chain_id1"] == df_aromatic["chain_id2"])]
        elif interchain:
            df_aromatic = df_aromatic[(df_aromatic["chain_id1"] != df_aromatic["chain_id2"])]

        # Check if the dataframe is empty and print accordingly
        if df_aromatic.empty:
            print("")
            print(f"NO {INTRAINTER} AROMATIC-AROMATIC INTERACTIONS FOUND\n\n\n".center(106))
        else:

            #Creates a beautiful table with tabulate :)
            header_aromatic = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "D(centroid-centroid)", "Dihedral Angle"]
            table_aromatic = tabulate(df_aromatic, headers = header_aromatic, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
            print(table_aromatic, "\n\n\n")
    else:
        print("")
        print(f"NO {INTRAINTER} AROMATIC-AROMATIC INTERACTIONS FOUND\n\n\n".center(106))



    # Aromatic-Sulphur Interactions ############################################
    print(f"{intrainter} Aromatic-Sulphur Interactions\n".center(106))
    print(f"Aromatic-Sulphur Interactions within {aromsulph_cutoff} Angstroms")

    # Calculates the distance matrix between aromatics centroids and slphur atoms
    res_1, res_2 = [], []
    chain_1, chain_2 = [], []
    list_dist = []

    # Check if the Dataframe is empty to not raise an exception
    try:
        df_coors_s = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"].isin(sulph_aa)) &    # Select the sulphides amino acids
                                                                     (df_all["atom_name1"].str.contains("S"))].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()    # Select only the sulphur atoms
        # Grabs the index (chain_id, res_num)
        index_s = list(df_coors_s.index)
        # Transforms the dataframe into a coordinate array
        arr_coors_s = df_coors_s.to_numpy()

        # Calculates the distance matrix between aromatic centroids and sulphur atoms
        distmat_centro_sulphur = distance_matrix(arr_centro, arr_coors_s)

        #Selects the residues names and chain ID of every pairs of centroids/sulphurs in the correct threshold
        for i in range(len(index_a)):
            for j in range(len(index_s)):
                if distmat_centro_sulphur[i,j] < aromsulph_cutoff:
                    res_1.append(index_a[i][1])
                    res_2.append(index_s[j][1])
                    chain_1.append(index_a[i][0])
                    chain_2.append(index_s[j][0])
                    list_dist.append(distmat_centro_sulphur[i,j])

    except pd.core.base.DataError:
        pass

    # Create a dataframe with aromatic-sulphur interactions
    df_aromatic_sulphur = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & 
                                                                                               ((df_all["res_num2"] == res_2[i]) & (df_all["chain_id2"] == chain_2[i]))) | 
                                                                                              (((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & 
                                                                                               ((df_all["res_num2"] == res_1[i]) & (df_all["chain_id2"] == chain_1[i])))].drop_duplicates()
        df_aromatic_sulphur = df_aromatic_sulphur.append(row)

    # Insert the dist column if the dataframe is not empty
    if not df_aromatic_sulphur.empty:
        df_aromatic_sulphur["D(centroid-centroid)"] = list_dist

        # Use an option to select the correct subset of data 
        if intrachain:
            df_aromatic_sulphur = df_aromatic_sulphur[(df_aromatic_sulphur["chain_id1"] == df_aromatic_sulphur["chain_id2"])]
        elif interchain:
            df_aromatic_sulphur = df_aromatic_sulphur[(df_aromatic_sulphur["chain_id1"] != df_aromatic_sulphur["chain_id2"])]


        # Check if the dataframe is empty after the option and print accordingly
        if df_aromatic_sulphur.empty:
            print("")
            print(f"NO {INTRAINTER} AROMATIC-SULPHUR INTERACTIONS FOUND\n\n\n".center(106))
        else:
            header_aromatic_s = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "D(Centroid-Sulphur)", "Angle"]
            table_aromatic_s = tabulate(df_aromatic_sulphur, headers = header_aromatic_s, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
            print(table_aromatic_s, "\n\n\n")

    else:
        print("")
        print(f"NO {INTRAINTER} AROMATIC-SULPHUR INTERACTIONS FOUND\n\n\n".center(106))



    # Cation-Pi Interactions ###################################################
    print(f"{intrainter} Cation-Pi Interactions\n".center(106))
    print(f"Cation-Pi Interactions within {aromcation_cutoff} Angstroms")

    # Calculates the distance matrix between aromatics centroids and cationics atoms
    res_1, res_2 = [], []
    chain_1, chain_2 = [], []
    list_dist = []

    # Check if the Dataframe is empty to not raise an exception
    try:
        df_coors_i = df_all[["chain_id1","res_num1","x1","y1","z1"]][(df_all["res_name1"].isin(cat_aa))    # Check that the residue is correct
                                                                   & (df_all["atom_name1"].isin(["NH1", "NH2", "NZ"]))].drop_duplicates().groupby(["chain_id1", "res_num1"]).mean()    # Select the correct atoms

        # Extract the residues names and chain id of all cationics amino acids
        index_i = list(df_coors_i.index)
        # Cast the coordinates into a NumPy array
        arr_coors_i = df_coors_i.to_numpy()
        # Calculate the distance matrix
        dist_mat_centro_i = distance_matrix(arr_centro, arr_coors_i)

        # Selects the residues names and chain ID of every pairs of centroids/cations in the correct threshold
        for i in range(len(index_a)):
            for j in range(len(index_i)):
                if dist_mat_centro_i[i,j] < aromcation_cutoff:
                    res_1.append(index_a[i][1])
                    res_2.append(index_i[j][1])
                    chain_1.append(index_a[i][0])
                    chain_2.append(index_i[j][0])
                    list_dist.append(dist_mat_centro_i[i,j])

    except pd.core.base.DataError:
        pass

    # Create a dataframe with cation-pi interactions
    df_arom_i = pd.DataFrame()
    for i in range(len(res_1)):
        row = df_all[["res_num1","res_name1","chain_id1","res_num2","res_name2","chain_id2"]][(((df_all["res_num1"] == res_1[i]) & (df_all["chain_id1"] == chain_1[i])) & 
                                                                                               ((df_all["res_num2"] == res_2[i]) & (df_all["chain_id2"] == chain_2[i]))) | 
                                                                                              (((df_all["res_num1"] == res_2[i]) & (df_all["chain_id1"] == chain_2[i])) & 
                                                                                               ((df_all["res_num2"] == res_1[i]) & (df_all["chain_id2"] == chain_1[i])))].drop_duplicates()
        df_arom_i = df_arom_i.append(row)

    if not df_arom_i.empty:
        #Insert the dist column
        df_arom_i["D(cation-Pi)"] = list_dist

        # Use an option to select the correct subset of data 
        if intrachain:
            df_arom_i = df_arom_i[(df_arom_i["chain_id1"] == df_arom_i["chain_id2"])]
        elif interchain:
            df_arom_i = df_arom_i[(df_arom_i["chain_id1"] != df_arom_i["chain_id2"])]

        # Check if the dataframe is empty and print accordingly
        if df_arom_i.empty:
            print("")
            print(f"NO {INTRAINTER} CATION-PI INTERACTIONS FOUND\n\n\n".center(106))
        else:
            header_arom_i = ["Position", "Residue", "Chain", "Position", "Residue", "Chain", "D(cation-Pi)", "Angle"]
            table_arom_i = tabulate(df_arom_i, headers = header_arom_i, showindex=False, numalign="left", floatfmt=".2f", tablefmt="rst")
            print(table_arom_i, "\n")
    else:
        print("")
        print(f"NO {INTRAINTER} CATION-PI INTERACTIONS FOUND\n\n\n".center(106))
