# Short Project: Protein Interactions Calculator
[![Geckodriver](https://img.shields.io/badge/Geckodriver-=0.26.0-brightgreen.svg)](https://github.com/mozilla/geckodriver/releases)
[![Selenium](https://img.shields.io/badge/Selenium-brightgreen.svg)](https://selenium-python.readthedocs.io/)
[![Firefox](https://img.shields.io/badge/Firefox-brightgreen.svg)](https://www.mozilla.org/en-US/)
[![Pandas](https://img.shields.io/badge/Pandas-brightgreen.svg)](https://pandas.pydata.org/)
[![ScyPy](https://img.shields.io/badge/SciPy-brightgreen.svg)](https://www.scipy.org/)



M2-BI short-project at Université de Paris.

## Contents 

  * [About this project](#about-this-project)
  * [Authors](#authors)
  * [Requirements](#requirements)
  * [Usage](#usage)

## About this project

The goal of this project is to re-implement calculations from the [Protein Interactions Calculator](http://pic.mbu.iisc.ernet.in/job.html) webserver, allowing to calculate differents interactions in a given protein structure. Our script is written in Python 3, and requier several packages included in a conda environment ([env.yml](env.yml)) for an user-friendly installation. It also requires an internet connection. 

## Authors

  * Debbah Nagi
  * Vander Meersche Yann 

We are students in 2<sup>nd</sup> year of bioinformatics master's Degree, and this work was realised for a short-project. 

## Requirements

To execute the script, you need the folowing Python 3 packages:
  - geckodriver (0.26 version)
  - Pandas
  - NumPy
  - SciPy
  - Firefox
  - Selenium
  - tabulate

You can find the corresponding conda environment file into the [repository](env.yml). 

You also need a protein structure in PDB format that will be used for calculation. 

## Usage 

First, you can get the whole repository by using git clone:

```bash
git clone https://github.com/Kambeiz/Projet-Court/
```

Then, once in the `Projet-Court` directory, you can install the requierd packages by creating and activating the conda environment.

```bash
conda env create -f env.yml
conda activate shortproject
```
Finally, you can launch the script with the following command line to calculate intra- and inter-chain interactions: 

```bash
python3 src/PIC.py -pdb data/1BTA.pdb
```
You can also use the help option *(-h)* to display every command line options. With these options, you can select the type of interaction *(intra- or inter-chain calculation)* and set your own thresholds, just like the PIC webserver.

```bash
usage: PIC.py [-h] -pdb PDBFILE [-intra | -inter] [-hydro HYDROPHOBIC]
              [-ion IONIC] [-AA A A] [-AS AROMSULPH] [-AC AROMCATION]

optional arguments:
  -h, --help            show this help message and exit
  -pdb PDBFILE, --pdbfile PDBFILE
                        Path to the PDB file.
  -intra, --intrachain  Use this option to only print the intra-chain
                        interactions.
  -inter, --interchain  Use this option to only print the inter-chain
                        interactions.
  -hydro HYDROPHOBIC, --hydrophobic HYDROPHOBIC
                        Enter the interaction cut-off value (Default 5A).
  -ion IONIC, --ionic IONIC
                        Enter the interaction cut-off value (Default 6A).
  -AA A A, --aromarom A A
                        Enter the interaction cut-off value (Default 4.5A to
                        7A).
  -AS AROMSULPH, --aromsulph AROMSULPH
                        Enter the interaction cut-off value (Default 5.3A).
  -AC AROMCATION, --aromcation AROMCATION
                        Enter the interaction cut-off value (Default 6A).
```

Which in practice, give you something like :

```bash
python3 src/PIC.py -p data/1BTA.pdb -hydro 6 -ion 5 -AA 5 7 -AS 6 -AC 7 -intra
```

The script will then print the results in your shell. If you want to save the results, don't forget to redirect the command output.
If the script is launched without an internet connection, there will be no results for *Hydrogen Bonds* calculation, and if the script take unusual time to run, PIC might be down and this script can not get the HBOND output file. Please wait, the script will still perform the  remaining calculations.
