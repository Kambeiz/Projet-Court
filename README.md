# Short Project: Protein Interactions Calculator
[![Geckodriver](https://img.shields.io/badge/Geckodriver-=0.26.0-brightgreen.svg)](https://github.com/mozilla/geckodriver/releases)
[![Selenium](https://img.shields.io/badge/Selenium-brightgreen.svg)](https://selenium-python.readthedocs.io/)
[![Firefox](https://img.shields.io/badge/firefox-brightgreen.svg)](https://www.mozilla.org/en-US/)
[![Pandas](https://img.shields.io/badge/Pandas-brightgreen.svg)](https://pandas.pydata.org/)
[![ScyPy](https://img.shields.io/badge/Pandas-brightgreen.svg)](https://www.scipy.org/)



Projet M2 BI for Université de Paris 

## Contents 

  * [About this project](#about-this-project)
  * [Authors](#authors)
  * [Requirements](#requirements)
  * [Usage](#usage)

## About this project

The goal of this project is to re-implement calculations from the [Protein Interactions Calculator](http://pic.mbu.iisc.ernet.in/job.html) webserver, allowing to calculate differents interactions in a given protein structure. Our script is written in Python 3, and requier several packages included in a conda environment (env.yml) for an user-friendly installation. It also requires an internet connection. 

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

You can find the corresponding conda environment file into the repository (à ajouter un lien direct). 

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
Finally, you can launch the script with the following command line: 

```bash
python3 src/PIC.py -p data/1BTA.pdb
```

