# Projet-Court : PIC
[![Geckodriver](https://img.shields.io/badge/snakemake-=0.26.0-brightgreen.svg)](https://github.com/mozilla/geckodriver/releases)
[![Firefox](https://img.shields.io/badge/firefox-brightgreen.svg)](https://www.mozilla.org/en-US/)
[![Pandas](https://img.shields.io/badge/firefox-brightgreen.svg)](https://pandas.pydata.org/)
[![Selenium]((https://img.shields.io/badge/firefox-brightgreen.svg)](https://selenium-python.readthedocs.io/)



Projet M2 BI for Université de Paris 

## Contents 

  * [About this project](#about-this-project)
  * [Authors](#authors)
  * [Requirements](#requirements)
  * [Usage](#usage)

## About this project

The goal of this project is to reproduce calculations from the [Protein Interactions Calculator](http://pic.mbu.iisc.ernet.in/job.html) allowing to infer interactions in a given protein. It is written in python language, and include a conda environment for making it user-friendly. 

## Authors

  * Debbah Nagi
  * Vander-Meersche Yann 

We are students in a 2nd year master's Degree of bioinformatics, and this work enter in a short time window project. 

## Requirements

To execute the script, you need those packages:
  - geckodriver (the 0.26 version)
  - pandas
  - python
  - scipy
  - firefox
  - selenium
  - tabulate

You can find the environment file corresponding into the repository (à ajouter un lien direct). 

You also need a pdb file that will be used for calculation. 

## Usage 

Firstly, you can get the whole repository by using git clone:

```bash
git clone https://github.com/Kambeiz/Projet-Court/
```

In your shell, you can launch the script by following this line of command, using one of pdb files provided into our data folder: 

```bash
cd Projet-Court && conda activate shortproject && python 3 src/PIC.py -p data/1BTA.pdb && conda deactivate
```

