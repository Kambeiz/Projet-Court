# Projet-Court : PIC

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

In your shell, you can launch the script by followings this line of command: 

```bash
cd Projet-Court && conda activate shortproject && python 3 src/PIC.py -p data/1BTA.pdb && conda deactivate
```
