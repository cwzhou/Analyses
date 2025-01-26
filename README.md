# Simulations and Application for itrSurv

This repository contains R code for reproducibility of the simulation and real-world data application analyses in the paper "Optimal individualized treatment regimes for survival data with competing risks", focusing on applying itrSurv to estimate optimal individualized treatment regimes in precision medicine using multi-utility optimization.

The estimator utilizes the R package "itrSurv" from https://github.com/cwzhou/itrSurv, to evaluate survival outcomes jointly with either 1) competing risk endpoints or 2) recurrent events. The analyses are designed to be run in R 4.1.3.

## Table of Contents
1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Requirements](#requirements)
4. [Usage Instructions and Reproducibility](#usage-reproduce)
5. [References](#references)
---

## Overview
This repository aims to ensure transparency and reproducibility for the analyses performed in the project. It includes:
- Scripts for generating and analyzing simulation data.
- Application of statistical methods to real-world data.
- Documentation of key findings.

## Directory Structure

Below is the simplified directory structure for relevant scripts.
```bash
├── RDA /                    # Scripts for real-data application analyses  
│   ├── Paper1_CR/           # Implementation Scripts for Competing Risk Endpoint  
│   │   ├── 3_output/        # RDS Output  
│   │   └── 3_figure/        # Figure Output  
│   ├── Paper2_multiCR/      # N/A for now 
│   └── Paper3_RE/           # Implementation Scripts for Recurrent Event Endpoint 
│   │   ├── output/          # RDS Output  
│   │   └── figure/          # Figure Output 
├── Simulations /            # Scripts for simulating CR/RE survival data and analyses  
│   ├── Paper1_CR/           # Implementation Scripts for Competing Risk Datasets
│   ├── Paper2_multiCR/      # N/A for now
│   └── Paper3_RE/           # Implementation Scripts for Recurrent Event
├── requirements.txt         # List of R package dependencies  
└── README.md                # This file  
```

## Requirements
To run the code, you need the following:
1. **Programming Languages**: R 4.1.3.
2. **Libraries/Packages**: Install dependencies
The following R packages should be installed before running:
install.packages("survival")
install.packages("itrSurv")
...
...
* need to update this list

## Usage Instructions and Reproducibility

This GitHub repository contains two folders: ~/Simulations and ~/RDA, which store reproducible codes for all computational work in the manuscript. Within each folder, Project1_CR is for the competing risk endpoint, Project2_multiCR is N/A for now, and Project3_RE is for the recurrent event endpoint. Please select accordingly.

1. Competing Risks (CR)
   #### Simulations
For cluster computing, run the bash script CR_S2value.sh to submit all jobs implementing each simulation setting on Slurm.

   #### RDA
The PAD observational cohort data is not publically available but the analysis code is provided here.

2. Recurrent Events (RE)
:folder: Simulations

:folder: RDA

## References
1. Zhou, Christina W., Nikki LB Freeman, Katharine L. McGinigle, and Michael R. Kosorok. "Optimal individualized treatment regimes for survival data with competing risks." arXiv preprint arXiv:2411.08315 (2024).
2. Fine-Gray
3. Ghosh-Lin
