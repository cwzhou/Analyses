# Simulations and Application for itrSurv:CR

This repository contains R code for reproducibility of the simulation and real-world data application analyses in the paper "Optimal individualized treatment regimes for survival data with competing risks", focusing on simulating survival data with competing risks as well as applying itrSurv for survival and competing risk outcomes to estimate optimal individualized treatment regimes in precision medicine.


The estimator utilizes the R package "itrSurv" from https://github.com/cwzhou/itrSurv, to evaluate the competing risk endpoint. The analyses are designed to be run in R 4.1.3.

## Table of Contents
1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Requirements](#requirements)
4. [Usage Instructions](#usage-instructions)
   - [Simulations](#simulations)
   - [Applications](#applications)
5. [Reproducibility](#reproducibility)
6. [Contact Information](#contact-information)
---

## Overview
This repository aims to ensure transparency and reproducibility for the analyses performed in the project. It includes:
- Scripts for generating and analyzing simulation data.
- Application of statistical methods to real-world data.
- Documentation of key findings.

## Directory Structure


## Requirements
To run the code, you need the following:
1. **Programming Languages**: R 4.1.3.
2. **Libraries/Packages**: Install dependencies




## Summary

This GitHub repository contains two folders: ~/Simulations and ~/RDA, which store reproducible codes for all computational work in the manuscript.

1. Simulations
The following R packages should be installed before running:
install.packages("survival")
install.packages("itrSurv")
...
...
* need to update this list

Use bash script S2value.sh to submit all cluster jobs implementing each simulation setting on Slurm.

2. PAD data analysis
The cohort data is not publically available but the analysis code is provided here.

## File Descriptions

## References
