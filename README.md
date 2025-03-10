# Simulations and Application for itrSurv

This repository contains R code for reproducibility of the simulation and real-world data application analyses in the paper "Optimal individualized treatment regimes for survival data with competing risks", focusing on applying itrSurv to estimate optimal individualized treatment regimes in precision medicine using multi-utility optimization.

The estimator utilizes the R package **`itrSurv`** to evaluate survival outcomes jointly with either 1) competing risk endpoints or 2) recurrent events. Please install `itrSurv` according to the README.md from **https://github.com/cwzhou/itrSurv** prior to running any of the scripts in this repo. The analyses are designed to be run in R 4.1.3.

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
**Please install `itrSurv` according to the README.md from https://github.com/cwzhou/itrSurv prior to running any of the scripts in this repo.**
The following R packages should be installed before running:
library(survival)
library(itrSurv)
...
...
* need to update this list

## Usage Instructions and Reproducibility

This repository contains two main folders: `~/Simulations` and `~/RDA`. These folders include the necessary code and scripts for reproducing all computational work presented in the manuscript. Below are the instructions for using the repository and ensuring reproducibility of the analyses. Within each folder, Project1_CR is for the competing risk endpoint, Project2_multiCR is N/A for now, and Project3_RE is for the recurrent event endpoint. Please select accordingly.

### 1. Competing Risks (CR)
#### Simulations

The `~/Simulations/Paper1_CR` folder contains the simulation scripts required to generate the results from the manuscript. There are two simulation settings: 1) generating independent failure tiems from an exponential distribution, and 2) simulating dependent failures times based on the simulation settings from Fine-Gray [2]. The simulation setting is specified by the parameter `generate_failure_method` (see below). To run the simulations:

1. **Edit Simulation Parameters**  
   Open the script `CR00.Simulation_Parameters.R` and modify the following parameters as necessary:
   - `generate_failure_method`: Set to either "simple_exp" or "fine_gray".
   - `local`: Set to `1` for running on a local machine or `0` for running on a cluster.
   - `parallel`: Set to `1` to enable parallel processing using the `parallel` package in R, or set to `0` for non-parallel execution.
   - Other parameters (e.g., `n.sim` for the number of simulations): Modify according to the simulation needs for various simulation scenarios for the for the specified `generate_failure_method` setting.

2. **Cluster Computing**  
   If running on a cluster, submit the jobs using the provided bash script:
   - **`CR_S2value.sh`**: This script submits all simulation jobs to a Slurm scheduler. It will automatically handle the execution of each simulation scenario for the specified `generate_failure_method` setting.

3. **Running Specific Simulations**  
   To run specific simulations manually, use the script **`CR01.Simulation_Run.R`**, which runs individual simulation settings that can be customized as needed for the specified `generate_failure_method` setting.
   
4. **Plotting Results**
   To plot the results from the simulation studies, use the script **`CR02.Simulation_Summary.R`**, which plots the results read from the sepcified date folder, using parametesr from `CR00.Simulation_Parameters.R` for the specified `generate_failure_method` setting. NOTE: this script works as is only if results from more than 2 simulation scenarios exists (since we use `facet_grid` in R). 

#### RDA (Analysis Code)

The `~/RDA/Paper1_CR` folder contains the analysis code for the PAD observational cohort study. Although the data from the PAD cohort is not publicly available, the code provided here allows users to reproduce the analysis. Follow these steps:

1. **Ensure Access to Data**  
   Since the PAD data is not publicly available, please ensure you have access to the dataset before running the analysis code.

2. **Run the Analysis Code**  
   The scripts in this folder will perform the analysis, using the cohort data and the parameters specified in the simulation settings. Refer to the code and modify any data paths or necessary parameters as needed for your environment. Specifically, run **`CR04.PAD_analysis.R`**.

3. **Plotting Results**
   Use the script **`CR05.PAD_summary.R`** to plot the results.

#### Script Descriptions
Below are descriptions of the scripts used but not directly called upon.

### 2. multi CR
Ignore for now.

### 3. Recurrent Events (RE)
#### Simulations
#### RDA
We analyze the public bladder recurrence dataset in the R package 'survival'.

#### Script Descriptions
Below are descriptions of the scripts used but not directly called upon.

---

### General Notes on Reproducibility

- All code is designed to be reproducible on both local machines and computing clusters.
- To ensure full reproducibility, please ensure that all required R packages and dependencies are installed. You can install the required packages by running the following in your R environment:

For Competing Risks scripts:
    ```r
    install.packages(c("parallel", "dplyr", "ggplot2", "itrSurv")) # need to update this list
    ```
For Recurrent Event scripts:
    ```r
    install.packages(c("parallel", "dplyr", "ggplot2", "itrSurv")) # need to update this list
    ```
    
- **Version Control**: If you are using a different version of R or any dependencies, please document them in your environment for better reproducibility.

## References
1. Zhou, C. W., Freeman, N. L., McGinigle, K. L., & Kosorok, M. R. (2024). Optimal individualized treatment regimes for survival data with competing risks. arXiv preprint arXiv:2411.08315.
2. Fine, J. P., & Gray, R. J. (1999). A Proportional Hazards Model for the Subdistribution of a Competing Risk. Journal of the American Statistical Association, 94(446), 496–509. https://doi.org/10.1080/01621459.1999.10474144
3. Ghosh, D., & Lin, D. Y. (2000). Nonparametric analysis of recurrent events and death. Biometrics, 56(2), 554–562. https://doi.org/10.1111/j.0006-341x.2000.00554.x
