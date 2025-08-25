#  Modelling approaches for meta-analyses with dependent effect sizes in ecology and evolution

This repository contains the scripts, data, and materials associated with the manuscript *"Modelling approaches for meta-analyses with dependent effect sizes in ecology and evolution: A simulation study"*.  

## Simulation study overview

We conducted three sub-studies based on different model specifications to assess alternative modelling approaches:  

- **Study 1**: Meta-analysis models  
- **Study 2**: Phylogenetic meta-analysis models  
- **Study 2.sub**: Phylogenetic meta-regression models  

All simulations were run on the UNSW High Performance Computing (HPC) cluster **Katana** (DOI: [10.26190/669X-A286](https://doi.org/10.26190/669X-A286)), which runs on RPM-based Linux OS (Rocky 8). The R version used was **4.3.1**.

## Repository structure

- **`/pbs/`**  
  PBS job scripts for running the simulations on HPC cluster.

- **`/logs/`**  
  HPC error and log files for all jobs.

- **`/R/`**  
  Simulation scripts and supporting R code.  
  Includes scripts to generate simulation data, fit models, and process results.

- **`/output/`**  
  Results of the simulation studies, including:  
  - Simulation parameter files  
  - Figures  

- **`/docs/`**  
  Quarto (`.qmd`) source code for the case study webpage.  
  The rendered case study HTML file is also included.

- **`/data/`**  
  Data used in the case studies, drawn from two published studies.

## Citation

> *add citation*    
> DOI: [Zenodo DOI]

