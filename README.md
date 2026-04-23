# Project title

Code and data associated with [brief description of the study, analysis, or manuscript].

## Overview

This repository contains the code required to reproduce the main analyses, figures, and tables from [title of the paper, preprint, or project]. 
All analyses can be reproduced by running the scripts in the `scripts/` directory, with outputs written to `results/`.


## Repository structure

-   `data/raw/`: original input data
-   `R/`: reusable functions, including the main analysis function
-   `scripts/`: analysis scripts to be run in sequence
-   `results/figures/`: final figures
-   `results/tables/`: exported tables
-   `results/objects/`: intermediate R objects

## Data

Original data are stored in `data/raw/`. Data are processed directly in scripts `02_case_studies_fig2.R` and `S1_index_BCI.R`.

-   `supramedioinfra_labelled_rev.xlsx`: consists in Catalan littoral habitat species (188) in the rows and, in the columns, 148 samples for each of the 3 zones of the littoral habitat: infra-, medio-, and supralittoral. More information on the system: <https://doi.org/10.1016%2Fj.ecss.2014.05.031> and <https://doi.org/10.1016/j.ecss.2021.107623> or contact [mariani\@ceab.csic.es](mailto:mariani@ceab.csic.es){.email}.
-   `bci.tree8.rdata`: BCI data from version Jun 07, 2019 is available as a Dryad dataset located at <https://datadryad.org/dataset/doi:10.15146/5xcp-0d46> . Specifically, we downloaded file bci.tree.zip and extracted file `bci.tree8.rdata` in `data/raw/`. This file contains a record for every tree ever recorded in the 2015 BCI census.
-   `aaz4797_ruger_data_s1.xlsx`: PCA data on demographic parameters for 282 BCI species, which corresponds to two trade-offs (growth-survival and stature-recruitment), included in <https://doi.org/10.1126/science.aaz4797>.

## Scripts

The manuscript analyses in `scripts/` are recommended to be run in order, though this is not strictly necessary in most cases:

-   `01_conceptual_figure.R`: generates the cooccurrence-occupancy curve for an *in silico* community of six species, as a first example of this neglected relationship.
-   `02_case_studies_fig2.R`: imports raw data, creates the cooccurrence-occupancy curve and association index for the two case studies in the main text, and plots both.
-   `03_beta_distribution.R`: simulates communities following two distinct Beta distributions for 1000 times, calculates the mean and 90% of the distribution of cooccurrences for each occupancy, and plots together the Beta distribution and these means and 90% distributions.
-   `S1_index_BCI.R`: relates the association index to demographic trade-offs for the BCI data, plotting it. The script requires prior execution of `02_case_studies_fig2.R`, or alternatively loading the saved object `resbci.RData`.

## Requirements

-   R \>= 4.x
-   Required packages: tidyverse, readxl, patchwork.
