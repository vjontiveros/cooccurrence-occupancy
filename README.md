# Code and data for: *The Co‐occurrence–Occupancy Curve in Ecological Communities: a novel approach to species association*.

This repository contains code and data to reproduce analyses for a community ecology study on the co-occurrence–occupancy relationship. The project examines how the tendency of species to co-occur with other species depends on their occupancy across sites, and develops a framework to distinguish non-random species association from patterns expected simply because some species are more common than others.

The code defines and analyses empirical co-occurrence–occupancy curves, derives expected curves under null models assuming site equivalence and species independence, and calculates a Species Association Index (SAI): an occupancy-standardized measure of whether species co-occur with others more or less often than expected given their frequency of occurrence.

The analyses are illustrated using two ecological case studies: tropical forest trees from Barro Colorado Island, and organisms from Mediterranean rocky shores.

## Overview

This repository contains the code required to reproduce the main analyses, figures, and tables from *The Co‐occurrence–Occupancy Curve in Ecological Communities: a novel approach to species association*. All analyses can be reproduced by running the scripts in the `scripts/` directory, with outputs written to `results/`.

## Repository structure

-   `data/raw/`: original input data
-   `R/`: reusable functions, including the main analysis function
-   `scripts/`: analysis scripts to be run in sequence
-   `results/figures/`: final figures
-   `results/objects/`: intermediate R objects

## Data

Original data are stored in `data/raw/`, except `aaz4797_ruger_data_s1.xlsx` (see below). Data are processed directly in scripts `02_case_studies_fig2.R` and `S1_index_BCI.R`.

-   `supramedioinfra_labelled_rev.xlsx`: consists in Catalan littoral habitat species (188) in the rows and, in the columns, 148 samples for each of the 3 zones of the littoral habitat: infra-, medio-, and supralittoral. More information on the system: <https://doi.org/10.1016%2Fj.ecss.2014.05.031> and <https://doi.org/10.1016/j.ecss.2021.107623> or contact [mariani\@ceab.csic.es](mailto:mariani@ceab.csic.es){.email}.
-   `bci.tree8.rdata`: BCI data from version Jun 07, 2019 is available as a Dryad dataset located at <https://datadryad.org/dataset/doi:10.15146/5xcp-0d46> . Specifically, we downloaded file bci.tree.zip and extracted file `bci.tree8.rdata` in `data/raw/`. This file contains a record for every tree ever recorded in the 2015 BCI census.
-   `aaz4797_ruger_data_s1.xlsx`: PCA data on demographic parameters for 282 BCI species, which corresponds to two trade-offs (growth-survival and stature-recruitment), included in <https://doi.org/10.1126/science.aaz4797>. This file is not redistributed in this repository; users should download it from the article supplementary materials and place it in `data/raw/`.

## Scripts

The manuscript analyses in `scripts/` are recommended to be run in order, though this is not strictly necessary in most cases:

-   `01_conceptual_figure.R`: generates the cooccurrence-occupancy curve for an *in silico* community of six species, as a first example of this neglected relationship. This script uses helper functions from `R/m_values.R`. Outputs figure `fig1_conceptual.png`.
-   `02_case_studies_fig2.R`: imports raw data, creates the cooccurrence-occupancy curve and association index for the two case studies in the main text, and plots both. This script uses helper functions from `R/m_values.R`. Outputs figure `fig2_case_studies.png`.
-   `03_beta_distribution.R`: simulates communities following two distinct Beta distributions for 1000 times, calculates the mean and 90% of the distribution of cooccurrences for each occupancy, and plots together the Beta distribution and these means and 90% distributions. This script uses helper functions from `R/m_values.R`. Outputs figure `fig3_betas.png`.
-   `04_empirical_occupancy.R`: compares the observed cooccurrence–occupancy relationship against expectations derived from (i) the empirical occupancy distribution and (ii) a fitted log-series occupancy distribution for the case studies. The script requires prior execution of `02_case_studies_fig2.R`, or alternatively loading the saved objects `resbci.RData` or `resmedio.RData`. This script uses helper functions from `R/expected_A_empirical.R` and `R/logseries_functions.R`. Outputs figure `fig4_empirical.png`.
-   `S1_index_BCI.R`: relates the association index to demographic trade-offs for the BCI data, plotting it. The script requires prior execution of `02_case_studies_fig2.R`, or alternatively loading the saved object `resbci.RData`. Outputs figure `fig_s1_tradeoffs.png`.
-   `S2_curve_for_beta.R`: compares the observed cooccurrence–occupancy relationship against expectations derived from (i) the empirical occupancy distribution assuming species equivalence and (ii) a fitted beta occupancy distribution for the case studies. It also produces a figure of the fit of the beta distribution to the data. The script requires prior execution of `02_case_studies_fig2.R`, or alternatively loading the saved objects `resbci.RData` or `resmedio.RData`. This script uses helper functions from `R/expected_A_sp_equivalence.R` and `R/beta_functions.R`. NOTE: This script requires a prior python installation. Outputs figures `fig_s2.png` and `fig_s3.png`.

## Requirements

-   R \>= 4.x
-   Required packages: tidyverse, readxl, patchwork.
-   Python \>= 3.x and package reticulate in R (needed to run script `S2_curve_for_beta.R`).


### License

The code in this repository is licensed under the MIT License. See `LICENSE` for details.

External data files are subject to their original licenses or reuse conditions. In particular, `aaz4797_Ruger_Data_S1.xlsx` is supplementary material from Rüger et al. (2020), Science, DOI: 10.1126/science.aaz4797, and is not covered by the MIT License of this repository.