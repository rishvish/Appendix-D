
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Theta Simulation Study

<!-- badges: start -->

[![Paper
DOI](https://img.shields.io/badge/Paper_DOI-10.1007/s10651--023--00563--w-green)](https://doi.org/10.1007/s10651-023-00563-w)
[![Code
DOI](https://img.shields.io/badge/Code_DOI-10.5281/zenodo.10827269-blue)](https://doi.org/10.5281/zenodo.10827269)
<!-- badges: end -->

<div style="text-align: justify">

This repository contains the source code for reproducing results of the
theta simulation study performed in [Vishwakarma et
al. 2023](https://doi.org/10.1007/s10651-023-00563-w) to test the
robustness of the non-linear parameter $\theta$ in Generalized
Diversity-Interactions (GDI; [Kirwan et
al. 2009](https://doi.org/10.1890/08-1684.1); [Connolly et
al. 2013](https://doi.org/10.1111/1365-2745.12052)) models across the
different interaction structures and suggest a model selection procedure
for fitting GDI models.

### The File/Folder descriptions are as follows:

- <u>**small example for verification.R:**</u> This is a small subset of
  the simulation to give a flavour of the results. Execute this file if
  you don’t want to run the entire simulation to verify the results. The
  code takes anywhere between 20-45 minutes to execute completely.

- <u>**4 species example.R:**</u> Code to run the 4 species example from
  the main paper. Also includes example with other interaction
  structures. Can take up to a couple of days to run the entire code.
  Pre-computed R objects have been provided which can be used to verify
  the results.

- <u>**9 species example.R:**</u> Code to run the 9 species example from
  the main paper. Also includes example with other interaction
  structures. Can take three or four days to run the entire code.
  Pre-computed R objects have been provided which can be used to verify
  the results.

- <u>**6 species example.R:**</u> Code to run the 6 species example from
  the appendix. Can take up to 7 hours to run the entire code.
  Pre-computed R objects have been provided which can be used to verify
  the results.

- <u>**16 species example.R:**</u> Code to run the 16 species example
  from the appendix to test different functional groupings. Can take up
  to 15 hours to run the entire code. Pre-computed R objects have been
  provided which can be used to verify the results.

- <u>**72 species example:**</u> Code to run the 72 species example from
  the appendix. Can take up to 2 days to run the entire code.
  Pre-computed R objects have been provided which can be used to verify
  the results.

- <u>**Datasets folder:**</u> This folder contains a list of folders
  each pertaining to one example code file and holds the datasets
  created by the code as well as the pre-computed R object which can be
  used to verify the results.

- <u>**Helpers folder:**</u>

  - <u>**Final_Pipeline.R:**</u> This file includes the main functions
    for simulating the raw data, testing the robustness of theta and the
    efficacy of the three model selection procedures. The main packages
    used are `readr` for reading and writing .csv files, `tidyverse` for
    data manipulation, `DImodels` for fitting Diversity-Interaction
    models and `ggpubr` for the themes of the plots. If any of these
    packages aren’t present in the system they’ll be installed
    automatically.

  - <u>**Pipeline_Helper_Code.R:**</u> Helper functions used for
    estimating theta and testing efficacy of model selection procedures.
    Used by `Final_Pipeline.R`

  - <u>**Theta_Reparameterization.R:**</u> Code showing examples which
    highlight that the reparameterisation of theta has no effect on its
    estimate.

</div>
