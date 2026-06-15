Code and data for household transmission analyses
This repository contains the R code used to reproduce the figures and analyses for the household transmission study in the article “A multiscale mathematical model integrating viral load improves reconstruction of transmission chains and epidemiological parameters in household studies”
Large simulation data and posterior result files are not stored directly on GitHub because of their size. They are provided separately as an external data archive.
Repository structure
src/
  Core model, path, helper, and visualization functions

Figures/
  Scripts used to generate the manuscript and supplementary figures

data/
  Small metadata files  

results/
  Directory for results files
After downloading the external data archive, the project should contain:
data/simulated_data/
results/chains/
results/ct/
Data and results
The full simulation data and fitted model outputs are available separately:
DOI: 10.5281/zenodo.20702104

Download and extract the archive at the root of this repository so that the folders data/simulated_data/, results/chains/, and results/ct/ are restored.
The cleaned results folder uses the following model names:
symptoms_only
constant
qualitative_vl_no_contact_symptoms
quantitative_vl_no_contact_symptoms
quantitative_vl_with_contact_symptoms
Reproducing the figures
Run the scripts from the repository root:
source("Figures/Figure2.R")
source("Figures/Figure3.R")
source("Figures/Fig4A.R")
source("Figures/Figure4BC.R")
source("Figures/FigureS1.R")
Generated figures are saved in:
Figures/
Requirements
The analyses were run in R. Required packages include:
dplyr
ggplot2
ggpubr
tibble
scales
here
HDInterval
transport
ggrepel
ggforce
purrr
Install missing packages with:
install.packages(c(
  "dplyr", "ggplot2", "ggpubr", "tibble", "scaleadapt ts",
  "here", "HDInterval", "transport", "ggrepel", "ggforce", "purrr"
))
Notes
The GitHub repository contains code and lightweight metadata only. Large .rds files are excluded from version control and should be obtained from the external archive.
