Household transmission model: code and reproducibility materials
This repository contains the R code used to reproduce the analyses and figures for a household transmission modelling study “A multiscale mathematical model integrating viral load improves reconstruction of transmission chains and epidemiological parameters in household studies”.
The repository is designed to contain code and lightweight metadata only. Large simulation datasets and fitted model outputs are provided separately through Zenodo.
Repository structure
.
├── src/
│   ├── main.R
│   ├── paths.R 
│   ├── helpers.R
│   ├── models/
│   └── visualization/
│
├── Figures/
│   ├── Figure2.R
│   ├── Figure3.R
│   ├── Fig4A.R
│   ├── Figure4BC.R
│   └── FigureS1.R
│
├── data/
│   └── Small metadata files  
│
├── results/
│   └── Directory for results files
│
└── README.md
Large data and fitted results
Large .rds files are not stored on GitHub.
The full simulation data and fitted model outputs are available on Zenodo:
Zenodo DOI: 10.5281/zenodo.20702104
After downloading and extracting the Zenodo archive, the repository should contain:
data/simulated_data/
results/chains/
results/ct/
The results/ directory should follow this structure:
results/
├── chains/
│   ├── symptoms_only/
│   ├── constant/
│   ├── qualitative_vl_no_contact_symptoms/
│   ├── quantitative_vl_no_contact_symptoms/
│   └── quantitative_vl_with_contact_symptoms/
│
└── ct/
    ├── symptoms_only/
    ├── constant/
    ├── qualitative_vl_no_contact_symptoms/
    ├── quantitative_vl_no_contact_symptoms/
    └── quantitative_vl_with_contact_symptoms/
Model names
The cleaned result system uses the following model names:
Folder name	Description
symptoms_only	Clinical/symptom-only model
constant	Constant transmission model
qualitative_vl_no_contact_symptoms	Qualitative viral-load model, with symptom information removed for infected contacts
quantitative_vl_no_contact_symptoms	Quantitative viral-load model, with symptom information removed for infected contacts
quantitative_vl_with_contact_symptoms	Quantitative viral-load model, with symptom information retained for infected contacts; used for supplementary analyses
Reproducing the figures
Run all scripts from the repository root.
source("Figures/Figure2.R")
source("Figures/Figure3.R")
source("Figures/Fig4A.R")
source("Figures/Figure4BC.R")
source("Figures/FigureS1.R")
Generated figures are saved in:
Figures/
Running inference
The cleaned inference script writes outputs directly into the new result structure:
source("src/Inference_French_HH_26082025_new_results_system.R")
It saves posterior samples and reconstructed transmission chains as:
results/chains/<model>/Households_...rds
results/ct/<model>/Households_...rds
R packages
The analyses require R and the following packages:
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
rstan
Install missing packages with:
install.packages(c(
  "dplyr",
  "ggplot2",
  "ggpubr",
  "tibble",
  "scales",
  "here",
  "HDInterval",
  "transport",
  "ggrepel",
  "ggforce",
  "purrr",
  "rstan"
))
GitHub and Zenodo usage
This GitHub repository contains:
code
figure scripts
helper functions
manifests
documentation
It does not contain:
large simulation data
posterior chains
reconstructed transmission chains
large .rds files
Those files should be obtained from the Zenodo archive and extracted into the repository root.
Notes
All paths are handled using the here package. Scripts should therefore be run from the repository root, or from an RStudio project opened at the repository root.

 
