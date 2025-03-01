# Replication package for "VAR Model with Sparse Group LASSO for Multi-population Mortality Forecasting"

Tim J. Boonen and Yuhuai Chen (The University of Hong Kong).

The reproducibility package was assembled at March 1, 2025.

## Overview & contents

The code in this replication material generates all numerical experiments, figures, and tables for the paper "VAR Model with Sparse Group LASSO for Multi-population Mortality Forecasting".
Each figure and table is generated separately by its corresponding script file `Figure_[xx]_*.R` or `Table_[xx]_*.R`, respectively.

The main contents of the repository are the following:

- `data/`: folder of raw data files and functions for processing them
- `out/`: folder of generated plots as PDF files and generated tables as csv files.
- `function/`: folder of all nested functions
- `param/`: folder of parameters of the SWVAR model
- `Online Resource/`: folder of all files for the online resource 1.
- `Figure_[xx]_*.R`: R scripts for creating the corresponding figures
- `Table_[xx]_*.R`: R scripts for creating the corresponding tables


## Instructions & computational requirements.

All file paths are relative to the root of the replication package. Please set your working directory accordingly, or open the `.Rproj` file using RStudio.

The analysis files `Figure_[xx]_*.R` and `Table_[xx]_*.R` can be run individually, in any order.

These analyses were run on R 4.3.3, and we explicitly use the following packages in the analysis files:  
abind(1.4-8), ggplot2(3.5.1), dplyr(1.1.4), psych(2.4.3), sparsevar(0.1.0), vars(1.6-1), MASS(7.3-60.0.1), forecast(8.23.0), vital(1.1.0), demography(2.0), tsibble(1.1.5), CVXR(1.0-14), glasso(1.11).

A comprehensive list of dependencies can be found in the `renv.lock` file. For a convenient setup in a (local) R session, we recommend using the `renv` package. The following steps are required once:
```
# install.packages("renv")
renv::activate()
renv::restore() # install dependencies
renv::status() # check environment
```

## Data availability and provenance
### Human Mortality Database
These data are downloaded from https://www.mortality.org/ and located at `data`. See Human Mortality Database. (2023, https://www.mortality.org/File/GetDocument/Public/Docs/MethodsProtocolV6.pdf) for more details.

### Geography centroid distance.
These data are manually obtained from https://www.distancefromto.net/ and located at `data`. 

## Files expected runtime.
| File Name  | Expected runtime |
| ------------- | ------------- |
| Table_01  | 11 min  |
| Table_04  | 12 min  |
| Table_07  | 12 min (Setting 1) + 120 min (Setting 2) |
| Figure_0304  | 15 min  |
| Figure_06  | 5 min  |
| Figure_07  | 30 min  |
All unlisted files will cost less than 3 minutes.

## Special case
If the STAR model cannot report or gives an error or a warning because of time constraints or insufficient CPU, the rest of the verification can still be carried out by adding a parameter `fitSTAR = FALSE` in the analysis function `func_table_[xx]_*()` in the analysis file `Table_07_*.R` and `Figure_07_*.R`.

## Reference
Human Mortality Database. (2023). Methods protocol for the Human Mortality Database (Version 6). Retrieved from https://www.mortality.org/File/GetDocument/Public/Docs/MethodsProtocolV6.pdf

 
