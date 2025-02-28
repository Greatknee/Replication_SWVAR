# Replication package for "VAR Model with Sparse Group LASSO for Multi-population Mortality Forecasting"

Tim Boonen and Yuhuai Chen
The reproducibility package was assembled at 03/01/2025

## Overview & contents

The code in this replication material generates all numerical experiments figures and tables for
the paper "VAR Model with Sparse Group LASSO for Multi-population Mortality Forecasting".
Each figure and table is generated separately by its corresponding script file
`Figure_[xx]_*.R` or `Table_[xx]_*.R`, respectively.

The main contents of the repository are the following:

- `data/`: folder of raw data files and the functions for processing them
- `out/`: folder of generated plots as PDF files and generated tables as csv files.
- `function/`: folder of all nested functions
- `Figure_[xx]_*.R`: R scripts to create the respective figures
- `Table_[xx]_*.R`: R scripts to create the respective tables


## Instructions & computational requirements.

All file paths are relative to the root of the replication package. Please set your working directory accordingly, or open the `.Rproj` file using RStudio.

The analysis files `Figure_[xx]_*.R` and `Table_[xx]_*.R` can be run individually, in any order.

These analyses were run on R 4.3.1, and we explicitly use the following packages in the analysis files:  
abind(1.4-8), ggplot2(3.5.1), dplyr(1.1.4), psych(2.4.3), sparsevar(0.1.0), vars(1.6-1), MASS(7.3-60.0.1), forecast(8.23.0), vital(1.1.0), demography(2.0), tsibble(1.1.5), dplyr(1.1.4), CVXR(1.0-14).

A comprehensive list of dependencies can be found in the `renv.lock` file. For a convenient setup in a (local) R session, we recommend using the `renv` package. The following steps are required once:
```
# install.packages("renv")
renv::activate()
renv::restore() # install dependencies
renv::status() # check environment
```

## Data availability and provenance
All data used in the paper 
### Human Mortality Database
The data are located at `data`. These files are downloaded from the https://www.mortality.org/.

### Geography Centroidd distance.
The data are located at `data`. These files are manually obtained from https://www.distancefromto.net/.

## Computer configuration and the expected runtime.
| File Name  | Expected runtime |
| ------------- | ------------- |
| Table_01  | 11 min  |
| Content Cell  | 12 min  |

## Special Case
If the STAR model cannot reproduce intermediary data files from the raw data because of time constraints or insufficient CPU, the rest of the verification can still be carried out by adding a parameter `fitSTAR = FALSE` in the analysis function `func_table_[xx]_*()` in the analysis file `Table_[xx]_*.R`.



 
