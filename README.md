# Replication package for "VAR Model with Sparse Group LASSO for Multi-population Mortality Forecasting"

Tim Boonen and Yuhuai Chen

## Overview & contents

The code in this replication material generates all numerical experiments figures and tables for
the paper "VAR Model with Sparse Group LASSO for Multi-population Mortality Forecasting".
Each figure and table is generated separately by its corresponding script file
`Figure_[xx]_*.R` or `Table_[xx]_*.R`, respectively.

The main contents of the repository are the following:

- `plots/`: folder of generated plots as PDF files
- `tables/`: folder of generated tables as txt files
- `data/`: folder of raw data files and the functions for processing them
- `out/`: folder of outpur
- `function/`: folder of all nested functions
- `Figure_[xx]_*.R`: R scripts to create the respective figures
- `Table_[xx]_*.R`: R scripts to create the respective tables


## Instructions & computational requirements.

All file paths are relative to the root of the replication package. Please set your working directory accordingly, or open the `.Rproj` file using RStudio.

The analysis files `Figure_[xx]_*.R` and `Table_[xx]_*.R` can be run individually, in any order.

These analyses were run on R 4.3.1, and we explicitly use the following packages in the analysis files: `triptych` (0.1.2), `ggplot2` (3.4.3), `patchwork` (1.1.3), `dplyr` (1.1.3), `tidyr` (1.3.0), `purrr` (1.0.2), `grid` (base R), `lubridate` (1.9.2).

A comprehensive list of dependencies can be found in the `renv.lock` file. For a convenient setup in a (local) R session, we recommend using the `renv` package. The following steps are required once:
```
# install.packages("renv")
renv::activate()
renv::restore() # install dependencies
renv::status() # check environment
```

## Data availability and provenance

### Human Mortality Database
The data are located at `data`. These files are downloaded from the https://www.mortality.org/.

### Geography Centroidd distance.
The data are located at `data`. These files are manually obtained from https://www.distancefromto.net/.

## Computer configuration and the expected runtime.

| File Name  | Expected runtime |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |
