# SSLB-examples

This repository contains R code to reproduce the results in the paper "Spike-and-Slab Lasso Biclustering" (Moran, Rockova and George, 2020).

The `SSLB` R package is available at the Github repository: 
`https://github.com/gemoran/SSLB`

Please follow the following instructions to run the code. 

1. Install the `SSLB` R package using the package `devtools` by running the following code in your R GUI:

```
install.packages("devtools")
library(devtools)
install_github("gemoran/SSLB")
```

2. Install the following R packages used in simulations and data analysis:

```
install.packages(c("mvtnorm", "isa2", "biclust", "pals", "reshape2", "tidyverse", "gridExtra", "clue"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("breastCancerNKI", "preprocessCore", "impute", "fabia", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "org.Mm.eg.db"))
```

3. We compare SSLB to the method BicMix (Gao et al. 2016).  We downloaded the BicMix software from [here](https://www.cs.princeton.edu/~bee/software/BicMix-Code-for-distribution.zip). Follow the instructions in the BicMix software to create the BicMix executable file. Change `bicmix_directory.txt` to reflect where the executable file is relative to the folder `SSLB_examples`. The R code will then call the BicMix exectuable file from within the R session.

4. We also compare SSLB to the method SSBiEM (Denitto et al 2017). We downloaded the SSBiEM software from [here](https://github.com/emme-di/SSBiEM/). To run SSBiEM, place the files `kronSpeye.m` and `SSBiEM.m` in each of the folders: `sim_study_1`, `sim_study_1`, `sim_study_1`, `sim_study_1`, `breastCancerNKI` and `zeisel`. Then, to run SSBiEM, run the Matlab code in `[folder-name]_ssbi.m`.


## Index

- `SSLB_functions.R`: helper functions required for all of the R scripts

- `K.R`: this code creates the tables of number of biclusters estimated by each method for all of the simulation studies.

- `sim_study_1`: folder contains code to reproduce the results in Simulation 1. 
    + `sim_study_1_ssbi.m`: this code runs the method SSBiEM
    + `sim_study_1.R`: this code runs all other biclustering methods. 
    + `plots.R`: this code produces the consensus, relevance, recovery and variance explained plots
    + `plot_matrices.R`: this code plots the factor and loadings matrices found by each method for one dataset from the simulation study.

- `sim_study_2`: folder contains code to reproduce the results in Simulation 2. Code structure the same as `sim_study_1`.

- `sim_study_3`: folder contains code to reproduce the results in Simulation 3. Code structure the same as `sim_study_1`.

- `sim_study_4`: folder contains code to reproduce the results in Simulation 4. Code structure the same as `sim_study_1`.

- `breastCancerNKI`: folder contains code to (i) process data from the `breastCancerNKI` R package for biclustering and 
(ii) run code to reproduce the results in Section 4 of Moran et al (2020). 

- `zeisel`: folder contains code to (i) process data from Zeisel et al (2015) and 
(ii) run code to reproduce the results in Section 5 of Moran et al (2020). 

NOTE: for both `breastCancerNKI` and `zeisel`, the script `0_process_data.R` does not need to be run as processed data is stored in `../data`. The script is included here for completeness.



## References

- Moran G, Rockova V, George E (2020) 
"Spike-and-Slab Lasso Biclustering" *Annals of Applied Statistics* (Accepted)

- Denitto, M., Bicego, M., Farinelli, A. and Figueiredo, M. A. (2017). "Spike and slab biclustering." *Pattern Recognition*

- Gao, C., McDowell, I. C., Zhao, S., Brown, C. D. and Engelhardt, B. E. (2016).
"Context Specific and Differential Gene Co-expression Networks via Bayesian Biclustering." *PLoS Comput Biol*.

- Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G, Quackenbush J (2019). 
breastCancerNKI: Genexpression dataset published by van't Veer et al. [2002] and van de Vijver et al. [2002] (NKI). R package version 1.22.0, http://compbio.dfci.harvard.edu/.

- Amit Zeisel, Ana B. Muñoz Manchado, Peter Lönnerberg, Gioele La Manno, Simone Codeluppi, Anna Juréus, Sueli Marques, Hermany Munguba, Liqun He, Christer Betsholtz, Charlotte Rolny, Gonçalo Castelo-Branco, Jens Hjerling-Leffler and Sten Linnarsson (2015)
"Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq" *Science* 
