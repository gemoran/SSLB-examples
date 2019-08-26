# SSLB-examples

This repository contains R code to reproduce the results in the paper "Spike-and-Slab Lasso Biclustering". 

It requires the R package `SSLB` (https://github.com/gemoran/SSLB). 

## Index

- `SSLB_functions.R`: helper functions required for all of the R scripts

- `sim_study_1`: folder contains code to reproduce the results in Simulation 1. 

- `sim_study_2`: folder contains code to reproduce the results in Simulation 2.

- `breastCancerNKI`: folder contains code to (i) process data from the `breastCancerNKI` R package for biclustering and 
(ii) run code to reproduce the results in Section 4 of Moran et al (2019). 

- `zeisel`: folder contains code to (i) process data from Zeisel et al (2015) and 
(ii) run code to reproduce the results in Section 5 of Moran et al (2019). 

NOTE: for both `breastCancerNKI` and `zeisel`, the script `0_process_data.R` does not need to be run as processed data is stored in `../data`. The script is included here for completeness.



## References

- Moran G, Rockova V, George E (2019) 
"Spike-and-Slab Lasso Biclustering"

- Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G, Quackenbush J (2019). 
breastCancerNKI: Genexpression dataset published by van't Veer et al. [2002] and van de Vijver et al. [2002] (NKI). R package version 1.22.0, http://compbio.dfci.harvard.edu/.

- Amit Zeisel, Ana B. Muñoz Manchado, Peter Lönnerberg, Gioele La Manno, Simone Codeluppi, Anna Juréus, Sueli Marques, Hermany Munguba, Liqun He, Christer Betsholtz, Charlotte Rolny, Gonçalo Castelo-Branco, Jens Hjerling-Leffler and Sten Linnarsson
"Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq"
*Science* (2015)
