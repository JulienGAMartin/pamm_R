# Power analysis for random effects in mixed models


### Description:
Simulation functions to assess or explore the power of
a dataset to estimates significant random effects in a mixed model.
The functions are based on the "lme4" package.

### Installation
The package is available on CRAN, in case of problem you can download the developper version directly from github using 
`
devtools::install_github("JulienGAMartin/pamm_R")
`


### References:
Martin, Nussey, Wilson and Reale Submitted Measuring between-individual
variation in reaction norms in field and experimental studies: a power
analysis of random regression models. Methods in Ecology and Evolution.

### Note 2020-01-20
Code as been fixed to work with the latest version of lme4 and lmerTest packages
Suggestiosn for new functionality are welcome
Code need to be cleaned a bit.

### ToDo
- restructure the package based on a simulate function and a power.eval function that PAMM, SSF and EAMM can call
- rewrite support information
