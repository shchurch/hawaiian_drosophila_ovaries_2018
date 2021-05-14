# hawaiian_drosophila_ovaries_2018
Code associated with the mansucript **Reproductive Capacity Evolves in Response to Ecology through Common Changes in Cell Number in Hawaiian _Drosophila_**, published in Current Biology: https://doi.org/10.1016/j.cub.2019.04.063 and BioRxiv: https://doi.org/10.1101/470898

# Contents

## Comparative analyses

`analyses.r`

This file contains the code required to reproduce the tables in our upcoming manuscript. The code performs a statistical comparison of ecological models, and tests the relationship between reproductive traits with a Phylogenetic Least Squares Regression (PGLS).

`plots.r`

This file contains the code required to reproduce the figures in our upcoming manuscript. The code performs an ancestral state reconstruction and plots the result of a PGLS regression of reproductive traits.

`data.txt` `avgs.txt` `logs.txt`

These files are the raw ovariole number, body size, and egg size data used in the comparative anlaysis. `data.txt` is required as input for `analyses.r`, while the averages `avgs.txt` and log transformed averages `logs.txt` are output.

## Ecological classification of taxa
`ecological_classification.txt`

This file contains the coding scheme used to classify the taxa in our study, as recorded in *Magnacca, Karl N., David Foote, and Patrick M. O’Grady. _A review of the endemic Hawaiian Drosophilidae and their host plants._ Zootaxa 1728.1 (2008): 58.*

## sequence_data_methods/
This directory contains all of the genetic sequence data, code, and resultant gene trees used to perform the species delineation analysis. It also contains the final matrices used in the phylgoenetic analysis with BEAST.

## Phylogenetic results
The following files are the results of our phylogenetic analysis of Hawaiian _Drosophila_. All are nexus format.

`nuc_mt_beast_posterior.trees`

`nuc_mt_beast_mcc.tre`

These are the posterior distribution and maximum clade credibility trees, estimated using both nuclear and mitochondrial markers.

`mt_beast_posterior.trees`

`mt_beast_mcc.tre`

These are the posterior distribution and maximum clade credibility trees, estimated using only mitochondrial markers.

# R Package versions

These are the package versions used to perform all analyses in our upcoming manuscript

```
R version 3.4.2 (2017-09-28)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Sierra 10.12.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] phytools_0.6-60 maps_3.3.0      dplyr_0.7.6     nlme_3.1-137   
 [5] surface_0.4-1   geiger_2.0.6    MASS_7.3-50     ouch_2.11-1    
 [9] subplex_1.5-4   OUwie_1.50      lattice_0.20-35 corHMM_1.22    
[13] GenSA_1.1.7     nloptr_1.0.4    ape_5.2        

loaded via a namespace (and not attached):
 [1] gmp_0.5-13.2            Rcpp_0.12.18            pillar_1.3.0           
 [4] bindr_0.1.1             compiler_3.4.2          digest_0.6.17          
 [7] tibble_1.4.2            rlang_0.2.2             pkgconfig_2.0.2        
[10] Matrix_1.2-14           fastmatch_1.1-0         igraph_1.2.2           
[13] parallel_3.4.2          mvtnorm_1.0-8           expm_0.999-3           
[16] bindrcpp_0.2.2          coda_0.19-1             Rmpfr_0.7-1            
[19] tidyselect_0.2.4        combinat_0.0-8          grid_3.4.2             
[22] nnet_7.3-12             scatterplot3d_0.3-41    glue_1.3.0             
[25] paleotree_3.1.0         deSolve_1.21            R6_2.2.2               
[28] plotrix_3.7-3           animation_2.5           phangorn_2.4.0         
[31] purrr_0.2.5             corpcor_1.6.9           magrittr_1.5           
[34] assertthat_0.2.0        mnormt_1.5-5            numDeriv_2016.8-1      
[37] quadprog_1.5-5          crayon_1.3.4            clusterGeneration_1.3.4
```
