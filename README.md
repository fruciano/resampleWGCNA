# resampleWGCNA - Resampling-based approaches to weighted gene coexpression analysis results

This package contains two functions which were developed as part an [empirical study](https://doi.org/10.1093/gbe/evz108).
If you use any of the functions in the package, kindly cite:

Fruciano C., Franchini P., Meyer A. *in press*. [Divergent Allometric Trajectories in Gene Expression and Coexpression Produce Species Differences in Sympatrically SpeciatingMidas Cichlid Fish](https://doi.org/10.1093/gbe/evz108). *Genome Biology and Evolution*  11(6):1644–1657.  [doi: 10.1093/gbe/evz108](https://doi.org/10.1093/gbe/evz108)

More specific citation details are provided in the help of each function.

The package currently contains two functions, which perform permutation tests ancillary to a normal WGCNA analysis:

- **test_eigengene_explained_variance**: performs a permutation test of the variance explained by the "eigengene" of a selected module
- **test_RV_target_module**: performs a permutation test based on the Escoufier RV coefficient of the multivariate association between the expression of genes in a given module and a univariate or multivariate array (e.g., phenotype)
