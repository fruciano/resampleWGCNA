# resampleWGCNA - Resampling-based approaches to weighted gene coexpression analysis results

This package contains two functions which are part of a still unpublished manuscript.
Please, kindly contact me ( carmelo.fruciano [at] ens.fr ) if you plan on using it or publishing the results.

It currently contains two functions, which perform permutation tests ancillary to a normal WGCNA analysis:

-**test_eigengene_explained_variance**: performs a permutation test of the variance explained by the "eigengene" of a selected module
- **test_RV_target_module**: performs a permutation test based on the Escoufier RV coefficient of the multivariate association between the expression of genes in a given module and a univariate or multivariate array (e.g., phenotype)
