#' Test the significance of a target eigengene
#'
#' Performs a permutation test of significance for the eigengene
#'
#' Given a dataset, a set of results from blockwiseModules, and a target module,
#' this function performs a permutation test of the variance explained by the eigengene
#' (first eigenvalue of the correlation matrix of the genes belonging to the module).
#' The target module should be expressed in the same way as when using blockwiseModules
#' (i.e., color or number)
#' Notice that option bicor=TRUE can be extremely computationally demanding
#'
#' @param OriginalData Matrix or data frame
#' containing the original data on which the blockwiseModules function has been run
#' (observations in rows, variables in columns).
#' @param Original_blockwiseModules output of the blockwiseModules function
#' @param target_module module whose largest eigenvalue/eigengene will be tested
#' @param permutations number of permutations to use
#' @param bicor whether one should use the standard eigengene function in WGCNA (default), or bicor
#' @param ... further parameters to be passed to either moduleEigengenes or bicor
#'
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{moduleEigengenes_ALL}{The eigengenes for all the modules (will be empty if bicor=TRUE)}
#'   \item{explained_variance_target_module}{The observed explained variance for the target module}
#'   \item{p_value}{The p value obtained through permutation}
#' }
#'
#' @section Citation:
#' If you use this function please cite Fruciano et al. in press
#'
#' @references Fruciano C., Franchini P., Meyer A. in press. Divergent Allometric Trajectories in Gene Expression and Coexpression Produce Species Differences in Sympatrically SpeciatingMidas Cichlid Fish. Genome Biology and Evolution
#'
#' @export
test_eigengene_explained_variance = function(OriginalData, Original_blockwiseModules, target_module, permutations = 999, bicor = FALSE, 
    ...) {
    
    extract_data_expr_target_module = OriginalData[, which(Original_blockwiseModules$colors == target_module)]
    # Extract only the original data corresponding to the selected target module
    
    
    if (bicor == FALSE) {
        moduleEigengenes_ALL = moduleEigengenes(OriginalData, colors = Original_blockwiseModules$colors, ...)
        target_eigengene_ME = which(colnames(moduleEigengenes_ALL$eigengenes) == paste("ME", target_module, sep = ""))
        observed_explained_variance = moduleEigengenes_ALL$varExplained[target_eigengene_ME]
        # Compute module eigengenes for all the modules in the Original_blockwiseModules and select the explained variance for the
        # target module
    } else {
        eigen_bicor_observed = eigen(bicor(extract_data_expr_target_module))$values
        observed_explained_variance = (eigen_bicor_observed/sum(eigen_bicor_observed))[1]
        moduleEigengenes_ALL = NA
        # use bicor and eigenvalue decomposition to compute the first eigenvalue
    }
    
    
    permuted_indices = lapply(seq_len(permutations), function(X) lapply(seq_len(ncol(extract_data_expr_target_module)), function(z) sample(seq_len(nrow(extract_data_expr_target_module)), 
        nrow(extract_data_expr_target_module), replace = FALSE)))
    permuted_datasets = lapply(permuted_indices, function(x) as.data.frame(matrix(NA, nrow(extract_data_expr_target_module), 
        ncol(extract_data_expr_target_module))))
    # preallocate list
    for (k in seq_len(length(permuted_datasets))) {
        for (j in seq_len(ncol(extract_data_expr_target_module))) {
            permuted_datasets[[k]][, j] = extract_data_expr_target_module[permuted_indices[[k]][[j]], j]
            
        }
    }
    # Generate permuted datasets with the genes in the target module
    
    colors_permutations = rep(1, ncol(permuted_datasets[[1]]))
    # simple vector of ones
    if (bicor == FALSE) {
        varExplained_permuted_datasets = unlist(lapply(permuted_datasets, function(X) moduleEigengenes(X, colors = colors_permutations, 
            ...)$varExplained))
    } else {
        varExplained_permuted_datasets = unlist(lapply(permuted_datasets, function(X) {
            k = eigen(bicor(X))$values
            (k/sum(k))[1]
        }))
    }
    
    names(varExplained_permuted_datasets) = paste("perm_", seq_len(permutations), sep = "")
    p_value_perm = (length(which(varExplained_permuted_datasets >= observed_explained_variance)) + 1)/(permutations + 1)
    # Compute for each permuted dataset variance explained and a p value
    
    Results = list(moduleEigengenes_ALL, explained_variance_target_module = observed_explained_variance, explained_variance_permuted_datasets = varExplained_permuted_datasets, 
        p_value = p_value_perm)
    # create results list to return
    return(Results)
}
