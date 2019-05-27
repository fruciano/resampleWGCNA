#' Analysis of the association between a target module and one or more variables
#'
#'
#' Given a dataset, a set of results from blockwiseModules, a target module, and one or more variables,
#' this function performs a permutation test based on Escoufier RV.
#' Optionally, it can also compute rarefied estimates of RV.
#' Notice that this function requires the package GeometricMorphometricsMix
#' and that the analysis can be computationally demanding for large modules
#'
#' @param OriginalData Matrix or data frame
#' containing the original data on which the blockwiseModules function has been run
#' (observations in rows, variables in columns).
#' @param Original_blockwiseModules output of the blockwiseModules function
#' @param target_module module whose association will be assessed
#' @param X one or more variables whose association with the target module will be tested
#' @param permutations number of permutations to use in the test of association
#' @param rarefied if TRUE, it also performs rarefaction-based estimation of Escoufier RV
#' @param reps_rarefaction number of replicated rarefied samples
#' @param samplesize_rarefaction sample size to which rarefy (required if rarefied=TRUE)
#'
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{ObsRV}{Observed RV coefficient}
#'   \item{p_value_perm}{p value obtained through permutation}
#'   \item{Rarefied_RV}{Results of the rarefaction analysis (empty if rarefied=FALSE)}
#' }
#'
#'
#' @section Citation:
#' If you use this function please cite
#' Fruciano et al. in press (development of the method for WGCNA results) and
#' Escoufier 1973 (for the use of the RV statistic).
#'
#' If you use the option to get rarefied estimates, in addition to the above, please kindly cite
#' Fruciano et al. 2013
#'
#'
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C., Franchini P., Meyer A. 2013. Resampling-Based Approaches to Study Variation in Morphological Modularity. PLoS ONE 8:e69376.
#' @references Fruciano C., Franchini P., Meyer A. in press. Divergent Allometric Trajectories in Gene Expression and Coexpression Produce Species Differences in Sympatrically SpeciatingMidas Cichlid Fish. Genome Biology and Evolution
#'
#'
#' @export
test_RV_target_module = function(OriginalData, Original_blockwiseModules, target_module, X, permutations = 999, rarefied = FALSE,
    reps_rarefaction = 1000, samplesize_rarefaction = NULL) {

    extract_data_expr_target_module = OriginalData[, which(Original_blockwiseModules$colors == target_module)]
    # Extract only the original data corresponding to the selected target module

    if (rarefied == TRUE) {
        if (length(samplesize_rarefaction) == 0) {
            stop(paste("rarefied set to TRUE but no samplesize_rarefaction provided"))
        }
        Rarefied_RV = GeometricMorphometricsMix::RVrarefied(Block1 = extract_data_expr_target_module, Block2 = X, reps = reps_rarefaction,
            samplesize = samplesize_rarefaction)
    } else {
        Rarefied_RV = NA
    }
    # Subroutine to compute rarefied estimates of Escoufier RV

    ObsRV = GeometricMorphometricsMix::EscoufierRV(extract_data_expr_target_module, X)
    # Observed value of RV

    n_obs = nrow(X)
    permuted_indices = lapply(seq_len(permutations), function(x) sample(seq_len(n_obs), n_obs, replace = FALSE))
    permuted_X = lapply(permuted_indices, function(x) cbind(X[x, ]))
    # Create permuted X

    permuted_RV = unlist(lapply(permuted_X, function(x) GeometricMorphometricsMix::EscoufierRV(extract_data_expr_target_module,
        x)))
    p_value_perm = (length(which(permuted_RV >= ObsRV)) + 1)/(permutations + 1)
    # Compute permuted RV and corresponding p value

    Results = list(Observed_RV = ObsRV, p_value_perm = p_value_perm, Rarefied_RV = Rarefied_RV)
    return(Results)
}
