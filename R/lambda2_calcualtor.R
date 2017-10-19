#' Lambda 2 Calculator.
#'
#' Calculates the lambda 2 values for use in the main SparseDC algorithm, the
#' lambda 2 value controls the number of genes that show condition-dependent
#' expression within each cell type. That is it controls the number of
#' different mean values across the conditions for each cluster. It is
#' calculated by estimating the value of lambda 2 that would result in no
#' difference in mean values across conditions when there are no meaningful
#' differences across between the conditions. For further details please see
#' the original manuscript.
#' @param pdat1 The centered data from condition 1, columns should be
#' samples (cells) and rows should be features (genes).
#' @param pdat2 The centered data from condition 2, columns should be
#' samples (cells) and rows should be features (genes). The number of genes
#' should be the same as \code{pdat1}.
#' as in pdat1.
#' @param ncluster The number of clusters present in the data.
#' @param nboot2 The number of bootstrap repetitions for estimating lambda 2,
#' the default value is 1000.
#' @param alpha2 The quantile of the bootstrapped lambda 2 values to use,
#' range is (0,1). The default value is 0.5, the median of the calculated
#' lambda 2 values.
#' @return The calculated value of lambda 2 to use in the main SparseDC
#' algorithm.
#' @examples
#'
#' set.seed(10)
#' # Select small dataset for example
#' data_test <- data_biase[1:100,]
#' # Split data into conditions A and B
#' data_A <- data_test[ , which(condition_biase == "A")]
#' data_B <- data_test[ , which(condition_biase == "B")]
#' # Pre-process the data
#' pre_data <- pre_proc_data(data_A, data_B, norm = FALSE, log = TRUE,
#' center = TRUE)
#' # Calculate the lambda 2 value
#' lambda2_calculator(pdat1 = pre_data[[1]], pdat2 = pre_data[[2]], ncluster = 3,
#'  alpha2 = 0.5, nboot2 = 1000)
#'
#'  # Can also run
#'  pdata_A <- pre_data[[1]]
#'  pdata_B <- pre_data[[2]]
#' lambda2_calculator(pdat1 = pdata_A, pdat2 = pdata_B, ncluster = 3,
#'  alpha2 = 0.5, nboot2 = 1000)
#'
#' @seealso  \code{\link{lambda1_calculator}}  \code{\link{sparsedc_cluster}}
lambda2_calculator <- function(pdat1, pdat2, ncluster, alpha2 = 0.5,
                               nboot2 = 1000) {
  pooled_dat <- cbind(pdat1, pdat2)  # Create pooled data
  l2_boot_vals <- rep(NA, nboot2)  #  Create vector to store bootstrap L2 values
  for  (b in 1:nboot2) {
    set.seed(b)
    X_1_p <- pooled_dat[, sample(ncol(pooled_dat), ncol(pdat1), replace=T)]  # Sample pooled data for X.1
    X_2_p <- pooled_dat[, sample(ncol(pooled_dat), ncol(pdat2), replace=T)]  # Sample pooled data for X.2
    clusters1 <- sample(ncluster, ncol(pdat1), replace=T) # Extract Cluster Assignments for X.1
    clusters2 <- sample(ncluster, ncol(pdat2), replace=T) # Extract Cluster Assignments for X.2
    l2_vals <- matrix(rep(NA,ncluster*nrow(X_1_p)),ncol=ncluster)
    for (k in 1:ncluster) {
      nk_1 <- sum(clusters1 == k) # Number of samples in cluster k for data 1
      nk_2 <- sum(clusters2 == k) # Number of samples in cluster k for data 2
      if (nk_1 > 0 & nk_2 > 0){
        Xk_1 <- X_1_p[, clusters1 == k, drop=FALSE]  # Create data frame of just samples from cluster k for data 1
        Xk_2 <- X_2_p[, clusters2 == k, drop=FALSE]  # Create data frame of just samples from cluster k for data 2
        mean_1 <- rowMeans(Xk_1)  # Create the mean vector for cluster k for data 1
        mean_2 <- rowMeans(Xk_2)  # Create the mean vector for cluster k for data 2
        m_1 <- mean_1
        m_2 <- mean_2
        v_1 <- (mean_1 - mean_2) / ((((sqrt(nk_1) + sqrt(nk_2))) / nk_1) +
                                      (((sqrt(nk_1) + sqrt(nk_2)))/ nk_2))
        u_1 <- mean_1 - (((sqrt(nk_1) + sqrt(nk_2))/ nk_1) * v_1)
        u_2 <- mean_2 + (((sqrt(nk_1) + sqrt(nk_2))/ nk_2) * v_1)
        v_2 <- (mean_2 - mean_1) / ((((sqrt(nk_1) + sqrt(nk_2))) / nk_1) +
                                      (((sqrt(nk_1) + sqrt(nk_2)))/ nk_2))
        u_3 <- mean_1 + (((sqrt(nk_1) + sqrt(nk_2)) / nk_1) * v_2)
        u_4 <- mean_2 - (((sqrt(nk_1) + sqrt(nk_2)) / nk_2) * v_2)
        #### Condtions for Solution Assignment ####
        cond_1 <- (u_1 > u_2)
        cond_2 <- (u_3 < u_4)
        #### Assign Solutions ####
        # mu.1 == mu.2
        l2_vals[,k] <- 0
        # mu.1 > mu.2
        l2_vals[cond_1,k] <- v_1[cond_1]
        l2_vals[cond_2,k] <- v_2[cond_2]
      } else if (nk_1 == 0 & nk_2 > 0) {
        l2_vals[,k] <- 0
      } else if (nk_1 > 0 & nk_2 == 0) {
        l2_vals[,k] <- 0
      } else if ( nk_1 == 0  & nk_2 == 0) {
        l2_vals[,k] <- 0
      }
    }
    l2_boot_vals[b] <- max(abs(l2_vals))
  }
  l2 <- unname(stats::quantile(abs(l2_boot_vals), alpha2))
  return(l2)
}
