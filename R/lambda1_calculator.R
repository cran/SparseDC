#' Lambda 1 Calculator.
#'
#' Calculates the lambda 1 value for the sparseDC algorithm. The lambda 1
#' value controls the number of marker genes selected for each
#' cluster in the output from SparseDC. It is calculated as the value of
#' lambda 1 that results in no marker genes being selected when then are no
#' meaningful clusters present in the data. Please see the original manuscript
#' for further details.
#' @param pdat1 The centered data from condition 1, columns should be
#' samples (cells) and rows should be features (genes).
#' @param pdat2 The centered data from condition 2, columns should be
#' samples (cells) and rows should be features (genes). The number of genes
#' should be the same as \code{pdat1}.
#' as in pdat1.
#' @param ncluster The number of clusters present in the data.
#' @param nboot1 The number of bootstrap repetitions used for estimating lambda 1,
#' the default value is 1000.
#' @param alpha1 The quantile of the bootstrapped lambda 1 values to use,
#' range is (0,1). The default value is 0.5, the median of the calculated
#' lambda 1 values.
#' @return The calculated value of lambda 1 to use in the main SparseDC
#' algorithm.
#' @examples
#' set.seed(10)
#' # Select small dataset for example
#' data_test <- data_biase[1:100,]
#' # Split data into conditions A and B
#' data_A <- data_test[ , which(condition_biase == "A")]
#' data_B <- data_test[ , which(condition_biase == "B")]
#' # Pre-process the data
#' pre_data <- pre_proc_data(data_A, data_B, norm = FALSE, log = TRUE,
#' center = TRUE)
#' # Calculate the lambda 1 value
#' lambda1_calculator(pdat1 = pre_data[[1]], pdat2 = pre_data[[2]], ncluster=3,
#'  alpha1 = 0.5, nboot1 = 1000)
#'
#' # Can also run
#'
#' # Pre-process the data
#' pre_data <- pre_proc_data(data_A, data_B, norm = FALSE, log = TRUE,
#' center = TRUE)
#' pdata_A <- pre_data[[1]]
#' pdata_B <- pre_data[[2]]
#' # Calculate the lambda 1 value
#' lambda1_calculator(pdat1 = pdata_A, pdat2 = pdata_B , ncluster=3,
#'  alpha1 = 0.5, nboot1 = 1000)
#'
#' @seealso  \code{\link{lambda2_calculator}} for how to calculate the lambda 2
#' parameter. \code{\link{sparsedc_cluster}} for the main sparse differential
#' clustering function.
lambda1_calculator <- function(pdat1, pdat2, ncluster, alpha1 = 0.5,
                               nboot1 = 1000) {
  cX_1 <- cbind(pdat1, pdat2)
  n <- ncol(cX_1)  # Calculate number of samples
  l1_bootvals <- rep(NA, nboot1)  # Create vector to store boot strap values
  for (b in 1:nboot1) {
    set.seed(b)
    clus_b <- sample(ncluster, n, replace=T)  # Assign random cluster assignments
    mean_1 <- matrix(NA, nrow=nrow(cX_1), ncol=ncluster)  # Create center vector
    for (k in 1 : ncluster) { # For each cluster
      nk_1 <- sum(clus_b == k)  # Number of samples in cluster k for data 1
      if (nk_1 > 0) {
        Xk_1_b <- cX_1[, clus_b == k, drop=FALSE]  # Create data frame of just samples from cluster k for data 1
        mean_1[,k] <- rowMeans(Xk_1_b) * sqrt(nk_1)  # Create the mean vector for cluster k for data 1
      } else if (nk_1 == 0) {
        mean_1[,k] <- 0
      }
    }
    l1_bootvals[b] <- max(abs(mean_1))
  }
  l1 <- unname(stats::quantile(l1_bootvals, alpha1))
  return(l1)
}
