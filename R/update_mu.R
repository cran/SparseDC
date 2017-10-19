#'Update the Center Values
#'
#' This fucntion updates the center values for each cluster for each iteration
#' of SparseDC. This function runs inside \code{sparse_dc_fun}
#' @param C_1 The cluster membership for samples in condition 1
#' @param C_2 The cluster membership for samples in condition 2
#' @param pdat1 The centered data from condition 1, columns should be
#' samples (cells) and rows should be features (genes).
#' @param pdat2 The centered data from condition 2, columns should be
#' samples (cells) and rows should be features (genes). The number of genes
#' should be the same as \code{pdat1}.
#' as in pdat1.
#' @param ncluster The number of clusters present in the data.
#' @param lambda1 The lambda 1 value to use in the SparseDC function. This
#' value controls the number of marker genes detected for each of the clusters
#' in the final result. This can be calculated using the "lambda1_calculator"
#' function or supplied by the user.
#' @param lambda2 The lambda 2 value to use in the SparseDC function. This
#' value controls the number of genes that show condition-dependent
#' expression within each cell type. This can be calculated using the
#' "lambda2_calculator" function or supplied by the user.
#' @return Returns a list containing the center values for each of the
#' clusters in condition 1 and 2.
#'
update_mu <- function(C_1, C_2, pdat1, pdat2, ncluster, lambda1, lambda2){
  mu_1 <- matrix(NA, nrow=nrow(pdat1), ncol=ncluster)  # Create Mean vector for data 1
  mu_2 <- matrix(NA, nrow=nrow(pdat2), ncol=ncluster)  # Create Mean vector for data 2
  for (k in 1 : ncluster) {  # For each cluster
    nk_1 <- sum(C_1 == k)  # Number of samples in cluster k for data 1
    nk_2 <- sum(C_2 == k)  # Number of samples in cluster k for data 2
    if (nk_1 > 0 & nk_2 > 0){
      Xk_1 <- pdat1[, C_1 == k, drop=FALSE]  # Create data frame of just samples from cluster k for data 1
      Xk_2 <- pdat2[, C_2 == k, drop=FALSE]  # Create data frame of just samples from cluster k for data 2
      mean_1 <- rowMeans(Xk_1)  # Create the mean vector for cluster k for data 1
      mean_2 <- rowMeans(Xk_2)  # Create the mean vector for cluster k for data 2
      s1 <- S_func(mean_1 - (((sqrt(nk_1) + sqrt(nk_2)) / nk_1) * lambda2), lambda1 / sqrt(nk_1))
      s2 <- S_func(mean_2 + (((sqrt(nk_1) + sqrt(nk_2)) / nk_2) * lambda2), lambda1 / sqrt(nk_2))
      s3 <- S_func(mean_1 + (((sqrt(nk_1) + sqrt(nk_2)) / nk_1) * lambda2), lambda1 / sqrt(nk_1))
      s4 <- S_func(mean_2 - (((sqrt(nk_1) + sqrt(nk_2)) / nk_2) * lambda2), lambda1 / sqrt(nk_2))
      s5 <- S_func(((nk_1 * mean_1) + (nk_2 * mean_2))/ (nk_1 + nk_2),
                   ((sqrt(nk_1) + sqrt(nk_2))/(nk_1 + nk_2)) * lambda1)
      # from clusters 1  & 2 for positive penalty 1
      cond_1 <- (s1 > s2)  # Check if S1 > S.2
      cond_2 <- (s3 < s4)
      mu_1[, k] <- s5
      mu_1[cond_1, k] <- s1[cond_1]
      mu_1[cond_2, k] <- s3[cond_2]
      mu_2[, k] <- s5
      mu_2[cond_1, k] <- s2[cond_1]
      mu_2[cond_2, k] <- s4[cond_2]
    } else if (nk_1 == 0 & nk_2 != 0) {
      Xk_2 <- pdat2[, C_2 == k, drop=FALSE]  # Create data frame of just samples from cluster k for data 2
      mean_2 <- rowMeans(Xk_2)  # Create the mean vector for cluster k for data 2
      s1 <- S_func(mean_2, lambda1 / sqrt(nk_2))
      mu_1[,k] <- mu_2[,k] <- s1
    } else if (nk_1 != 0 & nk_2 == 0){
      Xk_1 <- pdat1[, C_1 == k, drop=FALSE]  # Create data frame of just samples from cluster k for data 1
      mean_1 <- rowMeans(Xk_1)  # Create the mean vector for cluster k for data 1
      s1 <- S_func(mean_1, lambda1 / sqrt(nk_1))
      mu_1[,k] <- mu_2[,k] <- s1
    } else if (nk_1 == 0 & nk_2 == 0){
      mu_1[,k] <- 0
      mu_2[,k] <- 0
    }
  }
  return(list(mu_1=mu_1, mu_2=mu_2))
}
