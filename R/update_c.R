#' Update Clusters
#'
#' Updates the cluster membership for each iteration of SparseDC. Runs inside
#' \code{sparse_dc_fun}.
#' @param mu_1 The center values for each cluster in condition 1.
#' @param mu_2 The center values for each cluster in condition 2.
#' @param pdat1 The centered data from condition 1, columns should be
#' samples (cells) and rows should be features (genes).
#' @param pdat2 The centered data from condition 2, columns should be
#' samples (cells) and rows should be features (genes). The number of genes
#' should be the same as \code{pdat1}.
#' as in pdat1.
#' @param ncluster The number of clusters present in the data.
#' @return A list containing the cluster membership for condition 1 and
#' condition 2.
#'
update_c <- function(mu_1, mu_2, pdat1, pdat2, ncluster){
  # Create distance matrices p * k, number of variables times the number of clusters
  dist_mat_1 <- matrix(NA, nrow=ncol(pdat1), ncol=ncluster)
  dist_mat_2 <- matrix(NA, nrow=ncol(pdat2), ncol=ncluster)
  for (k in 1 : ncluster) {  # For each cluster do:
    dist_mat_1[, k] <- colSums((pdat1 - mu_1[, k]) ^ 2)  # Calculate distance between each sample and the cluster centers in data 1
    dist_mat_2[, k] <- colSums((pdat2 - mu_2[, k]) ^ 2)  # Calcualte distance between each sample and the cluster centers in data 2
  }
  C_1 <- apply(dist_mat_1, 1, which.min)  # Assign samples to cluster whose center they are closest too for data 1
  C_2 <- apply(dist_mat_2, 1, which.min)  # Assign samples to cluster whose center they are closest too for data 2
  return(list(C_1=C_1, C_2=C_2))  # Return new cluster assignments
}
