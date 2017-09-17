#' Sparse Differential Clustering
#'
#' The main SparseDC function. This function clustered the samples from the two
#' conditions and links the clusters across the conditions. It also identifies
#' marker genes for each of the clusters. There are three types of marker gene
#' which SparseDC identifies. Please see the original manuscript for further
#' details.
#' @param pdat1 The centered data from condition 1, columns should be
#' samples (cells) and rows should be features (genes).
#' @param pdat2 The centered data from condition 2, columns should be
#' samples (cells) and rows should be features (genes). The number of genes
#' should be the same as \code{pdat1}.
#' as in pdat1.
#' @param ncluster The number of clusters present in the data.
#' @param lambda1 The lambda 1 value to use in the SparseDC function. This
#' value controls the number of marker genes detected for each of the clusters
#' in the final result. This can be calcualted using the "l1_calculator"
#' function or supplied by the user.
#' @param lambda2 The lambda 2 value to use in the SparseDC function. This
#' value controls the number of genes that show condition-dependent
#' expression within each cell type.
#' @param nitter The max number of iterations for each of the start values, the
#' default value is 20.
#' @param nstarts The number of start values to use for SparseDC. The default
#' value is 50.
#' @return A list containing the clustering solution, cluster centers  and the
#' score of each of the starts.
#' @examples
#'
#' set.seed(10)
#' # Select small dataset for example
#' data_test <- data_biase[1:100,]
#' # Split data into conditions 1 and 2
#' data_1 <- data_test[ , which(condition_biase == "A")]
#' data_2 <- data_test[ , which(condition_biase == "B")]
#' # Preprocess data (log transform and center)
#' pre_data <- pre_proc_data(data_1, data_2, norm = FALSE, log = TRUE,
#' center = TRUE)
#' # Calculate lambda 1 parameter
#' lambda1 <- lambda1_calculator(pdat1 = pre_data[[1]], pdat2 = pre_data[[2]],
#'  ncluster=3, alpha1 = 0.5, nboot1 = 1000)
#' # Calculate lambda 2 parameter
#' lambda2 <- lambda2_calculator(pdat1 = pre_data[[1]], pdat2 = pre_data[[2]],
#'  ncluster = 3, alpha2 = 0.5, nboot2 = 1000)
#' # Run sparse DC
#' sdc_res <- sparsedc_cluster(pdat1 = pre_data[[1]], pdat2 = pre_data[[2]], ncluster = 3,
#' lambda1 = lambda1, lambda2 = lambda2, nitter = 20, nstarts =50)
#' # Extract results
#' clusters_1 <- sdc_res$clusters1  # Clusters for condition 1 data
#' clusters_2 <- sdc_res$clusters2  # Clusters for condition 2 data
#' centers_1 <- sdc_res$centers1  # Centers for condition 1 data
#' centers_2 <- sdc_res$centers2  # Centers for condition 2 data
#' # View clusters
#' summary(as.factor(clusters_1))
#' summary(as.factor(clusters_2))
#' # View Marker genes
#' gene_names <- row.names(data_test)
#' m_gene_c1_up1 <- gene_names[which(centers_1[,1] > 0)]
#' m_gene_c1_up2 <- gene_names[which(centers_2[,1] > 0)]
#' m_gene_c1_down1 <- gene_names[which(centers_1[,1] < 0)]
#' m_gene_c1_down2 <- gene_names[which(centers_2[,1] < 0)]
#' m_gene_c2_cond <- gene_names[which(centers_1[,2] != centers_2[,2])]
#'
#' # Can also run
#'
#' pre_data <- pre_proc_data(data_1, data_2, norm = FALSE, log = TRUE,
#' center = TRUE)
#' pdata_A <- pre_data[[1]]
#' pdata_B <- pre_data[[2]]
#' lambda1 <- lambda1_calculator(pdat1 = pdata_A , pdat2 = pdata_B,
#' ncluster=3, alpha1 = 0.5, nboot1 = 1000)
#' lambda2 <- lambda2_calculator(pdat1 = pdata_A, pdat2 = pdata_B,
#'  ncluster = 3, alpha2 = 0.5, nboot2 = 1000)
#' # Run sparse DC
#' sdc_res <- sparsedc_cluster(pdat1 = pdata_A, pdat2 = pdata_B, ncluster = 3,
#' lambda1 = lambda1, lambda2 = lambda2, nitter = 20, nstarts =50)
#'
#' @seealso  \code{\link{lambda1_calculator}}  \code{\link{lambda2_calculator}}
#' \code{\link{update_c}} \code{\link{update_mu}}
sparsedc_cluster <- function(pdat1, pdat2, ncluster, lambda1, lambda2,
                             nitter = 20, nstarts = 50) {
  clusterHistory1star <- clusterHistory2star <- vector(nitter, mode="list")  # Create best cluster assignment history list
  centerHistory1star <- centerHistory2star <- vector(nitter, mode="list")  # Create best center history list
  iter_star <- NA
  score_vec <- rep(NA, nstarts)  # Create vector to contain scores for each of the starts
  for( s in 1:nstarts) {  # For each start
    clusterHistory1 <- clusterHistory2 <- vector(nitter, mode="list")  # Create cluster assignment history list
    centerHistory1 <- centerHistory2 <- vector(nitter, mode="list")  # Create center history list
    c_dat <- cbind(pdat1, pdat2)
    fit_temp <- stats::kmeans(t(c_dat), ncluster, nstart=1)
    mu_1 <- mu_2 <- t(fit_temp$centers)
    centers <- list(mu_1=mu_1,mu_2=mu_2)  # Convert centers to list form
    iter_num <- 0
    clus_same <- FALSE
    while(iter_num < nitter & clus_same == FALSE) {
      iter_num <- iter_num + 1
      clusters <- update_c(mu_1=centers[[1]], mu_2=centers[[2]],
                           pdat1, pdat2, ncluster)  # Update the cluster assignments
      centers <- update_mu(C_1=clusters[[1]], C_2=clusters[[2]],
                           pdat1, pdat2, ncluster, lambda1, lambda2)  # Update centers
      clusterHistory1[[iter_num]] <- clusters[[1]]  # Record clusters for data1
      centerHistory1[[iter_num]] <- centers[[1]]  # Record centers for data1
      clusterHistory2[[iter_num]] <- clusters[[2]]  # Record clusters for data2
      centerHistory2[[iter_num]] <- centers[[2]]  # Record clusters for data2
      if ( iter_num > 1 ) {
        if (all(c(clusterHistory1[[iter_num]],
                  clusterHistory2[[iter_num]]) ==
                c(clusterHistory1[[iter_num-1]],
                  clusterHistory2[[iter_num-1]]))) {
          clus.same <- TRUE
        }
      }
    }
    mu_1 <- centers[[1]]
    mu_2 <- centers[[2]]
    dist_mat_1 <- matrix(NA, nrow=ncol(pdat1), ncol=ncluster)  # Create matrix of distance between each of the points and their cluster centers
    dist_mat_2 <- matrix(NA, nrow=ncol(pdat2), ncol=ncluster)
    mu_diff <- rep(NA,ncluster) # Create vector containg the penalities from the differences between cluster centers
    dist_penalty <- 0  # Create variable to store the distance between each cluster center and the elements contained in the cluster
    nk1 <- nk2 <- rep(NA,ncluster)
    for (k in 1 : ncluster) {   # For each cluster do:
      nk1[k] <- sum((clusterHistory1[[iter_num]] == k))
      nk2[k] <- sum((clusterHistory2[[iter_num]] == k))
      dist_mat_1[, k] <- colSums((pdat1 - mu_1[, k]) ^ 2) # Calculate distance between each sample and the cluster centers in data 1
      dist_mat_2[, k] <- colSums((pdat2 - mu_2[, k]) ^ 2) # Calcualte distance between each sample and the cluster centers in data 2
      mu_diff[k] <- lambda2*(sum(abs(mu_1[,k]-mu_2[,k])) *
                               ((sqrt(nk1[k]) + sqrt(nk2[k]))))  # Calculate the distance between cluster centers
      dist_penalty <- dist_penalty +
        sum(dist_mat_1[which(clusterHistory1[[iter_num]]== k),k]) +
        sum(dist_mat_2[which(clusterHistory2[[iter_num]]== k),k]) # Sum the distance between each cluster center and it's elements
    }
    mu_penalty <- lambda1*(sum(abs(mu_1) %*% diag(sqrt(nk1)))) +
      lambda1*(sum(abs(mu_2) %*% diag(sqrt(nk2))))  # Calculate penalties on the size of the mean vector
    score_vec[s] <- mu_penalty + sum(mu_diff) + dist_penalty  # Calculate the total score for the clustering
    if (score_vec[s] <= min(stats::na.omit(score_vec))) { # If the latest score is the lowest score
      clusterHistory1star <- clusterHistory1  # Replace best cluster history with new cluster history for data 1
      clusterHistory2star <- clusterHistory2  # Replace best cluster history with new cluster history for data 2
      centerHistory1star <- centerHistory1  # Replace best center history with new center history for data 1
      centerHistory2star <- centerHistory2  # Replace best center history with new center history for data 2
      iter_star <- iter_num
    }
  }
  clusters_1 <- unlist(clusterHistory1star[iter_star])
  clusters_2 <- unlist(clusterHistory2star[iter_star])
  centers_1 <- matrix(unlist(centerHistory1star[iter_star]),ncol=ncluster)
  centers_2 <- matrix(unlist(centerHistory2star[iter_star]),ncol=ncluster)

  return(list(clusters1=clusters_1, centers1=centers_1,
              clusters2=clusters_2, centers2=centers_2,
              scores=score_vec))  # Return lists of cluster assignments and centers
}
