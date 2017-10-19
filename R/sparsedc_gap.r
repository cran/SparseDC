#' Gap Statistic Calculator
#'
#' This function calculates the gap statistic for SparseDC. For use
#' when the number of clusters in the data is unknown. We recommend
#' using alternate methods to infer the number of clusters in the
#' data.
#'
#' @param pdat1 The centered data from condition 1, columns should be
#' samples (cells) and rows should be features (genes).
#' @param pdat2 The centered data from condition 2, columns should be
#' samples (cells) and rows should be features (genes). The number of genes
#' should be the same as \code{pdat1}.
#' as in pdat1.
#' @param min_clus The minimum number of clusters to try, minimum value is 2.
#' @param max_clus The maximum number of clusters to try.
#' @param nboots The number of bootstrap repetitions to use, default = 200.
#' @param nitter The max number of iterations for each of the start values, the
#' default value is 20.
#' @param nstarts The number of start values to use for SparseDC. The default
#' value is 10.
#' @param l1_boot The number of bootstrap repetitions used for estimating
#' lambda 1.
#' @param l2_boot The number of bootstrap repetitions used for estimating
#' lambda 2.
#' @return A list containing the optimal number of clusters, as well as gap
#' statistics and the calculated standard error for each number of clusters.
#' @examples
#' # load a small dataset
#' data_test <- data_biase[1:50,]
#' # Split data into conditions 1 and 2
#' data_1 <- data_test[ , which(condition_biase == "A")]
#' data_2 <- data_test[ , which(condition_biase == "B")]
#' # Preprocess data (log transform and center)
#' pre_data <- pre_proc_data(data_1, data_2, norm = FALSE, log = TRUE,
#' center = TRUE)
#' # Run with one bootstrap sample for example
#' gap_stat <- sparsedc_gap(pre_data[[1]], pre_data[[2]],
#'  min_clus <- 2, max_clus <- 3, nboots <- 2)
#'
#'
sparsedc_gap <- function(pdat1, pdat2, min_clus, max_clus,
                         nboots = 200, nitter = 20, nstarts = 10,
                         l1_boot = 50, l2_boot = 50) {
  clus_seq <- c(min_clus:(max_clus+1))
  clus_score_vec <- clus_gap_vec <- clus_se_vec <- rep(NA, length(clus_seq))
  score_gap_vec <- score_se <- rep(NA, length(clus_seq))
  for( i in 1:length(clus_seq)){
    l1 <- lambda1_calculator(pdat1, pdat2, ncluster = clus_seq[i],
                             nboot1 = l1_boot)
    l2 <- lambda2_calculator(pdat1, pdat2, ncluster = clus_seq[i],
                             nboot2 = l2_boot)
    fit <- sparsedc_cluster(pdat1, pdat2, ncluster = clus_seq[i], lambda1 = l1,
                            lambda2 = l2)
    dist_total <- 0
    for ( k in 1:clus_seq[i]){
      dist_total <- dist_total  + ((1/2) *
                                     sum((pdat1[,which(fit$clusters1 == k)] -
                                            fit$centers1[,k])^2)) +
        (sqrt(sum(fit$clusters1 == k)) * l1 * sum(abs(fit$centers1[,k]))) +
        ((1/2)*
           sum((pdat2[,which(fit$clusters2 == k)] - fit$centers2[,k])^2))+
        (sqrt(sum(fit$clusters2 == k)) * l1 * sum(abs(fit$centers2[,k]))) +
        ((sqrt(sum(fit$clusters1 == k)) + sqrt(sum(fit$clusters2 == k)))
         * l2 * sum(abs(fit$centers1[,k] - fit$centers2[,k])))
    }
    clus_score_vec[i] <- dist_total
    X.1b <- matrix(rep(NA, nrow(pdat1)*ncol(pdat1)), nrow=nrow(pdat1),
                   ncol=ncol(pdat1))
    X.2b <- matrix(rep(NA, nrow(pdat2)*ncol(pdat2)), nrow=nrow(pdat2),
                   ncol=ncol(pdat2))
    boot_score_vals <- boot_score_vals_s <- rep(NA, nboots)
    for ( j in 1:nboots){
      set.seed(j)
      X.1b <- generate_uni_dat(pdat1)
      X.2b <- generate_uni_dat(pdat2)
      l1_b <- lambda1_calculator(X.1b, X.2b, ncluster = clus_seq[i],
                                 nboot1 = l1_boot)
      l2_b <- lambda2_calculator(X.1b, X.2b, ncluster = clus_seq[i],
                                 nboot2 = l2_boot)
      fit_1 <- sparsedc_cluster(X.1b, X.2b,
                                ncluster = clus_seq[i],
                                lambda1 = l1_b,
                                lambda2 = l2_b)
      dist <- 0
      for ( k in 1:clus_seq[i]){
        dist <- dist + ((1/2) *
                          sum((X.1b[,which(fit_1$clusters1 == k)] -
                                 fit_1$centers1[,k])^2)) +
          (sqrt(sum(fit_1$clusters1 == k)) * l1_b *
             sum(abs(fit_1$centers1[,k]))) + ((1/2)*
                                                sum((X.2b[,which(fit_1$clusters2 == k)] - fit_1$centers2[,k])^2))+
          (sqrt(sum(fit_1$clusters2 == k)) * l1_b *
             sum(abs(fit_1$centers2[,k]))) +
          ((sqrt(sum(fit_1$clusters1 == k)) + sqrt(sum(fit_1$clusters2 == k)))
           * l2_b * sum(abs(fit_1$centers1[,k] - fit_1$centers2[,k])))
      }

      boot_score_vals[j] <- dist
    }
    clus_gap_vec[i] <- (1/nboots)*sum(log(boot_score_vals)) -
      log(clus_score_vec[i])
    clus_se_vec[i] <- stats::sd(log(boot_score_vals)) * sqrt(1 + 1/nboots)

  }
  diff_vec <- (clus_gap_vec ) - clus_se_vec
  test_vec <- rep(NA, length(clus_seq) - 1)
  for(l in 1:(length(clus_seq) - 1)){
    test_vec[l] <- clus_gap_vec[l] > diff_vec[l + 1]
  }

  if (sum(test_vec == TRUE) == 0){
    print("The gap statistic could not find a solution")
  } else {
    gap_clus <- which(test_vec == TRUE)[1]
  }
  gap_clus_res <- clus_seq[gap_clus]
  return(list(k = gap_clus_res, gap_stat = clus_gap_vec, gap_se = clus_se_vec))
}
#'
#' Uniform data generator
#' For use with the gap statistic. Generates datasets drawn from the reference
#' distribution where each reference feature is generated uniformly over the
#' range of observed values for that feature.
#'
#' @param data A dataset with rows as features and columns as samples.
#' @return A dataset drawn from the reference distribution for use internally
#' with the gap statistic.
#'
generate_uni_dat <- function(data){
  min_vals <- apply(data, 1, min)
  max_vals <- apply(data, 1, max)
  uniform_dat <- matrix(stats::runif(length(data), min = min_vals, max = max_vals),
                        ncol= ncol(data))
  return(uniform_dat)
}
