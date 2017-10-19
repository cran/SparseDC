#' Data Simulator
#'
#' Simulates two condition data for a range of conditions depending on the
#' parameters used.
#'
#' @param genes The number of genes to be simulated.
#' @param cells The number of cells to be simulated per condition.
#' @param sig.genes The number of marker genes for each cluster.
#' @param sig.genes.s The number of marker genes shared across conditions for
#' each cluster. Should be less than or equal to \code{sig.genes}
#' @param clus.t1 A vector of clusters present in the first condition. Start
#' at 1, e.g. \code{c(1,2,3,4)}
#' @param clus.t2 A vector of clusters present in the second condition. Does
#' not have to match clus.t1, e.g \code{c(3,4,5)}
#' @param same.sig TRUE or FALSE. Should each cluster have a unique set of
#' marker genes. default is FALSE.
#' @param u.l Lower bound for the cluster gene means, default is 1.
#' @param u.h Upper bound for the cluster gene means, default is 2.
#'
#' @return A list containing the two simulated data matrices, \code{dat.1} and
#' \code{dat.2}, true clusters for the cells in the first and second
#' conditions, \code{clusters1} and \code{clusters2}, a matrix indicating
#' marker genes for the first and second condition, \code{sig.gene.mat.1}
#' and \code{sig.gene.mat.2}, the base mean values for each gene
#' \code{gene.means} and the cluster specific additions for each gene
#' \code{clus.gene.means}
#'
#' @examples
#' set.seed(10)
#' genes <- 1000  # Simulate 1,000 genes
#' cells <- 100  # Simulate 100 cells per condition
#' clus.t1 <- c(1,2,3)  # Generate 3 clusters present in condition A
#' clus.t2 <- c(1,2,3)  # Generate 3 clusters present in condition B
#' sig.genes <- 30  # Generate 30 marker genes per cluster
#' sig.genes.s <- 15  # Let half of the 30 marker genes be shared.
#' temp_sim_dat <- sim_data(genes, cells, sig.genes, sig.genes.s,
#' clus.t1, clus.t2)
#'
sim_data <- function(genes, cells, sig.genes, sig.genes.s, clus.t1, clus.t2,
                     same.sig = FALSE, u.l = 1 , u.h = 2){
  gene.means <- stats::rnorm(genes, mean=0, sd=1)
  nclus <- length(unique(c(clus.t1, clus.t2)))
  clus.gene.means <- matrix(stats::runif(nclus*genes,u.l,u.h)*
                              sign(stats::rnorm(nclus*genes)),ncol=nclus)

  clusters1 <- sample(clus.t1, cells, replace=T)
  while(min(unname(summary(as.factor(clusters1)))) < 10 |
        length(unique(clusters1)) != length(unique(clus.t1))){
    clusters1 <- sample(clus.t1, cells, replace=T)
  }
  clusters2 <- sample(clus.t2, cells, replace=T)
  while(min(unname(summary(as.factor(clusters2)))) < 10 |
        length(unique(clusters2)) != length(unique(clus.t2))){
    clusters2 <- sample(clus.t2, cells, replace=T)
  }
  sig.gene.mat.1 <- sig.gene.mat.2 <- matrix(rep(0,nclus*genes), ncol=nclus)
  if (same.sig){
    s.sig.genes <- sample(genes, sig.genes.s)
    sig.gene.mat.1[s.sig.genes,] <- sig.gene.mat.2[s.sig.genes, ] <- 1
    if( sig.genes > sig.genes.s){
      sig.genes1 <- sample(setdiff(1:genes,s.sig.genes),sig.genes -
                             sig.genes.s)
      sig.genes2 <- sample(setdiff(1:genes,c(s.sig.genes,sig.genes1)),sig.genes
                           - sig.genes.s)
      sig.gene.mat.1[sig.genes1,] <- 1
      sig.gene.mat.2[sig.genes2,] <- 1
    }
  } else {
    for ( k in 1:nclus){
      s.sig.genes <- sample(genes, sig.genes.s)
      sig.gene.mat.1[s.sig.genes,k] <- sig.gene.mat.2[s.sig.genes, k] <- 1
      if( sig.genes > sig.genes.s){
        sig.genes1 <- sample(setdiff(1:genes,s.sig.genes),
                             sig.genes - sig.genes.s)
        sig.genes2 <- sample(setdiff(1:genes,c(s.sig.genes,sig.genes1))
                             ,sig.genes - sig.genes.s)
        sig.gene.mat.1[sig.genes1, k] <- 1
        sig.gene.mat.2[sig.genes2, k] <- 1
      }
    }
  }
  sig.gene.mat.o <- sig.gene.mat.1 + sig.gene.mat.2
  for ( j in 1:genes){
    if(sum(sig.gene.mat.o[j,] > 0) > 1)
      while(stats::sd(clus.gene.means[j,which(sig.gene.mat.o[j,] > 0)]) < 0.5){
        clus.gene.means[j,which(sig.gene.mat.o[j,] > 0)] <-
          stats::runif(sum(sig.gene.mat.o[j,] > 0), u.l, u.h)*
          sign(stats::rnorm(sum(sig.gene.mat.o[j,]>0)))
      }
  }
  dat.1 <- dat.2 <- matrix(rep(NA,genes*cells), ncol=cells)
  for ( j in 1:cells){
    dat.1[,j] <- stats::rnorm(genes, mean=(gene.means +
                                      clus.gene.means[,clusters1[j]]*
                                      sig.gene.mat.1[,clusters1[j]]), sd=1)
    dat.2[,j] <- stats::rnorm(genes, mean=(gene.means +
                                      clus.gene.means[,clusters2[j]]*
                                      sig.gene.mat.2[,clusters2[j]]), sd=1)
    }
  return( list(dat.1=dat.1, dat.2=dat.2, clusters1=clusters1, clusters2=clusters2,
               sig.gene.mat.1=sig.gene.mat.1, sig.gene.mat.2=sig.gene.mat.2,
               gene.means=gene.means, clus.gene.means=clus.gene.means))
}

