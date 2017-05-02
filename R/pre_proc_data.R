#' Pre-process Data
#'
#' This function pre-process the data so that SparseDC can be applied.
#' SparseDC requires data that have been normalized for sequencing depth,
#' log-transformed and centralized on a gene-by-gene basis. For the sequencing
#' depth normalization we recommend that users use one of the many methods
#' developed for normalizing scRNA-Seq data prior to using SparseDC and so
#' can set \code{norm = FALSE}. However, here we normalize the data by dividing
#' by the total number of reads. This function log transforms the data by
#' applying \code{log(x + 1)} to each of the data sets. By far the most
#' important pre-processing step for SparseDC is the centralization of the data.
#' Having centralized data is a core component of the SparseDC algorithm and is
#' necessary for both accurate clustering of the cells and identifying marker
#' genes. We therefore recommend that all users centralize their data using
#' this function and that only experienced users set \code{center = FALSE}.
#' @param dat1 The data for the first condition with samples (cells) as columns
#' and features (genes) as rows.
#' @param dat2 The data for the second condition with samples (cells) as
#' columns and features (genes) as rows.
#' @param log This parameter controls whether the data is transformed using
#' \code{log(x + 1)}. The default value is \code{TRUE}.
#' @param center This parameter controls whether the data is centered on a gene
#' by gene basis. We recommend all users center their data prior to applying
#' SparseDC and only experienced users should set this as \code{FALSE}. The
#' default value is \code{TRUE}.
#' @param norm This parameter controls whether the data is normalized for
#' sequencing depth by dividing each column by the total number of reads for
#' that sample. We recommend that user use one of the many methods for
#' normalizing scRNA-Seq data and so set this as \code{FALSE}. The default
#' value is \code{TRUE}
#' @return This function returns the two pre-processed datasets stored as a
#' list
#' @examples
#' set.seed(10)
#' # Select small dataset for example
#' data_test <- data_biase[1:100,]
#' # Split data into condition A and B
#' data_A <- data_test[ , which(condition_biase == "A")]
#' data_B <- data_test[ , which(condition_biase == "B")]
#' # Pre-process the data
#' pre_data <- pre_proc_data(data_A, data_B, norm = FALSE, log = TRUE,
#' center = TRUE)
#' # Extract Data
#' pdata_A <- pre_data[[1]]
#' pdata_B <- pre_data[[2]]
#'
pre_proc_data <- function(dat1, dat2, norm = TRUE, log = TRUE, center = TRUE){
  dat1_temp <- dat1  # Create temporary data 1
  dat2_temp <- dat2  # Create temporary data 2
  if (norm == TRUE){  # If data is to be normalized:
    data_1_tc <- colSums(dat1_temp)  # Calculate column sums for data 1
    data_2_tc <- colSums(dat2_temp)  # Calculate column sums for data 2
    dat1_temp <- dat1_temp %*% diag(1/data_1_tc)  # Normalize data 1
    dat2_temp <- dat2_temp %*% diag(1/data_2_tc)  # Normalize data 2
  }
  if ( log == TRUE){  # If data is to be log-transformed:
    dat1_temp <- log(dat1_temp + 1)  # Log-transform data 1
    dat2_temp <- log(dat2_temp + 1)  # Log-transform data 2
  }
  if ( center == TRUE){  # If data is to be centered:
    ncol_1 <- ncol(dat1_temp)  # Store number of samples in data 1
    c.dat <- t(scale(t(cbind(dat1_temp,dat2_temp)), scale=F)) # Center the data
    dat1_temp <- c.dat[,1:ncol_1]  # Extract data 1
    dat2_temp <- c.dat[,(ncol_1+1):ncol(c.dat)]  # Extract data 2
  }
  return(list(pdat1 = dat1_temp, pdat2=dat2_temp))  # Return preprocessed data
}




































