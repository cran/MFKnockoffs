#' Sample approximate second-order multivariate Gaussian knockoff variables
#' 
#' Samples approximate second-order multivariate Gaussian knockoff variables
#' for the original variables.
#' 
#' @param X normalized n-by-p realization of the design matrix
#' @param method either 'equi', 'sdp' or 'asdp' (default:'asdp')
#' This will be computed according to 'method', if not supplied 
#' @param shrink whether to shrink the estimated covariance matrix (default: FALSE)
#' @return n-by-p matrix of knockoff variables
#'  
#' @family methods for creating knockoffs
#' 
#' @details If the argument \code{shrink} is set to TRUE, a James-Stein-type shrinkage estimator for
#' the covariance matrix is used instead of the traditional maximum-likelihood estimate. This option
#' requires the package \code{corpcor}. Type \code{?corpcor::cov.shrink} for more details.
#' 
#' Even if the argument \code{shrink} is set to FALSE, in the case that the estimated covariance 
#' matrix is not positive-definite, this function will apply some shrinkage.
#' 
#' To use SDP knockoffs, you must have a Python installation with 
#' CVXPY. For more information, see the vignette on SDP knockoffs:
#' \code{vignette('sdp', package='MFKnockoffs')}
#' 
#' @references 
#'   Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
#'   arXiv:1610.02351 (2016).
#'   \href{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}
#'   
#' @export
MFKnockoffs.create.approximate_gaussian <- function(X, method=c("asdp","equi","sdp"), shrink=F) {
  method = match.arg(method)
  # Estimate the mean vectorand covariance matrix
  mu = colMeans(X)
  
  # Estimate the covariance matrix
  if(!shrink) {
    Sigma = cov(X)
    # Verify that the covariance matrix is positive-definite
    if(!is_posdef(Sigma)) {
      shrink=TRUE
    }
  }
  if(shrink) {
    if (!requireNamespace('corpcor', quietly=T))
      stop('corpcor is not installed', call.=F)
    Sigma = tryCatch({suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X,verbose=F)), nrow=ncol(X)))},
                     warning = function(w){}, error = function(e) {
                       stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
                     }, finally = {})
  }

  # Sample the Gaussian knockoffs
  MFKnockoffs.create.gaussian(X, mu, Sigma, method=method)
}