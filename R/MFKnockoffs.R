#' @docType package
#' @name MFKnockoffs
#' @import stats methods
NULL

#' Model-Free Knockoff Filter
#' 
#' Run the model-free knockoff procedure from start to finish, selecting variables
#' relevant for predicting the outcome of interest.
#' 
#' This function creates the knockoffs, computes the test statistics, and selects 
#' variables. 
#' It is the main entry point for the model-free knockoff package.
#' 
#' @param X matrix or data frame of predictors
#' @param y response vector
#' @param knockoffs the method used to construct knockoffs for the X variables.
#' It must be a function taking a n-by-p matrix X as input and returning a n-by-p matrix of knockoff variables
#' @param statistic the test statistic (by default, a lasso statistic with cross validation). See the
#'  Details section for more information.
#' @param q target FDR (false discovery rate)
#' @param threshold either 'knockoff+' or 'knockoff' (default: 'knockoff+').
#'
#' @return An object of class "MFKnockoffs.result". This object is a list 
#'  containing at least the following components:
#'  \item{X}{matrix of original variables}
#'  \item{X_k}{matrix of knockoff variables}
#'  \item{statistic}{computed test statistics}
#'  \item{threshold}{computed selection threshold}
#'  \item{selected}{named vector of selected variables}
#'
#' @details
#' 
#' The parameter \code{knockoffs} controls how knockoff variables are created.
#' By default, a multivariate normal distribution is fitted to the original
#' variables in X. The estimated mean vector and covariance matrix are used
#' to generate second-order approximate Gaussian knockoffs.
#' In general, \code{knockoffs} should be a function taking a n-by-p matrix of
#' observed variables X and returning a n-by-p matrix of knockoff variables.
#' Two optional functions for creating knockoffs are provided with this package.
#' 
#' If the rows of X are distributed as a multivariate Gaussian with known parameters,
#' then the function \code{MFKnockoffs.create.gaussian} can be used to generate
#' valid Gaussian knockoff variables, as shown in the examples below.
#' 
#' If the design matrix X is assumed to be fixed instead of random, one can create
#' knockoff variables using the function \code{MFKnockoffs.create.fixed}. This 
#' corresponds to the original framework of the (non Model-Free) knockoff filter.
#' 
#' For more information about creating knockoffs, type \code{??MFKnockoffs.create}.
#' 
#' The default test statistic is \link{MFKnockoffs.stat.glmnet_coef_difference}.
#' For a complete list of the statistics provided with this package, 
#' type \code{??MFKnockoffs.stat}.
#' 
#' It is also possible to provide custom test statistics.
#' An example can be found in the vignette.
#' 
#' @references 
#'   Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
#'   arXiv:1610.02351 (2016).
#'   \href{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}
#'   
#'   Barber and Candes,
#'   Controlling the false discovery rate via knockoffs. 
#'   Ann. Statist. 43 (2015), no. 5, 2055--2085.
#'   \href{https://projecteuclid.org/euclid.aos/1438606853}{https://projecteuclid.org/euclid.aos/1438606853}
#' 
#' @examples
#' p=100; n=200; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' # Basic usage with default arguments
#' result = MFKnockoffs.filter(X, y)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
#' k_stat = function(X, X_k, y) MFKnockoffs.stat.glmnet_coef_difference(X, X_k, y, nfolds=5)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' 
#' @export
MFKnockoffs.filter <- function(X, y,
                              knockoffs=MFKnockoffs.create.approximate_gaussian,
                              statistic=MFKnockoffs.stat.glmnet_coef_difference, 
                              q=0.10,
                              threshold=c('knockoff+','knockoff')
                              ) {
  
  # Validate input types.
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)  
  } else if (is.matrix(X))
    X.names = colnames(X)
  else
    stop('Input X must be a matrix or data frame')
  if (!is.factor(y) && !is.numeric(y)) {
    stop('Input y must be either of numeric or factor type')
  }
  if( is.numeric(y) ) y = as.vector(y)
  
  # Validate input dimensions
  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n)

  # If fixed-design knockoffs are being used, provive them with the response vector
  # in order to augment the data with new rows if necessary
  if( identical(knockoffs, MFKnockoffs.create.fixed) )
    knockoffs = function(x) MFKnockoffs.create.fixed(x, y=y)
  
  # Create knockoff variables
  knock_variables = knockoffs(X)
  
  # If fixed-design knockoffs are being used, update X and Y with the augmented observations (if present)
  if (is(knock_variables,"MFKnockoffs.variables")){
    X   = knock_variables$X
    X_k = knock_variables$X_k
    y   = knock_variables$y
    rm(knock_variables)
  }
  else if (is(knock_variables,"matrix")){
    X_k = knock_variables
    rm(knock_variables)
  }
  else {
    stop('Knockoff variables of incorrect type')
  }
  
  # Compute knockoff statistics
  W = statistic(X, X_k, y)
  
  # Run the knockoff filter
  t = MFKnockoffs.threshold(W, q, threshold)
  selected = which(W >= t)
  if (!is.null(X.names))
    names(selected) = X.names[selected]
  
  # Package up the results.
  structure(list(call = match.call(),
                 X = X,
                 X_k = X_k,
                 y = y,
                 statistic = W,
                 threshold = t,
                 selected = selected),
            class = 'MFKnockoffs.result')
}

#' Threshold for the knockoff filter
#' 
#' Computes the threshold for the knockoff filter.
#' 
#' @param W the test statistics
#' @param q target FDR (false discovery rate)
#' @param method either 'knockoff or 'knockoff+'
#' @return The threshold for variable selection.
#' 
#' @export
MFKnockoffs.threshold <- function(W, q=0.10, method=c('knockoff+','knockoff')) {
  offset = switch(match.arg(method),
                  'knockoff' = 0, 'knockoff+' = 1)
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= q)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

#' Print results for the knockoff filter
#' 
#' Prints the list of variables selected by the knockoff filter and the corresponding function call.
#' 
#' @param x the output of a call to MFKnockoffs.filter
#' @param ... unused
#' 
#' @method print MFKnockoffs.result
#' @export
print.MFKnockoffs.result <- function(x, ...) {
  cat('Call:\n')
  print(x$call)
  cat('\nSelected variables:\n')
  print(x$selected)
}

#' Verify dependencies for chosen statistics
#' 
#' @param statistic the knockoff statistic chosen by the user
#' 
#' @keywords internal
verify_stat_depends <- function(statistic) {
  
}