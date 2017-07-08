#' Cross-validated penalized logistic regression statistics for MFKnockoffs
#' 
#' Fit a logistic regression model via penalized maximum likelihood and cross-validation.
#' Then, compute the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the 
#' jth variable and its knockoff, respectively. The value of the regularization
#' parameter \eqn{\lambda} is selected by cross-validation and computed with glmnet.
#' 
#' @param X original design matrix (size n-by-p).
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). It should be either a factor with two levels, 
#' or a two-column matrix of counts or proportions 
#' (the second column is treated as the target class; for a factor, the last level 
#' in alphabetical order is the target class). If y is presented as a vector, 
#' it will be coerced into a factor.
#' @param cores Number of cores used to compute the knockoff statistics by running cv.glmnet.
#' If not specified, the number of cores is set to approximately half of the number of cores 
#' detected by the parallel package.
#' @param ... additional arguments specific to 'glmnet' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#' 
#' @details This function uses the \code{glmnet} package to fit the penalized logistic regression path.
#' 
#' This function is a wrapper around the more general \link{MFKnockoffs.stat.glmnet_coef_difference}.
#' 
#' The knockoff statistics \eqn{W_j} are constructed by taking the difference 
#' between the coefficient of the j-th variable and its knockoff.
#'  
#' By default, the value of the regularization parameter is chosen by 10-fold cross-validation.
#' 
#' The optional \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{100},
#' where \code{p} is the number of columns of \code{X}.
#' 
#' For a complete list of the available additional arguments, see \link[glmnet]{cv.glmnet}
#' and \link[glmnet]{glmnet}.
#' 
#' @family statistics for knockoffs
#' 
#' @examples
#' p=100; n=200; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' pr = 1/(1+exp(-X %*% beta))
#' y = rbinom(n,1,pr)
#' 
#' knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
#' # Basic usage with default arguments
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, 
#'                            statistic=MFKnockoffs.stat.lasso_coef_difference_bin)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.lasso_coef_difference_bin
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.lasso_coef_difference_bin
#' @export
MFKnockoffs.stat.lasso_coef_difference_bin <- function(X, X_k, y, cores=2, ...) {
  if (!is.factor(y) && !is.numeric(y)) {
    stop('Input y must be either of numeric or factor type')
  }
  MFKnockoffs.stat.glmnet_coef_difference(X, X_k, y, family='binomial', cores=cores, ...)
}
