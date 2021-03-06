#' Penalized linear regression statistics for MFKnockoffs
#' 
#' Fit the lasso path and computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of the 
#' regularization parameter \eqn{\lambda} at which the jth variable 
#' and its knockoff enter the penalized linear regression model, respectively.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). It should be numeric.
#' @param ... additional arguments specific to 'glmnet' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#' 
#' @details This function uses \code{glmnet} to compute the lasso path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{100}.
#' 
#' Unless a lambda sequence is provided by the user, this function generates it on a 
#' log-linear scale before calling 'glmnet' (default 'nlambda': 100).
#' 
#' This function is a wrapper around the more general \link{MFKnockoffs.stat.glmnet_lambda_difference}.
#' 
#' For a complete list of the available additional arguments, see \link[glmnet]{glmnet} 
#' or \link[lars]{lars}.
#' 
#' @family statistics for knockoffs
#' 
#' @examples
#' p=100; n=200; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
#' # Basic usage with default arguments
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, 
#'                            statistic=MFKnockoffs.stat.lasso_lambda_difference)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.lasso_lambda_difference
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.lasso_lambda_difference
#' @export
MFKnockoffs.stat.lasso_lambda_difference <- function(X, X_k, y, ...) {
  if( is.numeric(y) ){
    y = as.vector(y)
  } else {
    stop('Knockoff statistic MFKnockoffs.stat.lasso_lambda_difference requires the input y to be a numeric vector')
  }
  
  MFKnockoffs.stat.glmnet_lambda_difference(X, X_k, y, family='gaussian', ...)
}

#' Penalized linear regression statistics for MFKnockoffs
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the penalized linear regression model.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). It should be numeric.
#' @param ... additional arguments specific to 'glmnet' or 'lars' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#'   
#' @details This function uses \code{glmnet} to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The additional \code{nlambda} 
#' parameter can be used to control the granularity of the grid of \eqn{\lambda} values. 
#' The default value of \code{nlambda} is \code{100}.
#' 
#' Unless a lambda sequence is provided by the user, this function generates it on a 
#' log-linear scale before calling 'glmnet' (default 'nlambda': 100).
#' 
#' This function is a wrapper around the more general 
#' \link{MFKnockoffs.stat.glmnet_lambda_difference}.
#' 
#' For a complete list of the available additional arguments, see \link[glmnet]{glmnet}.
#' 
#' @examples
#' p=100; n=200; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
#' # Basic usage with default arguments
#' result = MFKnockoffs.filter(X, y, knockoff=knockoffs,
#'                            statistic=MFKnockoffs.stat.lasso_lambda_signed_max)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.lasso_lambda_signed_max
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.lasso_lambda_signed_max
#' @export
MFKnockoffs.stat.lasso_lambda_signed_max <- function(X, X_k, y, ...) {
  if( is.numeric(y) ){
    y = as.vector(y)
  } else {
    stop('Knockoff statistic MFKnockoffs.stat.lasso_lambda_signed_max requires the input y to be a numeric vector')
  }

  MFKnockoffs.stat.glmnet_lambda_signed_max(X, X_k, y, family='gaussian', ...)
}