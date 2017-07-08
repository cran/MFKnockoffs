#' Penalized logistic regression statistics for MFKnockoffs
#' 
#' Fit the lasso path and computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of the 
#' regularization parameter \eqn{\lambda} at which the jth variable 
#' and its knockoff enter the penalized logistic regression model, respectively.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). It should be either a factor with two levels, 
#' or a two-column matrix of counts or proportions 
#' (the second column is treated as the target class; for a factor, the last level 
#' in alphabetical order is the target class). If y is presented as a vector, 
#' it will be coerced into a factor.
#' @param ... additional arguments specific to 'glmnet' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#' 
#' @details This function uses \code{glmnet} to compute the lasso path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{100}.
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
#' pr = 1/(1+exp(-X %*% beta))
#' y = rbinom(n,1,pr)
#' 
#' knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
#' # Basic usage with default arguments
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, 
#'                            statistic=MFKnockoffs.stat.lasso_lambda_difference_bin)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.lasso_lambda_difference_bin
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.lasso_lambda_difference_bin
#' @export
MFKnockoffs.stat.lasso_lambda_difference_bin <- function(X, X_k, y, ...) {
  MFKnockoffs.stat.glmnet_lambda_difference(X, X_k, y, family='binomial', ...)
}

#' Penalized logistic regression statistics for MFKnockoffs
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the penalized logistic regression model.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). It should be either a factor with two levels, 
#' or a two-column matrix of counts or proportions 
#' (the second column is treated as the target class; for a factor, the last level 
#' in alphabetical order is the target class). If y is presented as a vector, 
#' it will be coerced into a factor.
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
#' pr = 1/(1+exp(-X %*% beta))
#' y = rbinom(n,1,pr)
#' 
#' knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
#' # Basic usage with default arguments
#' result = MFKnockoffs.filter(X, y, knockoff=knockoffs,
#'                            statistic=MFKnockoffs.stat.lasso_lambda_signed_max_bin)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.lasso_lambda_signed_max_bin
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.lasso_lambda_signed_max_bin
#' @export
MFKnockoffs.stat.lasso_lambda_signed_max_bin <- function(X, X_k, y, ...) {
  MFKnockoffs.stat.glmnet_lambda_signed_max(X, X_k, y, family='binomial', ...)
}