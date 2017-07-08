#' GLM statistics for MFKnockoffs
#' 
#' Fit a generalized linear model via penalized maximum likelihood and
#' computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of the 
#' regularization parameter \eqn{\lambda} at which the jth variable 
#' and its knockoff enter the model, respectively.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). Quantitative for family="gaussian", 
#' or family="poisson" (non-negative counts). For family="binomial" 
#' should be either a factor with two levels, or a two-column matrix of counts 
#' or proportions (the second column is treated as the target class; for a factor, 
#' the last level in alphabetical order is the target class). For family="multinomial", 
#' can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. 
#' For either "binomial" or "multinomial", if y is presented as a vector, it will 
#' be coerced into a factor. For family="cox", y should be a two-column matrix with 
#' columns named 'time' and 'status'. The latter is a binary variable, with '1' 
#' indicating death, and '0' indicating right censored. The function Surv() in 
#' package survival produces such a matrix. For family="mgaussian", y is a matrix 
#' of quantitative responses.
#' @param family Response type (see above)
#' @param ... additional arguments specific to 'glmnet' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#' 
#' @details This function uses \code{glmnet} to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{100}.
#' 
#' If the family is 'binomial' and a lambda sequence is not provided by the user, 
#' this function generates it on a log-linear scale before calling 'glmnet'.
#' 
#' The default response family is 'gaussian', for a linear regression model.
#' Different response families (e.g. 'binomial') can be specified by passing an
#' optional parameter 'family'.
#' 
#' For a complete list of the available additional arguments, see \link[glmnet]{glmnet}.
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
#'                            statistic=MFKnockoffs.stat.glmnet_lambda_difference)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.glmnet_lambda_difference
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.glmnet_lambda_difference
#' @export
MFKnockoffs.stat.glmnet_lambda_difference <- function(X, X_k, y, family='gaussian', ...) {
  Z = lasso_max_lambda(cbind(X, X_k), y, method='glmnet', family=family, ...)
  p = ncol(X)
  orig = 1:p
  Z[orig] - Z[orig+p]
}

#' GLM statistics for MFKnockoffs
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the generalized linear model.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). Quantitative for family="gaussian", 
#' or family="poisson" (non-negative counts). For family="binomial" 
#' should be either a factor with two levels, or a two-column matrix of counts 
#' or proportions (the second column is treated as the target class; for a factor, 
#' the last level in alphabetical order is the target class). For family="multinomial", 
#' can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. 
#' For either "binomial" or "multinomial", if y is presented as a vector, it will 
#' be coerced into a factor. For family="cox", y should be a two-column matrix with 
#' columns named 'time' and 'status'. The latter is a binary variable, with '1' 
#' indicating death, and '0' indicating right censored. The function Surv() in 
#' package survival produces such a matrix. For family="mgaussian", y is a matrix 
#' of quantitative responses.
#' @param family Response type (see above)
#' @param ... additional arguments specific to 'glmnet' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#'   
#' @details This function uses \code{glmnet} to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The additional \code{nlambda} 
#' parameter can be used to control the granularity of the grid of \eqn{\lambda} values. 
#' The default value of \code{nlambda} is \code{100}.
#' 
#' If the family is 'binomial' and a lambda sequence is not provided by the user, 
#' this function generates it on a log-linear scale before calling 'glmnet'.
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
#'                            statistic=MFKnockoffs.stat.glmnet_lambda_signed_max)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.glmnet_lambda_signed_max
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.glmnet_lambda_signed_max
#' @export
MFKnockoffs.stat.glmnet_lambda_signed_max <- function(X, X_k, y, family='gaussian', ...) {
  Z = lasso_max_lambda(cbind(X, X_k), y, method='glmnet', family=family, ...)
  p = ncol(X)
  orig = 1:p
  pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
}

#' @rdname lasso_max_lambda
#' @keywords internal
lasso_max_lambda_lars <- function(X, y, ...) {
  if (!requireNamespace('lars', quietly=T))
    stop('lars is not installed', call.=F)
  
  fit <- lars::lars(X, y, normalize=T, intercept=F, ...)
  lambda <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    entry <- fit$entry[j]
    if (entry > 0) lambda[j] <- fit$lambda[entry]
  }
  return(lambda)
}

#' @rdname lasso_max_lambda
#' @keywords internal
lasso_max_lambda_glmnet <- function(X, y, nlambda=200, intercept=T, standardize=T, ...) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  
  # Standardize the variables
  if( standardize ){
    X = normc(X)
  }
    
  n = nrow(X); p = ncol(X)
  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family
  
  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }

  fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=intercept, 
                        standardize=F, standardize.response=F, ...)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices] * n)
}

#' Maximum lambda in lasso model
#' 
#' Computes the earliest (largest) lambda's for which predictors enter the
#' lasso model.
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @param method either 'glmnet' or 'lars'
#' @return vector of maximum lambda's
#' 
#' @keywords internal
lasso_max_lambda <- function(X, y, method=c('glmnet','lars'), ...) {
  switch(match.arg(method), 
         glmnet = lasso_max_lambda_glmnet(X,y,...),
         lars = lasso_max_lambda_lars(X,y,...)
         )
}