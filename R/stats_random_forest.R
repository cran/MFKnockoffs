#' Random forest statistics for MFKnockoffs
#' 
#' Computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the random forest feature importances
#' of the jth variable and its knockoff, respectively.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). If a factor, classification is assumed, 
#' otherwise regression is assumed.
#' @param ... additional arguments specific to 'ranger' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#' 
#' @details This function uses the \code{ranger} package to compute variable 
#' importance measures. The importance of a variable is measured as the total decrease
#' in node impurities from splitting on that variable, averaged over all trees. 
#' For regression, the node impurity is measured by residual sum of squares.
#' For classification, it is measured by the Gini index.
#' 
#' For a complete list of the available additional arguments, see \link[ranger]{ranger}. 
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
#'                            statistic=MFKnockoffs.stat.random_forest)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.random_forest
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nodesize=5)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.random_forest
#' @export
MFKnockoffs.stat.random_forest <- function(X, X_k, y, ...) {
  if (!requireNamespace('ranger', quietly=T))
    stop('ranger is not installed', call.=F)
  
  Z = random_forest_importance(cbind(X, X_k), y) 
  p = ncol(X)
  orig = 1:p
  abs(Z[orig]) - abs(Z[orig+p])
}

#' @keywords internal
random_forest_importance <- function(X, y, ...) {
  df = data.frame(y=y, X=X)
  rfFit = ranger::ranger(y~., data=df, importance="impurity", write.forest=F, ...)
  as.vector(rfFit$variable.importance)
}