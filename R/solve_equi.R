#' Optimization for equi-correlated knockoffs
#' 
#' Solves the optimization problem needed to create equi-correlated knockoffs
#' 
#' @param Sigma A positive-definite covariance matrix
#' @return The solution \eqn{s} to the semidefinite programming problem defined above
#' 
#' @details Computes the closed-form solution to the semidefinite programming problem:
#'  \deqn{ \mathrm{maximize}  \; s \quad
#'        \mathrm{subject} \; \mathrm{to:}   \; 0 <= s <= 1, \;
#'        2\Sigma - sI >= 0 }
#' used to generate equi-correlated knockoffs.
#' 
#' The closed form-solution to this problem is \eqn{s = 2\lambda_{\mathrm{min}}(\Sigma) \land 1}
#' 
#' @family Optimize knockoffs
#' 
#' @export
MFKnockoffs.knocks.solve_equi <- function(Sigma) {
  p = nrow(Sigma)
  tol = 1e-10
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  
  if (requireNamespace('rARPACK', quietly=T)) {
    converged=FALSE
    maxitr=10000
    while (!converged) {
      lambda_min = rARPACK::eigs(G, 1, which = "SR", opts = list(retvec = FALSE, maxitr=maxitr, tol=tol))$values
      if (length(lambda_min)==1) {
        converged = TRUE
      } else {
        if (maxitr>1e8) {
          warning('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                covariance matrix. rARPACK::eigs did not converge. Giving up and computing full SVD with built-in R function.',immediate.=T)
          lambda_min = eigen(G, symmetric=T, only.values = T)$values[p]
          converged=TRUE
        } else {
          warning('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                covariance matrix. rARPACK::eigs did not converge. Trying again with increased number of iterations.',immediate.=T)
          maxitr = maxitr*10
        }
      }
    }
  } else {
    warning('Package rARPACK is not installed. Fast creation of equi-correlated knockoffs is not available.', call.=F,immediate.=T)
    lambda_min = eigen(G, symmetric=T, only.values = T)$values[p]
  }
  
  if (lambda_min<0) {
    stop('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                covariance matrix. The covariance matrix is not positive-definite.')
  }
  
  s = rep(1, nrow(Sigma)) * min(2*lambda_min, 1)
  
  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    diag_s = diag(s*(1-s_eps))
    GInv_s = solve(G,diag_s)
    Sigma_k = round(2*diag_s - diag_s %*% GInv_s,10)
    if (matrixcalc::is.positive.definite(Sigma_k)) {
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}