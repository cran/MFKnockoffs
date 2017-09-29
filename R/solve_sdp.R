#' Optimization for SDP knockoffs
#'
#' Solves the optimization problem needed to create SDP knockoffs using an interior point method
#' 
#' @param Sigma A positive-definite correlation matrix
#' @param maxit The maximum number of iterations for the solver (default: 1000)
#' @param gaptol Tolerance for duality gap as a fraction of the value of the objective functions (default 1e-6)
#' @return The solution \eqn{s} to the semidefinite programming problem defined above
#'
#' @details
#' Solves the semidefinite programming problem:
#'
#'   \deqn{ \mathrm{maximize}      \; \mathrm{sum}(s) \quad
#'           \mathrm{subject} \; \mathrm{to}    0 <= s <= 1, \;
#'                                  2\Sigma - \mathrm{diag}(s) >= 0}
#'
#' If the matrix Sigma supplied by the user is a non-scaled covariance matrix 
#' (i.e. its diagonal entries are not all equal to 1), then the appropriate scaling is applied before
#' solving the SDP defined above. The result is then scaled back before being returned, as to match 
#' the original scaling of the covariance matrix supplied by the user.
#' 
#' @family Optimize knockoffs
#' 
#' @export
MFKnockoffs.knocks.solve_sdp <- function(Sigma, gaptol=1e-6, maxit=1000) {
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]

  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    stop('The covariance matrix is not positive-definite: cannot solve SDP',immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),] 
  Cs = c(2*G)
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution')
  }
  
  # # Call solver with increasing number of iterations until it converges
  # This may no longer be necessary as we switched to a better solver
  # converged = F
  # converged_tried = 1
  # while (! converged ) {
  #   max_iters_run = max_iters * converged_tried
  #   control <- list(eps = eps, max_iters = max_iters_run, verbose=F, normalize=F)
  #   sol <- scs::scs(A, b, obj, cone, control)
  #   if ( sol$info$statusVal == 1 ) {
  #     converged = T
  #   }
  #   else {
  #     converged_tried = converged_tried+1
  #     warning(paste('Conic solver (scs) for SDP knockoffs did not converge after', max_iters_run,
  #                   'iterations. Trying again with double number of iterations'),immediate.=T)
  #   }
  # }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  psd = 0
  s_eps = 1e-8
  while (psd==0) {
    if (is_posdef(2*G-diag(s*(1-s_eps),length(s)))) {
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  
  # Verify that the solution is correct
  if (max(s)==0) {
   warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}

#' Vectorize a matrix into the SCS format
#'  
#' @rdname vectorize_matrix
#' @keywords internal
MFKnockoffs.knocks.vectorize_matrix = function(M) {
  # Scale the off-diagonal entries by sqrt(2)
  vectorized_matrix = M
  vectorized_matrix[lower.tri(M,diag=FALSE)] = M[lower.tri(M,diag=FALSE)] * sqrt(2)
  # Stack the lower triangular elements column-wise
  vectorized_matrix = vectorized_matrix[lower.tri(vectorized_matrix,diag=TRUE)]
}