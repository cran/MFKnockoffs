#' Create fixed-design knockoff variables
#' 
#' Creates fixed-design knockoff variables for the original variables.
#' 
#' @param X normalized n-by-p design matrix (n >= 2p)
#' @param method either 'equi' or 'sdp'  (default:'sdp')
#' @param sigma noise level, used to augment the data with extra rows if necessary (default: NULL)
#' @param y vector of observed responses, used to estimate the noise level if 'sigma' is not provided (default: NULL)
#' @param randomize whether the knockoffs are deterministic or randomized (default:False)
#' @return An object of class "MFKnockoffs.variables". This object is a list 
#'  containing at least the following components:
#'  \item{X}{The n-by-p matrix of original variables (possibly augmented or transformed)}
#'  \item{X_k}{The n-by-p matrix of knockoff variables}
#'  \item{y}{The vector of observed responses (possibly augmented) }
#' 
#' @family methods for creating knockoffs
#' 
#' @references 
#'   Barber and Candes,
#'   Controlling the false discovery rate via knockoffs. 
#'   Ann. Statist. 43 (2015), no. 5, 2055--2085.
#'   \href{https://projecteuclid.org/euclid.aos/1438606853}{https://projecteuclid.org/euclid.aos/1438606853}
#' 
#' 
#' Fixed-design knockoff assume a linear regression model for Y|X. Moreover, they only guarantee
#' FDR control with statistics satisfying the "sufficiency" property. In particular, the default
#' statistics with cross-validated lasso does not satisfy this property and should not be used
#' with fixed-design knockoffs.
#' 
#' @export
MFKnockoffs.create.fixed <- function(X, method=c('sdp','equi'), sigma=NULL, y=NULL, randomize=F) {
  method = match.arg(method)
  
  # Validate dimensions, if using fixed-design knockoffs
  n = nrow(X); p = ncol(X)
  if (n <= p)
    stop('Input X must have dimensions n > p')
  else if (n < 2*p) {
    warning('Input X has dimensions p < n < 2p. ',
            'Augmenting the model with extra rows.',immediate.=T)
    X.svd = svd(X, nu=n, nv=0)
    u2 = X.svd$u[,(p+1):n]
    X = rbind(X, matrix(0, 2*p-n, p))
    if( is.null(sigma) ) {
      if( is.null(y) ) {
        stop('Either the noise level "sigma" or the response variables "y" must
             be provided in order to augment the data with extra rows')
      }
      else{
        sigma = sqrt(mean((t(u2) %*% y)^2)) # = sqrt(RSS/(n-p))
      }
    }
    if (randomize)
      y.extra = rnorm(2*p-n, sd=sigma)
    else
      y.extra = with_seed(0, rnorm(2*p-n, sd=sigma))
    y = c(y, y.extra)
  }
  # Normalize X, if using fixed-design knockoffs
  X = normc(X)
  
  X_k = switch(match.arg(method), 
               "equi" = create_equicorrelated(X,randomize),
               "sdp"  = create_sdp(X,randomize)
              )
  structure(list(X=X, X_k=X_k, y=y), class='MFKnockoffs.variables')
}

# Create equicorrelated knockoffs.
create_equicorrelated <- function(X, randomize) {
  # Compute SVD and U_perp.
  X.svd = decompose(X, randomize)
  
  # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
  # have the same value 1-s.
  if (any(X.svd$d <= 1e-5 * max(X.svd$d)))
    stop(paste('Data matrix is rank deficient.',
               'Equicorrelated knockoffs will have no power.'))
  lambda_min = min(X.svd$d)^2
  s = min(2*lambda_min, 1)
  
  # Construct the knockoff according to Equation 1.4.
  s_diff = pmax(0, 2*s - (s/X.svd$d)^2) # can be negative due to numerical error
  X_ko = (X.svd$u %*diag% (X.svd$d - s / X.svd$d) +
          X.svd$u_perp %*diag% sqrt(s_diff)) %*% t(X.svd$v)
}

# Create SDP knockoffs.
create_sdp <- function(X, randomize) {
  # Compute SVD and U_perp.
  X.svd = decompose(X, randomize)
  
  # Check for rank deficiency.
  tol = 1e-5
  d = X.svd$d
  d_inv = 1 / d
  d_zeros = d <= tol*max(d)
  if (any(d_zeros)) {
    warning(paste('Data matrix is rank deficient.',
                  'Model is not identifiable, but proceeding with SDP knockoffs'),immediate.=T)
    d_inv[d_zeros] = 0
  }
  
  # Compute the Gram matrix and its (pseudo)inverse.
  G = (X.svd$v %*diag% d^2) %*% t(X.svd$v)
  G_inv = (X.svd$v %*diag% d_inv^2) %*% t(X.svd$v)
  
  # Optimize the parameter s of Equation 1.3 using SDP.
  s = MFKnockoffs.knocks.solve_sdp(G)
  s[s <= tol] = 0
  
  # Construct the knockoff according to Equation 1.4:
  C.svd = canonical_svd(2*diag(s) - (s %diag*% G_inv %*diag% s))
  X_ko = X - (X %*% G_inv %*diag% s) + 
    (X.svd$u_perp %*diag% sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}

# Compute the SVD of X and construct an orthogonal matrix U_perp such that
# U_perp * U = 0.
decompose <- function(X, randomize) {
  n = nrow(X); p = ncol(X)
  stopifnot(n >= 2*p)
  
  result = canonical_svd(X)
  Q = qr.Q(qr(cbind(result$u, matrix(0,n,p))))
  u_perp = Q[,(p+1):(2*p)]
  if (randomize) {
      Q = qr.Q(qr(rnorm_matrix(p,p)))
      u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}