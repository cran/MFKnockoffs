% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_glmnet_cv.R
\name{MFKnockoffs.stat.glmnet_coef_difference}
\alias{MFKnockoffs.stat.glmnet_coef_difference}
\title{Cross-validated GLM statistics for MFKnockoffs}
\usage{
MFKnockoffs.stat.glmnet_coef_difference(X, X_k, y, family = "gaussian",
  cores = 2, ...)
}
\arguments{
\item{X}{original design matrix (size n-by-p)}

\item{X_k}{knockoff matrix (size n-by-p)}

\item{y}{response vector (length n). Quantitative for family="gaussian", 
or family="poisson" (non-negative counts). For family="binomial" 
should be either a factor with two levels, or a two-column matrix of counts 
or proportions (the second column is treated as the target class; for a factor, 
the last level in alphabetical order is the target class). For family="multinomial", 
can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. 
For either "binomial" or "multinomial", if y is presented as a vector, it will 
be coerced into a factor. For family="cox", y should be a two-column matrix with 
columns named 'time' and 'status'. The latter is a binary variable, with '1' 
indicating death, and '0' indicating right censored. The function Surv() in 
package survival produces such a matrix. For family="mgaussian", y is a matrix 
of quantitative responses.}

\item{family}{Response type (see above)}

\item{cores}{Number of cores used to compute the knockoff statistics by running cv.glmnet.
Unless otherwise specified, the number of cores is set equal to two (if available).}

\item{...}{additional arguments specific to 'cv.glmnet' (see Details)}
}
\value{
A vector of statistics \eqn{W} (length p)
}
\description{
Fit a generalized linear model via penalized maximum likelihood and cross-validation.
Then, compute the difference statistic
  \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the 
jth variable and its knockoff, respectively. The value of the regularization
parameter \eqn{\lambda} is selected by cross-validation and computed with glmnet.
}
\details{
This function uses the \code{glmnet} package to fit a generalized linear model
via penalized maximum likelihood.

The knockoff statistics \eqn{W_j} are constructed by taking the difference 
between the coefficient of the j-th variable and its knockoff.
 
By default, the value of the regularization parameter is chosen by 10-fold cross-validation.

The default response family is 'gaussian', for a linear regression model.
Different response families (e.g. 'binomial') can be specified by passing an
optional parameter 'family'.

The optional \code{nlambda} parameter can be used to control the granularity of the 
grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{100},
where \code{p} is the number of columns of \code{X}.

If the family is 'binomial' and a lambda sequence is not provided by the user, 
this function generates it on a log-linear scale before calling 'glmnet'.

For a complete list of the available additional arguments, see \link[glmnet]{cv.glmnet}
and \link[glmnet]{glmnet}.
}
\examples{
p=100; n=200; k=15
mu = rep(0,p); Sigma = diag(p)
X = matrix(rnorm(n*p),n)
nonzero = sample(p, k)
beta = 3.5 * (1:p \%in\% nonzero)
y = X \%*\% beta + rnorm(n)

knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
# Basic usage with default arguments
result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, 
                           statistic=MFKnockoffs.stat.glmnet_coef_difference)
print(result$selected)

# Advanced usage with custom arguments
foo = MFKnockoffs.stat.glmnet_coef_difference
k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
print(result$selected)

}
\seealso{
Other statistics for knockoffs: \code{\link{MFKnockoffs.stat.forward_selection}},
  \code{\link{MFKnockoffs.stat.glmnet_lambda_difference}},
  \code{\link{MFKnockoffs.stat.lasso_coef_difference_bin}},
  \code{\link{MFKnockoffs.stat.lasso_coef_difference}},
  \code{\link{MFKnockoffs.stat.lasso_lambda_difference_bin}},
  \code{\link{MFKnockoffs.stat.lasso_lambda_difference}},
  \code{\link{MFKnockoffs.stat.random_forest}},
  \code{\link{MFKnockoffs.stat.sqrt_lasso}},
  \code{\link{MFKnockoffs.stat.stability_selection}}
}
