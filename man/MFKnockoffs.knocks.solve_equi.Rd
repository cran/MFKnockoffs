% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_equi.R
\name{MFKnockoffs.knocks.solve_equi}
\alias{MFKnockoffs.knocks.solve_equi}
\title{Optimization for equi-correlated knockoffs}
\usage{
MFKnockoffs.knocks.solve_equi(Sigma)
}
\arguments{
\item{Sigma}{A positive-definite covariance matrix}
}
\value{
The solution \eqn{s} to the semidefinite programming problem defined above
}
\description{
Solves the optimization problem needed to create equi-correlated knockoffs
}
\details{
Computes the closed-form solution to the semidefinite programming problem:
 \deqn{ \mathrm{maximize}  \; s \quad
       \mathrm{subject} \; \mathrm{to:}   \; 0 <= s <= 1, \;
       2\Sigma - sI >= 0 }
used to generate equi-correlated knockoffs.

The closed form-solution to this problem is \eqn{s = 2\lambda_{\mathrm{min}}(\Sigma) \land 1}
}
\seealso{
Other Optimize knockoffs: \code{\link{MFKnockoffs.knocks.solve_asdp}},
  \code{\link{MFKnockoffs.knocks.solve_sdp}}
}
