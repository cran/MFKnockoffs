library(R.matlab)
library(MFKnockoffs)
library(testthat)

setwd("~/Workspace/mf_knockoffs/mf_knockoff_R/MFKnockoffs/tests/testthat")
source("../../R/util.R")

# Load MFKnockoffs for MATLAB
evaluate(matlab, 'addpath ~/Workspace/mf_knockoffs/mf_knockoffs_matlab/ ')
evaluate(matlab, 'addpath ~/Workspace/mf_knockoffs/mf_knockoffs_matlab/glmnet_matlab ')

# Note: This command requires that the 'matlab' binary be on the PATH.
Matlab$startServer(port=9998, options=c("nodesktop", "nodisplay", "nosplash"))

matlab <- Matlab(port=9998)
if (!isOpen(matlab))
  open(matlab)

call_function <- function(matlab, fn, ...) {
  args = list(...)
  names(args) = paste0('R_arg', 1:length(args))
  do.call(setVariable, c(list(matlab), args))
  
  result_name = 'R_result'
  arg_str = paste0(names(args), collapse=', ')
  call_str = paste0(result_name, '=', fn, '(', arg_str, ');')
  evaluate(matlab, call_str)
  getVariable(matlab, result_name)[[1]]
}

test_that('SVDs in R and MATLAB are identical', {
  X = normc(rnorm_matrix(15, 10))
  X.svd = canonical_svd(X)
  
  setVariable(matlab, X=X)
  evaluate(matlab, '[U,S,V] = knockoff.fixed.canonicalSVD(X);')
  X.svd.matlab = getVariable(matlab, c('U','S','V'))
  
  expect_equal(X.svd$u, X.svd.matlab$U)
  expect_equal(X.svd$d, diag(X.svd.matlab$S))
  expect_equal(X.svd$v, X.svd.matlab$V)
})

test_that('equicorrelated knockoffs in R and MATLAB are identical', {
  X = normc(rnorm_matrix(20, 10))
  X_ko = MFKnockoffs.create.fixed(X, method='equi')$X_k
  X_ko.matlab = call_function(matlab, 'knockoff.fixed.create', X, 'equi')
  expect_equal(X_ko, X_ko.matlab)
})

test_that('SDP knockoffs in R and MATLAB are identical', {
  X = normc(rnorm_matrix(20, 10))
  X_ko = MFKnockoffs.create.fixed(X, method='sdp')$X_k
  X_ko.matlab = call_function(matlab, 'knockoff.fixed.create', X, 'sdp')
  # FIXME: MATLAB and R are using different SDP solvers, but is this
  # precision too low?
  expect_equal(X_ko, X_ko.matlab, tol=0.01)
})

test_that('forward selection in R and MATLAB are identical', {
  n = 200; p = 100
  X = normc(rnorm_matrix(n, p))
  y = X %*% rnorm(p) + 0.1 * rnorm(n)
  X_ko = MFKnockoffs.create.fixed(X, method='equi')$X_k
  
  path = MFKnockoffs.stat.forward_selection(X, X_ko, y, omp=FALSE)
  path.matlab = as.vector(call_function(
    matlab, 'knockoff.stats.forwardSelection', X, X_ko, y))
  expect_equal(path, path.matlab)
  
  path = MFKnockoffs.stat.forward_selection(X, X_ko, y, omp=TRUE)
  path.matlab = as.vector(call_function(
    matlab, 'knockoff.stats.forwardSelectionOMP', X, X_ko, y))
  expect_equal(path, path.matlab)
})

test_that('test statistics in R and MATLAB are identical', {
  prob = random_problem(100, 50)
  X = prob$X; y = prob$y
  knock.vars = MFKnockoffs.create.fixed(X, method='equi')
  X = knock.vars$X
  X_ko = knock.vars$X_k
  
  expect_stats_equal <- function(stat.r, stat.matlab) {
    W = stat.r(X, X_ko, y)
    W.matlab = call_function(matlab, stat.matlab, X, X_ko, y)
    expect_equal(W, c(W.matlab), tol=0.01)
  }
  expect_stats_equal(MFKnockoffs.stat.forward_selection, 'knockoff.stats.forwardSelection')
  
  expect_stats_equal <- function(stat.r, stat.matlab) {
    W = stat.r(X, X_ko, y, nlambda=100000)
    W.matlab = call_function(matlab, stat.matlab, X, X_ko, y, 100000)
    expect_equal(W, c(W.matlab), tol=0.01)
  }
  expect_stats_equal(MFKnockoffs.stat.glmnet_lambda_signed_max, 'knockoff.stats.lassoLambdaSignedMax')
})

close(matlab)