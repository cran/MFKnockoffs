#' Model-free controlled variable selection
#'
#' Model-free knockoffs provide a general and powerful tool to perform high-dimensional 
#' controlled variable selection.
#' 
#' The method constructs artificial 'knockoff copies' of the variables 
#' in a statistical model and then selects those variables that are clearly better 
#' than their corresponding fake copies.
#' A wide range of statistics and machine learning tools can be exploited to estimate the 
#' importance of each feature, while guaranteeing finite-sample control of the false
#' discovery rate (FDR). 
#' This model-free approach makes it possible to use knockoffs for data originating 
#' from any conditional model (Y|X), no matter how high-dimensional, provided that the 
#' distribution of the covariates (X) is known.
#' 
#' For more information, see the website below and the accompanying paper.
#' 
#' \url{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}
#' 
#' @name MFKnockoffs
#' @docType package
NULL
