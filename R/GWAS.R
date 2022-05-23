#' Computes casusal SNP's
#'
#' This function is heavily build onto bigstatsr::big_univLinReg. For more documentation check \cr https://privefl.github.io/bigstatsr/reference/big_univLinReg.html
#' @param G A file backed matrix with SNP data.
#' @param y A vector containing targets for each position.
#' @param ncores Amount of cores to use.
#' @param p The significance level.
#' @return A matrix containing: \cr
#'  - The slopes of each regression \cr
#'  - The standard errors of each slope \cr
#'  - The t-scores associated with each slope. \cr
#'  - The p-values for each slope \cr
#'  - And a binary vector saying if the position is causal or not given p
#' @export
GWAS <- function(G, y, ncores = 1, p){
  lm = bigstatsr::big_univLinReg(X = G, y.train = y, ncores = ncores)
  p_vals = predict(lm, log10 = FALSE)
  causal_estimate = ifelse(p_vals <=p, 1, 0)
  return(cbind(lm, p_vals, causal_estimate))
}

