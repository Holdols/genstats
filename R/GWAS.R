#' Computes causal SNP's
#'
#' This function is heavily build onto bigstatsr::big_univLinReg. For more documentation check \cr https://privefl.github.io/bigstatsr/reference/big_univLinReg.html
#' @param G A file backed matrix with SNP data.
#' @param y A vector containing targets for each position.
#' @param p The significance level.
#' @param logreg If TRUE do logistic regression, else Linear regression.
#' @param ncores Amount of cores to use.
#' @return A matrix containing: \cr
#'  - The slopes of each regression \cr
#'  - The standard errors of each slope \cr
#'  - If logreg = T: the number of iterations for each slope. If NA, it means that the algorithm didn't converge, and glm was used instead.
#'  - The t-scores associated with each slope. \cr
#'  - The p-values for each slope \cr
#'  - And a binary vector saying if the position is causal or not given p
#' @export
GWAS <- function(G, y, p, logreg = FALSE, ncores = 1){
  if (logreg == TRUE){
    lm = bigstatsr::big_univLogReg(X = G, y01.train = y, ncores = ncores)
  }
  else{
    lm = bigstatsr::big_univLinReg(X = G, y.train = y, ncores = ncores)
  }
  p_vals = predict(lm, log10 = FALSE)
  causal_estimate = ifelse(p_vals <=p, 1, 0)
  return(cbind(lm, p_vals, causal_estimate))
}

