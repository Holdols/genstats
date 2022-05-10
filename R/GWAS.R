#' Computes casusal SNP's
#'
#' This function is heavily build onto bigstatsr::big_univLinReg. For more documentation check \cr https://privefl.github.io/bigstatsr/reference/big_univLinReg.html
#' @param G A FBM with SNP data
#' @param y A vector containing targets for each position.
#' @param ncores Amount of cores to use
#' @param p The significance level standard is 0.05 with bonferroni correction
#' @return a matrix containing the slopes of each regression, \cr
#'  the standard errors of each slope \cr
#'  the t-scores associated with each slope. \cr
#'  the p-values for each slope \cr
#'  and a binary vector saying if the position is causal or not given p
#' @export
GWAS <- function(G, y, ncores = 1, p = 0.05/1000000){ #jeg får hovedpine af at gætte på værdien af p her. hvis i ikke vil bruge 5e-8, brug istedet 0.05/5e-6. det andet er ikke til at læse-
  lm = big_univLinReg(X = G, y.train = y, ncores = ncores)
  p_vals = predict(lm, log10 = FALSE)
  causal_estimate = ifelse(p_vals <=p, 1, 0) # når det bare er 0 eller 1 som skal returneres, er det fint bare at bruge p_vals <= p
  return(cbind(lm, p_vals, causal_estimate)) # jeg vil anbefale jer bare at bruge en tibble i stedet. de er generelt bedre at arbejde med.
}

