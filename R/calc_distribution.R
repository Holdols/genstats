#' Calculates the constants parameters for the conditional distribution belonging to each
#' liabilities of the multivariate gaussian.
#'
#' Of the form \cr
#' mu_mult = \eqn{\Sigma_{12}\Sigma_{22}^{-1}\mu _2} \cr
#' sigma_bar = \eqn{\Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}} \cr
#' It is used as a helper function for LTFH
#' @param sigma The covariance matrix
#' @return a list containing sigma_bar and mu_mult for each liability
#' @examples calc_distribution(get_cov(0.5))
#' @export
calc_distribution = function(sigma){
  out = list()
  for (i in 1:nrow(sigma)){
    mu_mult_bar = solve(sigma[ -i, -i], sigma[i , -i])
    sigma_bar = sigma[i,i] - mu_mult_bar %*% sigma[-i , i]
    out[[i]] = list('sigma'=sigma_bar, 'mu_mult'=mu_mult_bar)
  }
  return(out)
}

