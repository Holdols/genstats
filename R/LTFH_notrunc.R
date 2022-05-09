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
#' @noRd
calc_distribution = function(sigma){
  out = list()
  for (i in 1:nrow(sigma)){
    mu_mult_bar = solve(sigma[ -i, -i], sigma[i , -i])
    sigma_bar = sigma[i,i] - mu_mult_bar %*% sigma[-i , i]
    out[[i]] = list('sigma'=sigma_bar, 'mu_mult'=mu_mult_bar)
  }
  return(out)
}


#' Create LTFH estimations without truncation.
#'
#' @param covmat The covariance matrix
#' @param N Amount of rows in the matrix that the estimations are saved in
#' @param K The prevalance of trait
#' @param s_val The starting value of estimates
#' @return A vector containing LTFH estimate of liabilities of the form \cr
#' c(genetic_liability_subject, liability_subject, liability_parent1, \cr
#' liability_parent2, liability_sibling1, ..., liability_siblingN)
#' @examples
#' LTFH(get_cov(0.5, n_sib = 1))
#' @export
LTFH_notrunc <- function(covmat, N=3000, K = 0.05, s_val = 0){

  k = nrow(covmat)

  const_list = calc_distribution(covmat)

  liabil = matrix(s_val, nrow = N, ncol = k) # i formen c(l_g, l, l_p1, l_p2)
  mu_vec = matrix(s_val, nrow = N-1, ncol = k)
  current_liabil = rep(s_val, k)

  for (i in 1:(N)) {
    for (j in 1:k) {
      # Udregner parametre
      sigma = const_list[[j]]$sigma
      mu_mult = const_list[[j]]$mu_mult
      mu = mu_mult %*% current_liabil[-j] # Udelukker den liability vi er kommet til
      mu_vec[i,j] = mu
      current_liabil[j] = rnorm(1, mu, sqrt(sigma))
    }
    liabil[i,] = current_liabil
  }


  return(list('liabilities'=liabil, 'mu'=mu_vec))
}
