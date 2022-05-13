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





no_trunc_estimate = function(mu, sigma){ return(rnorm(1, mu, sqrt(sigma))) }






trunc_estimate = function(j, phenos, K=0.05, mu, sigma){
  if (j == 1) {return(rnorm(1, mu, sqrt(sigma)))
  }
  else {
    if (phenos[j-1] == 'numeric') {phen = phenos[j-1]} else {phen = phenos[j-1][1]}

    crit_bound = pnorm(qnorm(1-K), mu, sqrt(sigma))
    interval = phen * c(crit_bound, 1) + (1-phen)*c(0, crit_bound)

    U = runif(1,interval[1], interval[2])
    return(qnorm(U, mu, sqrt(sigma)))
  }
}








#' Create LTFH estimations for a configuration.
#'
#' @param covmat The covariance matrix
#' @param phenos A binary vector containing the phenotype for each family member
#' of the form \cr c(p_subject, p_parent1, p_parent2, p_sibling1, ... ,p_siblingN) \cr
#' where p_familymember is a binary value (1 or 2) \cr
#' @param K The prevalance of trait
#' @param s_val The starting value of estimates
#' @param start_run Number of iterations before convergence is expected
#' @param all_est If TRUE return the value for each iteration after burn in
#' @return A vector containing LTFH estimate of liabilities of the form \cr
#' c(genetic_liability_subject, liability_subject, liability_parent1, \cr
#' liability_parent2, liability_sibling1, ..., liability_siblingN)
#' @examples
#' LTFH(get_cov(0.5, n_sib = 1), c(1, 1, 0, 0))
#' @export
gibbs_sampl <- function(covmat, phenos, K = 0.05, s_val = 0, start_run=500, all_est=FALSE){
  k = nrow(covmat)
  const_list = calc_distribution(covmat)

  # i formen c(l_g, l, l_p1, l_p2)
  current_liabil = rep(s_val, k)
  liabil_list = list(current_liabil)
  liabil = matrix(s_val, nrow=1, ncol=k)
  total_runs = 1

  while(TRUE) {
    if (total_runs>2 && all(sapply(1:k, function(i) sd(liabil[,i])) / sqrt(total_runs) < 0.01)) break
    for (i in 1:start_run) {
      for (j in 1:k) {

        # Udregner parametre
        sigma = const_list[[j]]$sigma
        mu_mult_ = const_list[[j]]$mu_mult
        mu = mu_mult_ %*% current_liabil[-j] # Udelukker den liability vi er kommet til

        if (K == FALSE) {current_liabil[j] = no_trunc_estimate(mu, sigma)
        } else {current_liabil[j] = trunc_estimate(j, phenos, K, mu, sigma)}


      liabil_list[[length(liabil_list)+1]] = current_liabil
      total_runs = total_runs + 1
    }
    liabil = liabil_list %>% do.call('rbind',.)
  }
  if (all_est){
    return(liabil[-(1:(start_run+1)),])
  }
  return(colMeans(liabil[-(1:(start_run+1)),]))
  }
}


