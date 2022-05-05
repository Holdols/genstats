#' Create LTFH estimations
#'
#' @param covmat The covariance matrix
#' @param phenos A binary vector containing the phenotype for each family member
#' of the form \cr c(p_subject, p_parent1, p_parent2, p_sibling1, ... ,p_siblingN) \cr
#' where p_familymember is a binary value (1 or 2) \cr
#' @param N Amount of rows in the matrix that the estimations are saved in
#' @param K The prevalance of trait
#' @param s_val The starting value of estimates
#' @param min_run The minimum amount of runs to do after burn in
#' @param all_est If TRUE return the value for each iteration after burn in
#' @return A vector containing LTFH estimate of liabilities of the form \cr
#' c(genetic_liability_subject, liability_subject, liability_parent1, \cr
#' liability_parent2, liability_sibling1, ..., liability_siblingN)
#' @examples
#' LTFH(get_cov(0.5, n_sib = 1), c(1, 1, 0, 0))
#' @export
LTFH <- function(covmat, phenos, N=10000, K = 0.05, s_val = 0, min_run=400, all_est=FALSE){


  k = nrow(covmat)

  const_list = calc_distribution(covmat)

  liabil = matrix(s_val, nrow = N, ncol = k) # i formen c(l_g, l, l_p1, l_p2)

  current_liabil = rep(s_val, k)

  start_run = 500
  s = 0
  p = 0

  while(TRUE) {
    p = p +1
    for (i in 2:start_run) {
      s = s + 1
      for (j in 1:k) {
        # Udregner parametre
        sigma = const_list[[j]]$sigma
        mu_mult_ = const_list[[j]]$mu_mult
        mu = mu_mult_ %*% current_liabil[-j] # Udelukker den liability vi er kommet til


        if (j == 1) {
          current_liabil[j] = rnorm(1, mu, sqrt(sigma))
        }
        else {
          crit <- qnorm(1-K)
          crit_bound = pnorm(crit, mu, sqrt(sigma))

          interval_2 = c(0, crit_bound)
          interval_1 = c(crit_bound, 1)

          phen = phenos[j-1]
          interval = phen * interval_1 + (1-phen)*interval_2
          U = runif(1,interval[1], interval[2])
          current_liabil[j] = qnorm(U, mu, sqrt(sigma))
        }


        if(s > start_run){
          liabil[i+start_run*(p-2),] = current_liabil
        }
      }
    }

    if(s>(start_run+min_run)){
      end_val = s-start_run
      if(sd(liabil[1:end_val,1])/sqrt(length(liabil[1:end_val,1])) < 0.1){
        if (all_est){
          return(liabil[1:end_val,])
        }
        return(colMeans(liabil[1:end_val,]))

      }
    }

  }
}

