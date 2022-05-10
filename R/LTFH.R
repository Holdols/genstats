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
#' Create LTFH estimations for a configuration.
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

  # det er muligvis en præference ting, men jeg er ikke så stor fan af at hardcode jeres burn-in på den her måde.
  # jeg synes egentligt i skulle sætte den som argument i LTFH.
  # i kunne faktisk også overveje bare at ændre navnet på denne funktion, da den, i al sin simpelhed, bare sampler
  # fra en trunkeret flerdimensionel normal fordeling. Når i så estimere den genetiske liability osv,
  # kan i kalde den funktion LTFH
  start_run = 500
  s = 0
  p = 0

  while(TRUE) {
    p = p +1 # jeg er ikke lige helt med på formålet med p.
    for (i in 2:start_run) {
      s = s + 1
      for (j in 1:k) {
        # Udregner parametre
        sigma = const_list[[j]]$sigma
        mu_mult_ = const_list[[j]]$mu_mult
        mu = mu_mult_ %*% current_liabil[-j] # Udelukker den liability vi er kommet til


        if (j == 1) {
          current_liabil[j] = rnorm(1, mu, sqrt(sigma))
        } else { #jeg er ikke sikker på at den læser dette ordentligt, hvis } afslutter if uden else er efter.
          crit <- qnorm(1-K)
          crit_bound = pnorm(crit, mu, sqrt(sigma))

          interval_2 = c(0, crit_bound)
          interval_1 = c(crit_bound, 1)

          phen = phenos[j-1] # i burde nok lave et tjek for at sikre jer, at denne værdi er en numerisk værdi og ikke en 1x1 matrix
          interval = phen * interval_1 + (1-phen)*interval_2
          U = runif(1,interval[1], interval[2])
          current_liabil[j] = qnorm(U, mu, sqrt(sigma))
        }


        if(s > start_run){
          liabil[i+start_run*(p-2),] = current_liabil
        }
      }
    }

    # hvis jeg forstår dette korrekt, så kunne i nok have gjort det lettere for jer selv ved bare at få
    # alle værdier gemt i liabil, og så bruge jeres start_run værdi til at fjerne det antal fra begyndelsen.
    # med fx liabil[-(1:start_run),]



    if(s>(start_run+min_run)){
      print(s)
      end_val = s-start_run
      if(sd(liabil[1:end_val,1])/sqrt(length(liabil[1:end_val,1])) < 0.1){ # fx sæt all(sapply(1:ncol(liabil), function(i) sd(liabil[,i])) / sqrt(end_val) < 0.01)
        if (all_est){
          return(liabil[1:end_val,])
        }
        return(colMeans(liabil[1:end_val,]))

      }
    }
    #i returnere alle indgange, men tjekker kun én. Det skal i passe lidt på med. I ved jo ikke om
    #nogle af de andre indgange er meget ustabile stadigvæk.

  }
  # post kig: nu er jeg mere med på hvad i laver. Jeg synes sådan set i har fat i den lange end, men jeg synes også, at i virker til at gøre det mere besværligt end det behøver at være.
  # fx behøver i ikke tjekke om værdierne er efter burn-in, bare hav en burn-in periode, og så fjern dem senere (når i tjekke konvergens og returnere)
  # I vil kunne lave et while loop på fx  all(sapply(1:ncol(liabil), function(i) sd(liabil[,i])) / sqrt(end_val) < 0.01),
  # og så bare fylde op i en liste med nye batches. brug her do.call ind i bind_row eller lignende inden i tjekker værdierne. så skal i ikke på forhånd gøre jer nogen tanker om
  # i ender ud over jeres predefineret dimensioner.
  # at vokse en liste i et loop er ikke et problem på samme måde som at vokse en vektor eller matrix (da en masse kopiering sker internt i R her, men det sker ikke for lister)
}


