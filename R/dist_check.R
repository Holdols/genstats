#' Check distribution of liabilities
#'
#' The function calculates the mean and variance of each full liability returned by gen_sim().
#' The assumption which the simulation is based on is that each full liability will be normaly distributed.
#' The function therefore also returns the qqplot of each full liability. The purpose of this function is to provide
#' certainty that the simulated data are correct.
#' @param data List generated from gen_sim.
#' @return Creates a qqplot, mean and sample variance for the liability of each family.
#' @examples
#' dist_check(snp_attach("test.rds"))
#' @export
dist_check = function(data){
  l = data$fam %>% dplyr::select(., contains('l_f'))
  par(mfrow=c(1,ncol(l)))
  for (j in l){
    qqnorm(j)
  }
  c = 0
  for (j in l){
    out = c(mean(j), var(j))
    names(out) = c(sprintf("mean_l%s", c), sprintf("variance_l%s", c))
    print(out)
    c = c + 1
  }
}
