

#' Estimates liabilities for every configuration in data
#'
#' @param unique_comb Combinations to estimate
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Tibble with estimated liabilities
estimate_conf = function(unique_comb, K=0.5, h2=0.5, n_sib=0){
  k = ncol(unique_comb)
  covmat = get_cov(h2, n_sib)


  estimates = matrix(NA, nrow = nrow(unique_comb), ncol=4+n_sib)
  for (i in 1:nrow(unique_comb)){
    estimates[i,] = LTFH(covmat, unique_comb[i,], K=K)
  }

  out = unique_comb %>% bind_cols('l_g_est_0'=estimates[,1])

  return(out)
}



#' Estimates liabilities for every subject
#'
#' @param pheno Phenotypes
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Estimated liabilities
#' @export
LTFH = function(data, n_sib=0, K=0.05, h2=0.5){

  unique_comb = data %>% select(., contains('pheno')) %>%
    relocate(order = get_names(c("phen"), n_sib)) %>%
    unique(.)


  estimated_liabil = estimate_conf(unique_comb, K, n_sib, h2)
  full_data = left_join(data, estimated_liabil,
                        by=get_names(c("phen"),n_sib))
  return(as_tibble(full_data))
}


