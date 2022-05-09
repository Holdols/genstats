

#' Estimates liabilities for every configuration in data
#'
#' @param unique_comb Combinations to estimate
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Tibble with estimated liabilities
#' @noRd
estimate_conf = function(unique_comb, K=0.5, h2=0.5, n_sib=0){
  k = ncol(unique_comb)
  covmat = get_cov(h2, n_sib)


  estimates = matrix(NA, nrow = nrow(unique_comb), ncol=4+n_sib)
  for (i in 1:nrow(unique_comb)){
    estimates[i,] = LTFH(covmat, unique_comb[i,], K=K)
  }
  col_names = get_names(c("phen"), n_sib)
  pheno = as_tibble(unique_comb)
  colnames(pheno) = col_names
  esti = as_tibble(estimates[,1])
  colnames(esti) = 'l_g_est_0'

  return(bind_cols(pheno, esti))
}



#' Estimates liabilities for every subject
#'
#' @param pheno Phenotypes
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Estimated liabilities
#' @export
estimate_liabilities = function(pheno, K=0.05, n_sib=0, h2=0.5){

  unique_comb = unique(pheno)
  estimated_liabil = estimate_conf(unique_comb, K, n_sib, h2)
  key = get_names(c("phen"),n_sib)
  full_data = left_join(as_tibble(pheno),as_tibble(estimated_liabil),by=key)
  return(as_tibble(full_data))
}



#' Creates tibble containing information and estimates for every subject
#'
#' @param data Liabilities of subject and family
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Tibble with data and estimates
#' @export
get_all_data = function(data, K=0.05, n_sib=0, h2=0.5){
  full_liabil = calc_full_liabil(data, n_sib)
  pheno = calc_phen(data, n_sib, K)
  estimate_liabil = estimate_liabilities(pheno, K, n_sib, h2)
  order = get_names(c('l', 'l_g', 'l_e', "phen"), n_sib)
  print(order)
  full_data = bind_cols(as_tibble(full_liabil), estimate_liabil, data) %>% relocate(order)

  return(full_data)
}
