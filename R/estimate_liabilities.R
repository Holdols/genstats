#'  Using phentoypes the function estimates full liabilities using LT-FH
#'
#' @param pheno matrix containing phenotypes
#' @param K The prevalence of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Tibble containing phenotypes and estimates of the individual and family
#' @export
estimate_liabilities = function(pheno, K=0.05, n_sib=0, h2=0.5){
  unique_comb = unique(pheno)
  estimated_liabil = estimate_conf(unique_comb, h2, n_sib, K)
  key = get_names(c("phen"),n_sib)
  full_data = left_join(as_tibble(pheno),as_tibble(estimated_liabil),by=key)
  order = get_names(c("phen", 'est_l_g'), n_sib)
  return(as_tibble(full_data) %>% relocate(order))
}
