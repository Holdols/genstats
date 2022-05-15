

#' Estimates liabilities for every configuration in data
#'
#' @param unique_comb Combinations to estimate
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Tibble with estimated liabilities
#' @importFrom magrittr "%>%"
estimate_conf = function(unique_comb, n_sib=0, K=0.5, h2=0.5){
  k = ncol(unique_comb)
  covmat = get_cov(h2, n_sib)
  estimates = matrix(NA, nrow = nrow(unique_comb), ncol=4+n_sib)


  for (i in 1:nrow(unique_comb)){
    estimates[i,] = gibbs_sampl(covmat, unlist(unique_comb[i,], use.names = FALSE), K)
  }

  out = unique_comb %>% dplyr::bind_cols('l_g_est_0'=estimates[,1])

  return(out)
}



#' Estimates liabilities for every subject
#'
#' @param pheno Phenotypes
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Estimated liabilities
#' @importFrom magrittr "%>%"
#' @export
LTFH = function(data, n_sib=0, K=0.05, h2=0.5){


  unique_comb = data %>% dplyr::select(., dplyr::contains('pheno')) %>%
    dplyr::relocate(get_names(c("pheno"), n_sib)) %>%
    unique(.)


  estimated_liabil = estimate_conf(unique_comb, n_sib, K, h2)


  full_data = dplyr::left_join(data, estimated_liabil,
                        by=get_names(c("pheno"),n_sib))
  return(full_data)
}


