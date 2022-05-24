

#' Estimates liabilities for every configuration in data
#'
#' @param unique_comb Configuration of phenotypes to estimate genetic liability for.
#' @param n_sib Amount of siblings.
#' @param K The prevalance of trait.
#' @param h2 The heritability of trait.
#' @return Tibble with estimated liabilities for each configuration.
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
#' @param data List generated from gen_sim.
#' @param n_sib Amount of siblings.
#' @param K The prevalance of trait.
#' @param h2 The heritability of trait.
#' @return Tibble containing data$fam and estimated liabilities.
#' @importFrom magrittr "%>%"
#' @export
LTFH = function(data, n_sib=0, K=0.05, h2=0.5){
  fam = data$fam
  stopifnot(tibble::is_tibble(fam) || is.data.frame(fam))
  stopifnot(is.double(n_sib), n_sib >= 0, n_sib%%1==0)
  stopifnot(is.double(K), K < 1 || K > 0)
  stopifnot(is.double(h2), h2 < 1 || h2 > 0)


  unique_comb = fam %>% dplyr::select(., dplyr::contains('pheno')) %>%
    dplyr::relocate(get_names(c("pheno"), n_sib)) %>%
    unique(.)


  estimated_liabil = estimate_conf(unique_comb, n_sib, K, h2)


  full_data = dplyr::left_join(fam, estimated_liabil,
                        by=get_names(c("pheno"),n_sib))
  return(full_data)
}


