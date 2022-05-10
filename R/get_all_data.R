#'  Using genetic and enviromental liabilities the function calculates phenotypes, full liabilites and estimates full liabilities using LT-FH
#'
#' @param data Tibble containing information about each indiviual and family's genetic and enviromental liabilities
#' @param K The prevalence of trait
#' @param n_sib Amount of siblings
#' @param h2 The heritability
#' @return Tibble containing all information and estimates of the individual and family
#' @export
get_all_data = function(data, K=0.05, n_sib=0, h2=0.5){
  full_liabil = calc_full_liabil(data, n_sib) # i kunne bare udregne jeres fulde liabtilities tidligere (fx når i simulere genotyperne), så er dette ikke nødvendigt.
  pheno = calc_phen(data, n_sib, K) # det samme gælder her. Så er det hele også ordnet under paralleliseringen
  estimate_liabil = estimate_liabilities(pheno, K, n_sib, h2)
  order = get_names(c('l', 'l_g', 'l_e', "phen", 'est_l_g'), n_sib)
  full_data = bind_cols(as_tibble(full_liabil), estimate_liabil, data) %>% relocate(order)
  return(full_data)
}
