
#' Calculate phentypes from from environmental and genetic liabilities
#'
#' @param liabilities Liabilities of subject and family
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @return Matrix containing phenotypes for each subject
#' @export
calc_phen = function(liabilities, K=0.05, n_sib=0){
  T_ = qnorm(1-K)
  k = 3+n_sib
  phen = matrix(0, nrow = nrow(liabilities), ncol = k)
  search = get_names(c(""), n_sib)
  for (i in 1:k){
    phen[,(i)] = rowSums(select(liabilities, contains(search[i]))) >=T_ # hvorfor laver i rowsums? har i environment og genetic liabs? eller genetic og full liabs ? umiddelbart burde det være nok bare at se på full liabs til status
  }
  #måske jeg er lidt for dplyr minded for tiden, men jeg ville nok have haft den fulde liab for hver individ i en tibble. Her er "full" en substring, som alle individer har for an angive at det er den fulde liab
  # så fx p1_full og p1_gen for at adskille de to liabs:

  # liabilities = mvtnorm::rmvnorm(n = 1000, mean = rep(0, k), sigma = covmat)
  # colnames(liabilities) = c("child_gen", paste0(c("child", "p1", "p2", "s1"), "_full"))
  # T_ = qnorm(1-0.05)
  # as_tibble(liabilities) %>% mutate(across(.cols = contains("full"),
  #                                          .fns = ~ .x > T_,
  #                                          .names = "{stringr::str_replace(.col, 'full','status')}"))

  colnames(phen) = get_names(c("phen"), n_sib)

  return(phen)
}


#' Calculate full liabilities from envirromental and genetic liabilities
#'
#' @param liabilities Liabilities of subject and family
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @return Matrix containing full lliabilities for each subject
#' @export
calc_full_liabil = function(liabilities, n_sib=0){
  N = nrow(liabilities)
  k = 3+n_sib
  full_l = matrix(0, nrow=N, ncol=k)
  search = get_names(c(""), n_sib)
  for (i in 1:k){
    full_l[,(i)] = rowSums(select(liabilities, contains(search[i]))) # i kunne vel bare udregne den fulde liabs mens i laver det genetiske data. Så skal i ikke tænker på at have en fkt til det.
  }

  col_names=get_names(c("l"), n_sib)
  colnames(full_l) = col_names

  return(full_l)
}
