
#' Calculate phentypes from from envirromental and genetic liabilities
#'
#' @param liabilities Liabilities of subject and family
#' @param K The prevalance of trait
#' @param n_sib Amount of siblings
#' @return Matrix containing phenotypes for each subject
#' @export
calc_phen = function(liabilities, K=0.05, n_sib=0){
  T_ = qnorm(1-K)
  N = nrow(liabilities)
  k = 3+n_sib
  phen = matrix(0, nrow=N, ncol=k)
  search = get_names(c(""), n_sib)
  for (i in 1:k){
    phen[,(i)] = rowSums(select(liabilities, contains(search[i]))) >=T_
  }

  col_names=get_names(c("phen"), n_sib)
  colnames(phen) = col_names

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
    full_l[,(i)] = rowSums(select(liabilities, contains(search[i])))
  }

  col_names=get_names(c("l"), n_sib)
  colnames(full_l) = col_names

  return(full_l)
}
