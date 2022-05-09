#' Calculate phenotypes
#'
#' @param l_e_0 vector containing enviromental liabilities
#' @param l_g_0 vector containing genetic liabilities
#' @param k prevalance of trait
#' @return binary vector contatining phenotypes
#' @export
y_func = function(l_e_0, l_g_0, K=0.05) {
  l_0 = l_g_0 + l_e_0
  T_ = qnorm(1-K)
  y = ifelse(l_0>=T_, 1, 0)
  return(y)
}
