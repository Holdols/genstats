#'  Uses cross validation to train and estimate PRS for each fold
#'
#' @param data List generated from gen_sim.
#' @param y The target vector. Could either be estimated liabilities from LTFH or phenotypes.
#' @param cross_folds Number of folds in cross validation.
#' @param LogReg Boolean indicating if logistic regression should be used to estimate the casual effect.
#' @return List with estimated PRS for each fold.
#' @importFrom magrittr "%>%"
#' @export
PRS_cross <- function(data, y01, cross_folds, LogReg = FALSE){
  folds = list()
  G = data$genotypes
  indexes <- rows_along(G)
  test_size = length(indexes)%/%cross_folds

  for (i in 1:cross_folds){
    folds[[i]] <- sort(sample(indexes, test_size))
    indexes <- setdiff(indexes, folds[[i]])
  }

  out = list()
  N = nrow(G)
  for (i in 1:cross_folds){


    ind_test <- folds[[i]]
    ind_train <- setdiff(rows_along(G), ind_test)
    if (LogReg == TRUE){
      bigstatsr::big_univLogRegg(G, y.train = y01[ind_train], ind.train = ind_train)
    }
    else {
      gwas_train <- bigstatsr::big_univLinReg(G, y.train = y01[ind_train], ind.train = ind_train)
    }

    pval <- -predict(gwas_train)
    pval_thrs <- seq(0, 4, by = 0.5)
    prs <- bigsnpr::snp_PRS(G, betas.keep = gwas_train$estim,
                            ind.test = ind_test,
                            lpS.keep = pval,
                            thr.list = pval_thrs)
    out[[i]] = prs
  }
  return(out)
}
