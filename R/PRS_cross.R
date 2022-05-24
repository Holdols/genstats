PRS_cross <- function(data, y01, cross_folds){
  folds = list()
  G = data$genotypes
  indexes <- rows_along(G)
  test_size = length(indexes)%/%cross_folds

  for (i in 1:cross_folds){
    folds[[i]] <- sort(sample(indexes, test_size))
    indexes <- setdiff(indexes, folds[[i]])
  }

  out = list()
  for (i in 1:cross_folds){

    N = nrow(G)
    ind_test <- folds[[i]]
    ind_train <- setdiff(rows_along(G), ind_test)

    gwas_train <- bigstatsr::big_univLinReg(G, y.train = y01[ind_train], ind.train = ind_train)
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
