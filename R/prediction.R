#'  Uses cross validation to train and estimate PRS for each fold
#'
#' The funcion returns list of the PRS estimated of each fold for each threshold. The function will train on the remaining data and test for each fold.
#' The function uses GWAS to estimate the casual SNP and uses these the find the PRS. The out will be a list containing matrices with a column for each threshold.
#' The interpretation of this treshold is that a given effect of a SNP will not be included if the -log10 transformation of the p-value is is smaller than the threshold.
#' Fx. a value of 3 will correspond to the p value being smaller than 0.001.
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
    } else {
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




#'  Trains model on train data given the target vector.
#'
#' The funcion will first use GWAS to find the casual effect of each SNP.
#' It will then calculate the PRS and use a linear or logistic regression to estimate the phenotypes.
#' The final output is list containing the model from GWAS and the regression. This can be used to predict the phenotypes of test data.
#' @param train_data List generated from gen_sim.
#' @param y The target vector. Could either be estimated liabilities from LTFH or phenotypes.
#' @param thr Treshold for p-value to be used in calculating PRS.
#' @param LogReg Boolean indicating if logistic regression should be used to estimate the casual effect.
#' @return List containing output from GWAS and Linear regression of PRS on phenotype of the subject.
#' @importFrom magrittr "%>%"
#' @export
pred_model = function(train_data, y, thr, LogReg_g = FALSE, LogReg_prs = FALSE){
  if (!LogReg_g) {gwas <- bigstatsr::big_univLinReg(train_data$genotypes, y.train = y)}
  else {gwas <- bigstatsr::big_univLogReg(train_data$genotypes, y.train = y)}
  prs_ <- bigstatsr::snp_PRS(G = train_data$genotypes, betas.keep =gwas$estim, lpS.keep = -predict(gwas), thr.list = thr)
  prs_ = prs_[,1]


  if (!LogReg_prs) {model = lm(train_data$fam$pheno_0 ~ prs_)}
  else {model <- glm(train_data$fam$pheno_0 ~ prs_, family = binomial(link = "probit"))}

  return(list('gwas'=gwas, 'model_prs'=model))
}



#' Predicts the probability of having the given trait
#'
#' The function will take outcome from pred_model() and use it an to find the estimated probability of each subject having the trait/sickness.
#' @param test_data List generated from gen_sim.
#' @param model Output from pred_model().
#' @param thr Treshold for p-value to be used in calculating PRS.
#' @return Vector of estimated probabilities.
#' @importFrom magrittr "%>%"
#' @export
prediction = function(test_data, model, thr) {
  gwas = model$gwas
  model_prs = model$model_prs

  prs_ = snp_PRS(G = test_data$genotypes, betas.keep = gwas$estim, lpS.keep = -predict(gwas), thr.list = thr)
  probs <- predict(model_prs, data.frame('prs_' = prs_[,1]), type = "response")
  return(probs)
}
