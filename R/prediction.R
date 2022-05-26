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
  Else {gwas <- bigstatsr::big_univLogReg(train_data$genotypes, y.train = y)}
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
