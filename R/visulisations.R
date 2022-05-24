#'  Plotting estimated liabilities agianst distribution
#'
#' @param pheno matrix containing phenotypes
#' @param h2 The heritability
#' @return Plot of estimated liabilities
#' @export
control_plot = function(phenos, h2, col="black"){
  covmat = get_cov(h2)
  estimates = gibbs_sampl(covmat, phenos, all_est = TRUE)
  sd = sd(estimates[,1])
  mean = mean(estimates[,1])


  ggplot(data.frame(x = c(-2, 3)), aes(x = x)) +
    ggplot2::stat_function(fun = dnorm, args = list(mean = mean, sd = sd), mapping = aes(x),color=col) +
    ggplot2::stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(h2)), mapping = aes(x)) +
    ggplot2::geom_vline(xintercept=mean, linetype="dashed",color=col) +
    ggplot2::geom_text(aes(x=mean+0.2, label=round(mean,2), y=0), colour=col)
}


MSE_prs_plot = function(PRS, data){
  pval_thrs = seq(0, 4, by = 0.5)
  lapply(1:length(test), function(i) {
    targ = data$fam$pheno_0[as.numeric(row.names(te))]
    mse = apply(test[[i]], 2 , function(pred,target) mean((pred-target)^2), target = targ)
    tibble('MSE'=mse, 'Threshold'=pval_thrs, 'fold' = paste0('Fold ', i))
  }) %>% bind_rows %>% ggplot(aes(y=MSE, x=Threshold)) + geom_point() + geom_line()+ xlab('Threshold for p value') + facet_wrap(~fold)
}


AUC_prs_plot = function(PRS, data){
  pval_thrs = seq(0, 4, by = 0.5)
  lapply(1:length(test), function(i) {
    targ = data$fam$pheno_0[as.numeric(row.names(te))]
    auc = apply(test[[i]], 2 , AUC, target = targ)
    tibble('AUC'=auc, 'Threshold'=pval_thrs, 'fold' = paste0('Fold ', i))
  }) %>% bind_rows %>% ggplot(aes(y=AUC, x=Threshold)) + geom_point() + geom_line()+ xlab('Threshold for p value') + facet_wrap(~fold)
}
