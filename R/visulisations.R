#'  Plotting estimated liabilities agianst distribution
#'
#' A recreation of figure 1. b in Hujoel et al. (2020)
#'
#' @param pheno Vector containing phenotypes.
#' @param h2 The heritability of trait.
#' @return Plot of estimated liabilities.
#' @importFrom magrittr "%>%"
#' @export
control_plot = function(phenos, h2, col="black"){
  covmat = get_cov(h2)
  estimates = gibbs_sampl(covmat, phenos, all_est = TRUE)
  sd = sd(estimates[,1])
  mean = mean(estimates[,1])


  ggplot2::ggplot(data.frame(x = c(-2, 3)), ggplot2::aes(x = x)) +
    ggplot2::stat_function(fun = dnorm, args = list(mean = mean, sd = sd), mapping = ggplot2::aes(x),color=col) +
    ggplot2::stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(h2)), mapping = ggplot2::aes(x)) +
    ggplot2::geom_vline(xintercept=mean, linetype="dashed",color=col) +
    ggplot2::geom_text(ggplot2::aes(x=mean+0.2, label=round(mean,2), y=0), colour=col)
}



#' Creates evaluation plot for MSE, AUC or r squared
#'
#'The function creates one plot for each fold in PRS.
#'This plot will show the chosen method evaluated for each threshold.
#'Lastly the function will plot the mean for each threshold.
#'The function can be used to find which threshold. of p-values best can describe the data.
#' @param PRS Output from PRS_cross.
#' @param data List generated from gen_sim.
#' @param method Indicates Which method each fold should be evaluated on.
#' @return Plots method for each threshold for each fold and mean of all folds.
#' @importFrom magrittr "%>%"
#' @export
prs_plot = function(PRS, data, method='MSE'){
  pval_thrs = seq(0, 4, by = 0.5)

  df = lapply(1:length(PRS), function(i) {
    targ = data$fam$pheno_0[as.numeric(row.names(PRS[[i]]))]

    if (method=='AUC') {eval = apply(PRS[[i]], 2, AUC, target = targ)}
    else if (method=='R2') {
      eval = apply(PRS[[i]], 2, function(pred,target) summary(lm(target~pred))$r.squared,  target = targ)}
    else {eval = apply(PRS[[i]], 2, function(pred,target) mean((pred-target)^2), target = targ)}

    tibble('Eval'=eval, 'Threshold'=pval_thrs, 'fold' = paste0('Fold ', i))
  }) %>% bind_rows

  print(df %>%
          ggplot(ggplot2::aes(y=Eval, x=Threshold)) +
          ggplot2::geom_point() +
          ggplot2::geom_line()+
          ggplot2::xlab('Threshold for p value') +
          ggplot2::ylab(method) +
          ggplot2::facet_wrap(~fold))

  df %>%
    dplyr::group_by(Threshold) %>%
    dplyr::summarise('mean_Eval'=mean(Eval), 'sd'=sd(Eval)) %>%
    ggplot2::ggplot(ggplot2::aes(y=mean_Eval, x=Threshold)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::xlab('Threshold for p value') +
    ggplot2::ylab(paste0('Mean ', method)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_Eval-sd, ymax=mean_Eval+sd), width=.2, position=position_dodge(0.05))
}



#' Creates a manhattan plot
#'
#' This function returns a manhattan plot that vizualises the which
#'
#' @param gwas_summary Output from GWAS.
#' @param beta A vector containing the actual casual effect of each SNP.
#' @return A manhattan plot.
#' @importFrom magrittr "%>%"
#' @export
manhattan_plot = function(gwas_summary, beta, thresholds = FALSE){
  plot = tibble::tibble('p_vals' = gwas_summary$p_vals) %>%
    dplyr::mutate('causal' = (beta!=0)-0) %>%
    ggplot2::ggplot(ggplot2::aes(x=(1:length(p_vals)), y=-log10(p_vals))) +
    ggplot2::geom_point(ggplot2::aes(colour=as.character(causal))) +
    ggplot2::ylab('-log10(P-value)') +
    ggplot2::xlab('Chromosome') +
    ggplot2::labs(color='True casual SNPs')

  if (any(thresholds != FALSE)){
    plot = plot + ggplot2::geom_hline(yintercept = -log10(thresholds), linetype='dashed')
  }
  plot
}


#' Plots true causal effects agianst estimated effects
#'
#' @param gwas_summary Output from GWAS.
#' @param beta A vector containing the actual casual effect of each SNP.
#' @return Scatterplot of true causal effects agianst estimated effects
#' @importFrom magrittr "%>%"
#' @export
scatter_plot = function(gwas_summary, beta){
  plot = tibble::tibble('beta'= beta, 'estim'=gwas_summary$estim, 'causal_est'=gwas_summary$causal_estimate) %>%
    dplyr::mutate('causal' = (beta!=0)-0) %>%
    ggplot2::ggplot(ggplot2::aes(x=beta, y=estim)) +
    ggplot2::geom_point(ggplot2::aes(colour=as.character(causal_est))) +
    ggplot2::geom_abline(slope = 1) +
    ggplot2::labs(color='Estimated to have casual effect') +
    ggplot2::ylab('Estimated effect') +
    ggplot2::xlab('Actual effect')

  plot
}



#' Plot the convergence of gibb_sampl
#'
#' @param ests Output from gibb_sampl with all_est=TRUE.
#' @return Plots showing values at each iteration.
#' @importFrom magrittr "%>%"
#' @export
plot_gibbs <- function(ests){
  means <- tibble::as_tibble(ests) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "pheno") %>%
    dplyr::group_by(pheno) %>%
    dplyr::summarise(MN = mean(value), MN_x = length(value)/2, .groups = "keep")

  tibble::as_tibble(ests) %>%
    dplyr::mutate(iterations = (1:length(ests[,1]))) %>%
    tidyr::pivot_longer(!iterations, names_to = "pheno") %>%
    ggplot2::ggplot(ggplot2::aes(x = iterations, y = value)) +
    ggplot2::geom_line(ggplot2::aes(col = "brick")) +
    ggplot2::geom_hline(data = means, ggplot2::aes(yintercept = MN)) +
    ggplot2::facet_wrap(~pheno) +
    ggplot2::geom_label(data = means, ggplot2::aes(x = MN_x, label=round(MN,digits = 3), y=MN -0.1), size = 2.5, vjust = "top") +
    ggplot2::theme(legend.position="none")
}

#' Creates a power plot
#'
#' @param gwas_summary Output from GWAS.
#' @param beta A vector containing the actual casual effect of each SNP.
#' @return A power plot.
#' @importFrom magrittr "%>%"
#' @export
power_plot <- function(gwas_summary, beta){
  tibble::as_tibble(gwas_summary) %>%
    dplyr::mutate(true_causal = (beta != 0) - 0) %>%
    dplyr::filter(true_causal == 1) %>%
    dplyr::arrange(abs(estim)) %>%
    dplyr::mutate(power = cumsum(causal_estimate)/sum(true_causal)) %>%
    ggplot2::ggplot(ggplot2::aes(x = estim, y = power)) + ggplot2::geom_line()
}




#' Plots estimated liabilities agianst true liabilities
#'
#' @param LTFH_est Output from LTFH.
#' @param gaps Varaible for distance between labels
#' @return Plot of estimated liabilities agianst true liabilities.
#' @importFrom magrittr "%>%"
#' @export
LTFH_plot = function(LTFH_est, gaps=1){
  new_df = LTFH_est %>% mutate('conf' = apply(dplyr::select(., contains('phen')), 1, paste0, collapse=', '))

  plot = new_df %>%
    ggplot2::ggplot(ggplot2::aes(x=l_g_est_0, y=l_g_0)) +
    ggplot2::geom_point(ggplot2::aes(color=conf)) +
    ggplot2::geom_abline(slope=1)  +
    guides(color = 'none') + xlab('Estimated liabilities') + ylab('True liabilities')

  min_df = new_df %>% group_by(conf) %>%
    summarise('min_' = min(l_g_0), 'mean' = mean(l_g_est_0), n = n()) %>%
    mutate('sequence' = seq(0,gaps, gaps/(length(conf)-1)))

  plot + geom_text(data=min_df,
                   ggplot2::aes(x=mean+sequence/5, y=min(min_)+sequence-0.5, label=conf, color=conf),  size = 3) +
    geom_segment(data=min_df, ggplot2::aes(x=mean+sequence/5-0.09, y=min(min_)+sequence-0.49, xend=mean, yend=min_, color=conf))

}




#' Creates evaluation plot for decision boundary
#'
#'The function creates one plot for each given bound. It shows the confusion matrix for mean of each bound.
#'The function can be used to find which boundary for PRS best describes the data.
#' @param train_data List generated from gen_sim.
#' @param y The target vector. Could either be estimated liabilities from LTFH or phenotypes.
#' @param cross_folds Number of folds in cross validation.
#' @param bounds Decision boundaries to plot outcome for.
#' @param thr Treshold for p-value to be used in calculating PRS.
#' @param ncores Amount of cores to be used.
#' @param LogReg Boolean indicating if logistic regression should be used to estimate the casual effect.
#' @return List containing output from GWAS and Linear regression of PRS on phenotype of the subject. It
#' @importFrom magrittr "%>%"
#' @export
decision_cross <- function(train_data, y, cross_folds, bounds, thr, ncores = 1, LogReg = FALSE){
  folds = list()
  G = train_data$genotypes
  target = train_data$fam$pheno_0
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
      gwas_train <- bigstatsr::big_univLogReg(G, y01.train = y[ind_train], ind.train = ind_train, ncores = ncores)
    } else {
      gwas_train <- bigstatsr::big_univLinReg(G, y.train = y[ind_train], ind.train = ind_train, ncores = ncores)
    }

    pval <- -predict(gwas_train)
    prs <- bigsnpr::snp_PRS(G, betas.keep = gwas_train$estim,
                            ind.test = ind_test,
                            lpS.keep = pval,
                            thr.list = thr)



    true = target[ind_test]
    temp = list()

    out[[i]] = lapply(1:length(bounds), function(j) {
      pred = (prs>bounds[j])-0

      n = c(sum(pred == 1 & true == 1), sum(pred == 1 & true ==0), sum(pred == 0 & true == 1), sum(pred == 0 & true ==0))


      confusion_matrix = tibble::tibble('Predicted'=as.character(c(1,1,0,0)), 'Actual'=as.character(c(1,0,1,0)), n)


      info = tibble::tibble(bound = bounds[j], fold = paste0('fold ', i))
      dplyr::bind_cols(confusion_matrix, info)
    }) %>% dplyr::bind_rows(.)
  }



  temp = out %>% dplyr::bind_rows(.) %>% dplyr::group_by(bound, Predicted, Actual) %>% dplyr::summarise('mean_n'=mean(n), .groups='keep')

  plt = temp %>% dplyr::bind_rows(.) %>% ggplot2::ggplot(mapping = ggplot2::aes(x = Predicted, y = Actual)) +
    ggplot2::geom_tile(ggplot2::aes(fill = mean_n), show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%1.0f", mean_n)), vjust = 1) +
    ggplot2::scale_fill_gradient(high = "firebrick", low = 'dodgerblue3', trans='pseudo_log') +
    ggplot2::facet_wrap(~bound)

  print(plt)
}
