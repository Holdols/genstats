#'  Plotting estimated liabilities agianst distribution
#'
#' @param pheno Vector containing phenotypes
#' @param h2 The heritability of trait
#' @return Plot of estimated liabilities
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


manhatten_plot = function(summary, beta, thresholds = FALSE){
  plot = tibble('p_vals' = summary$p_vals) %>% mutate('causal' = (beta!=0)-0) %>%
    ggplot(aes(x=(1:length(p_vals)), y=-log10(p_vals))) +
    geom_point(aes(colour=as.character(causal))) +
    ylab('-log10(P-value)') + xlab('Chromosome') + labs(color='True casual SNPs')
  if (any(thresholds != FALSE)){
    plot = plot + geom_hline(yintercept = -log10(thresholds), linetype='dashed')
  }
  plot
}

scatter_plot = function(summary, beta){
  plot = tibble('beta'= beta, 'estim'=summary$estim, 'causal_est'=summary$causal_estimate) %>%
    mutate('causal' = (beta!=0)-0) %>% ggplot(aes(x=beta, y=estim)) + geom_point(aes(colour=as.character(causal_est))) + geom_abline(slope = 1) + labs(color='Estimated to have casual effect') + yab('Estimated effect') + xlab('Actual effect')
  plot
}

#' Plot the iterations of gibb_sampl
#'
#' @param ests output from gibb_sampl
#' @return plots showing values at each iteration
#' @export
plot_gibbs <- function(ests){
  colnames(ests) = c(c("l_g_0"), get_names("l_f", n_sib = ncol(ests)-4))
  means <- tibble::as_tibble(ests) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "pheno") %>%
    dplyr::group_by(pheno) %>%
    dpyr::summarise(MN = mean(value), MN_x = length(value)/2, .groups = "keep")

  tibble::as_tibble(ests) %>%
    dplyr::mutate(iterations = (1:length(l_g_0))) %>%
    tidyr::pivot_longer(!iterations, names_to = "pheno") %>%
    ggplot2::ggplot(ggplot2::aes(x = iterations, y = value)) +
    ggplot2::geom_line(ggplot2::aes(col = "brick")) +
    ggplot2::geom_hline(data = means, ggplot2::aes(yintercept = MN)) +
    ggplot2::facet_wrap(~pheno) +
    ggplot2::geom_label(data = means, ggplot2::aes(x = MN_x, label=round(MN,digits = 3), y=MN -0.1), size = 2.5, vjust = "top")
}

#' Create power plot
#'
#' @param gwas_summary output from GWAS
#' @param beta A vector containing the actual casual effect of each SNP.
#' @return A power plot
#' @export
power_plot <- function(gwas_summary, beta){
  tibble::as_tibble(gwas_df) %>%
    dplyr::mutate(true_causal = (beta != 0) - 0) %>%
    dplyr::filter(true_causal == 1) %>%
    dplyr::arrange(abs(estim)) %>%
    dplyr::mutate(power = cumsum(causal_estimate)/100) %>%
    ggplot2::ggplot(ggplot2::aes(x = estim, y = power)) + ggplot2::geom_line()
}
