#'  Plotting estimated liabilities agianst distribution
#'
#' @param pheno matrix containing phenotypes
#' @param h2 The heritability
#' @return Plot of estimated liabilities
#' @export
control_plot = function(phenos, h2, col="black"){
  covmat = get_cov(h2)
  estimates = LTFH(covmat, phenos, all_est = TRUE, min_run = 8000)
  sd = sd(estimates[,1])
  mean = mean(estimates[,1])


  ggplot(data.frame(x = c(-2, 3)), aes(x = x)) +
    stat_function(fun = dnorm, args = list(mean = mean, sd = sd), mapping = aes(x),color=col) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(h2)), mapping = aes(x)) +
    geom_vline(xintercept=mean, linetype="dashed",color=col) +
    geom_text(aes(x=mean+0.2, label=round(mean,2), y=0), colour=col)
}

