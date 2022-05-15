#'  Plotting estimated liabilities agianst distribution
#'
#' @param pheno matrix containing phenotypes
#' @param h2 The heritability
#' @return Plot of estimated liabilities
#' @export
control_plot = function(phenos, h2, col="black"){
  covmat = get_cov(h2)
  estimates = LTFH(covmat, phenos, all_est = TRUE)
  sd = sd(estimates[,1])
  mean = mean(estimates[,1])


  ggplot(data.frame(x = c(-2, 3)), aes(x = x)) +
    ggplot2::stat_function(fun = dnorm, args = list(mean = mean, sd = sd), mapping = aes(x),color=col) +
    ggplot2::stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(h2)), mapping = aes(x)) +
    ggplot2::geom_vline(xintercept=mean, linetype="dashed",color=col) +
    ggplot2::geom_text(aes(x=mean+0.2, label=round(mean,2), y=0), colour=col)
}

