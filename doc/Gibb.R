## -----------------------------------------------------------------------------
library(genstats)

covariance = get_cov(0.5, n_sib = 1)
estimated = gibbs_sampl(covariance, phenos=c(1,0,1,0), K=FALSE, all_est = TRUE)
round(cov(estimated), digits = 2)

## -----------------------------------------------------------------------------
ests = gibbs_sampl(covmat = get_cov(0.5, n_sib = 2), phenos = c(1,1,0,1,0), K = 0.05, start_run = 500, all_est = TRUE)
colnames(ests) = c("Genetic Liability subject",
                   "Liability subject, CASE",
                   "Liability parent1, CASE",
                   "Liability parent2, CONTROL",
                   "Liability sibling1, CASE",
                   "Liability sibling2, CONTROL")
# means <- as_tibble(ests) %>%
#           pivot_longer(cols = everything(), names_to = "pheno") %>%
#             group_by(pheno) %>%
#               summarise(MN = mean(value), MN_x = length(value)/2, .groups = "keep")
# 
# as_tibble(ests) %>%
#   mutate(iterations = (1:length(l_g_0))) %>%
#   pivot_longer(!iterations, names_to = "pheno") %>%
#   ggplot(aes(x = iterations, y = value)) + geom_line(aes(col = "brick")) + geom_hline(data = means, aes(yintercept = MN)) + facet_wrap(~pheno) + geom_label(data = means, aes(x = MN_x, label=round(MN,digits = 3), y=MN -0.1), size = 2.5, vjust = "top")
plot_gibbs(ests)

