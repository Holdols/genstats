## -----------------------------------------------------------------------------
library(genstats)
library(bigsnpr)
library(tidyverse)

## -----------------------------------------------------------------------------
gridExtra::grid.arrange(control_plot(phenos = c(0,0,0), h2 = 0.5, col = "royalblue"),
                        control_plot(phenos = c(0,0,1), h2 = 0.5, col = "mediumorchid"),
                        control_plot(phenos = c(1,1,1), h2 = 0.5, col = "firebrick"), heights = 0.5, name = c("Plots for configuration (0,0,0), (0,0,1), (1,1,1)"))

## -----------------------------------------------------------------------------
# The rds file containing simulated data is loaded
genetic_data = snp_attach('genetic_data.rds')

# We save the genotypes and beta for later use
G = genetic_data$genotypes
beta = genetic_data$map$beta


## -----------------------------------------------------------------------------
LTFH_est = LTFH(genetic_data, n_sib = 2)
head(LTFH_est)

## -----------------------------------------------------------------------------
LTFH_summary = GWAS(G=G, y=LTFH_est$l_g_0, p=0.05/100000, logreg = FALSE, ncores=3)
head(LTFH_summary)

## -----------------------------------------------------------------------------
manhattan_plot(LTFH_summary, beta, c(0.05/1000, 0.05))

## -----------------------------------------------------------------------------
scatter_plot(LTFH_summary, beta)

## -----------------------------------------------------------------------------
power_plot(LTFH_summary, beta)

