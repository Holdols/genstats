## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(genstats)
library(bigsnpr)
library(dplyr)

genetic_data = snp_attach('genetic_data.rds')
G=genetic_data$genotypes
y =  genetic_data$fam %>% select('pheno_0') 
beta = genetic_data$map$beta


## -----------------------------------------------------------------------------
gwas_summary = GWAS(G=G, y=y[[1]], p=0.05/100000, logreg = FALSE, ncores=3)
gwas_summary

## -----------------------------------------------------------------------------
manhattan_plot(gwas_summary, beta, c(0.05/100000, 0.05))

## -----------------------------------------------------------------------------
scatter_plot(gwas_summary, beta)

## -----------------------------------------------------------------------------
power_plot(gwas_summary, beta)

