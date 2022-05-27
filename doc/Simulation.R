## -----------------------------------------------------------------------------
library(genstats)
library(bigsnpr)

## ----eval=FALSE---------------------------------------------------------------
#  gen_sim('simple_sim_example', N=1000, M=1000, C=50, block_size = 100, fam=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  simple_data = snp_attach('simple_sim_example.rds')
#  simple_data

## ----eval=FALSE---------------------------------------------------------------
#  simple_data$G[(1:5),(1:5)]
#  simple_data$fam
#  simple_data$map

## ----eval=FALSE---------------------------------------------------------------
#  dist_check(simple_data)

## ----eval=FALSE---------------------------------------------------------------
#  gen_sim('family_sim_example', N=1000, M=1000, C=50, block_size = 100, n_sib=2)
#  family_data = snp_attach('family_sim_example.rds')
#  family_data$fam

## ----eval=FALSE---------------------------------------------------------------
#  dist_check(simple_data)

