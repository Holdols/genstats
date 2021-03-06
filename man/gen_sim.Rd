% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_sim.R
\name{gen_sim}
\alias{gen_sim}
\title{Simulate genetic data}
\usage{
gen_sim(
  filename,
  beta = NULL,
  MAF = NULL,
  N = 1e+05,
  M = 1e+05,
  n_sib = 0,
  K = 0.05,
  h2 = 0.5,
  C = 1000,
  block_size = 1000,
  fam = TRUE,
  parallel_plan = "multisession"
)
}
\arguments{
\item{filename}{Filename for file backed matrix(FBM) and rds file.}

\item{beta}{Vector containing casual effect of each SNP. If NULL a vector will be simulated.}

\item{MAF}{Vector containing minor allele frequencies (MAF). If NULL a vector will be simulated.}

\item{N}{Amount of subjects.}

\item{M}{Amount of SNPs.}

\item{K}{The prevalance of trait.}

\item{h2}{The heritability of trait.}

\item{C}{Amount of causal SNP's.}

\item{block_size}{Size of FBM to be processed in each iteration.}

\item{fam}{Boolean deciding if simulation should include a family structure.}

\item{parallel_plan}{Plan for parallelization. See ?future::plan.}

\item{n_sibs}{Amount of siblings.}
}
\value{
rds file referring to a list containing FBM with genotypes, generated SNP info and information about subject and family.
}
\description{
This function is used to simulate genetic data. Using parallelization the function can simulate genotypes,
liabilities both with and without family structure. The function will return a list where all information about
the simulation can be found. This includes which SNP are casual and what their effect is, the genetic and full liabilities
og subjects, parents and siblings and their phenotypes. This object can be used to further develop or test statistical analysis
on genetic data.
}
