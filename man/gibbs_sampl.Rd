% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbs_sampler.R
\name{gibbs_sampl}
\alias{gibbs_sampl}
\title{Create LTFH estimations for a configuration.}
\usage{
gibbs_sampl(
  covmat,
  phenos,
  K = 0.05,
  s_val = 0,
  start_run = 500,
  all_est = FALSE
)
}
\arguments{
\item{covmat}{The covariance matrix.}

\item{phenos}{A binary vector containing the phenotype for each family member
of the form c(p_subject, p_parent1, p_parent2, p_sibling1, ... ,p_siblingN)
where p_familymember is a binary value (1 or 2).}

\item{K}{The prevalance of trait. If False no truncation is applied.}

\item{s_val}{The starting value of liabilities.}

\item{start_run}{Number of iterations before convergence is expected.}

\item{all_est}{If TRUE return the value for each iteration after burn in, else return mean of values.}
}
\value{
A vector containing LTFH estimate of liabilities of the form
c(genetic_liability_subject, liability_subject, liability_parent1,
liability_parent2, liability_sibling1, ..., liability_siblingN).
}
\description{
This function estimates the genetic liability for a subject given the configuration of each familiy members phenotype.
The the phenotype is given by either a 0 or 1 indicating if the family member has the trait/sickness.
The function uses Monte Carlo integration to estimate the conditional distribution.
}
\examples{
gibbs_sampl(get_cov(0.5, n_sib = 1), c(1, 1, 0, 0))
}
