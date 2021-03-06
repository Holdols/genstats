% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_LTFH.R
\name{estimate_conf}
\alias{estimate_conf}
\title{Estimates liabilities for every configuration in data}
\usage{
estimate_conf(unique_comb, n_sib = 0, K = 0.5, h2 = 0.5)
}
\arguments{
\item{unique_comb}{Configuration of phenotypes to estimate genetic liability for.}

\item{n_sib}{Amount of siblings.}

\item{K}{The prevalance of trait.}

\item{h2}{The heritability of trait.}
}
\value{
Tibble with estimated liabilities for each configuration.
}
\description{
Helper function for LTFH()
}
