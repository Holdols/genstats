
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genstats

<!-- badges: start -->
<!-- badges: end -->

genstats provide easy access to tools for statistical analysis of
genetic data. Genetic data is difficult to come by, but provides a very
informative basis for statistical use. We have therefore through some
assumptions provided method for simulating genetic data which can be
used for both association analysis and predictive statistics. The
package provides standard methods for association analysis such as GWAS
and LT-FH. Though the structure of simulated data is easy to use in
comnination with dplyr and it is therefor easy to implement your own
methods.

The packages uses the method of LT-FH as described by Hujoel, M.L.A.,
Gazal, S., Loh, PR. et al.. Link to article is found in Biblography.

The package is developed as a part of our bachelor’s program in data
science at Aarhus University.

## Installation

You can install the development version of genstats from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Holdols/genstats")
```

## Biblography

-   Hujoel, M.L.A., Gazal, S., Loh, PR. et al. Liability threshold
    modeling of case–control status and family history of disease
    increases association power. Nat Genet 52, 541–547 (2020).
    <https://doi.org/10.1038/s41588-020-0613-6>"

The package uses functions and methods from both bigsnpr:
<https://privefl.github.io/bigsnpr/> and bigstatsr:
<https://privefl.github.io/bigstatsr/index.html>
