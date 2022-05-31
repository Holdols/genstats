
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genstats

<!-- badges: start -->
<!-- badges: end -->

genstats provide easy access to tools for statistical analysis of
genetic data. Genetic data is difficult to come by, but provides a very
informative basis for statistical anaysis. We have therefore through
some assumptions provided method for simulating genetic data which can
be used for both association analysis and predictive statistics. The
package provides standard methods for association analysis such as GWAS
and LT-FH. Furthermore, we have provided a method which can be used to
predict the phenotypes of individuals, based on the their genotypes. The
structure of simulated data is easy to use in combination with dplyr and
it is therefore easy to implement your own methods. This also means that
the functions are very dependent on column names. It is therefore
important if this package is used with data not simulated by genstats,
that the columns match the following column names.

``` r
genstats::get_names(c('l', 'l_g', 'pheno'), n_sib = 1)
#>  [1] "l_0"      "l_g_0"    "pheno_0"  "l_p1"     "l_g_p1"   "pheno_p1"
#>  [7] "l_p2"     "l_g_p2"   "pheno_p2" "l_s3"     "l_g_s3"   "pheno_s3"
```

With the amount of siblings being variable from 1 to what the computer
can handle.

Examples and details about each method can be found in the articles. Due
to a large amount of data, many of the plots are included as images.
Source code for all plots can be found in
<https://github.com/Holdols/genstats/blob/main/src_plots.Rmd>.

The packages uses the method of LT-FH as described by Hujoel, M.L.A.,
Gazal, S., Loh, PR. et al.. Link to article is found in Bibliography.

The package is developed as a part of our bachelor’s program in data
science at Aarhus University.

We would like to thank our supervisor Emil Michael Pedersen for patience
and great supervision!

## Installation

You can install the development version of genstats from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Holdols/genstats")
```

## Bibliography

-   Hujoel, M.L.A., Gazal, S., Loh, PR. et al. Liability threshold
    modeling of case–control status and family history of disease
    increases association power. Nat Genet 52, 541–547 (2020).
    <https://doi.org/10.1038/s41588-020-0613-6>"

The package uses functions and methods from both bigsnpr:
<https://privefl.github.io/bigsnpr/> and bigstatsr:
<https://privefl.github.io/bigstatsr/index.html>

Please visit those site for more information about how data is stored
