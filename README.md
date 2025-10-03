
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rPSAQ

<!-- badges: start -->

<!-- badges: end -->

The goal of rPSAQ is to provide an R package to deconvoluate PSAQ+1
(Protein Standard Absolute Quantification) signals from MS2
isotopologues abundances. The package uses the MS2 fragments abundances
to estimate the respective contributions of PSAQ and endogenous peptides
and calculate the abundance ratios between PSAQ and endogenous peptides,
based on fragment intensities or areas.

The package provides a sample MS2 abundances dataset named
`sample_10ng_r4`

## Installation

You can install the development version of rPSAQ from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("edyp-lab/rPSAQ")
```

## Example

This is a basic example which shows you deconvolute signals from the
example dataset:

``` r
library(rPSAQ)
## basic example code using the provided dataset
deconvoluted_data = deconvolute_peptides_abundances(sample_10ng_r4)
```

The computed ratios of each fragment can be graphically represented:

``` r
# plot the previously deconvoluted data
plot_psaq_ratios(deconvoluted_data)
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Removed 4 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```

<img src="man/figures/README-figure-1.png" width="100%" />

To write the result, use the openxlsx package

``` r
openxlsx::write.xlsx(deconvoluted_data, "outputfile.xlsx")
```

The deconvolution can be performed from a list of Excel files stored in
a single folder. For example to perform deconvolution on Excel files
starting by “7843” and ending with “.xlsx” from the the “inst/extdata”
folder:

``` r
batch_psaq_analysis("inst/extdata", "^7843.*\\.xlsx$")
```
