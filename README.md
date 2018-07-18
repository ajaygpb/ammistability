
## `AMMIStbP`: Additive Main Effects and Multiplicative Interaction Model Stability Parameters

###### Version : [0.0.0.9000](https://ajaygpb.github.io/AMMIStbP/); License: [GPL-2|GPL-3](https://www.r-project.org/Licenses/)

##### *Ajay, B. C.<sup>1</sup>, Aravind, J.<sup>2</sup>, Abdul Fiyaz, R.<sup>3</sup> and S. K. Bera<sup>4</sup>*

1.  RRS, ICAR-Directorate of Groundnut Research, Anantapur
2.  ICAR-National Bureau of Plant Genetic Resources, New Delhi
3.  ICAR-Indian Institute of Rice Research, Hyderabad
4.  ICAR-Directorate of Groundnut Research, Junagadh

-----

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.0.2-6666ff.svg)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version-last-release/AMMIStbp)](https://cran.r-project.org/package=AMMIStbp)
<!-- [![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/AMMIStbp?color=green)](https://CRAN.R-project.org/package=AMMIStbp) -->
<!-- [![packageversion](https://img.shields.io/badge/Package%20version-0.2.3.3-orange.svg)](https://github.com/ajaygpb/AMMIStbp) -->
[![develVersion](https://img.shields.io/badge/devel%20version-0.0.0.9000-orange.svg)](https://github.com/ajaygpb/AMMIStbp)
<!-- [![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/ajaygpb/AMMIStbp/total.svg)] -->
[![Project Status:
WIP](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Last-changedate](https://img.shields.io/badge/last%20change-2018--07--19-yellowgreen.svg)](/commits/master)
<!-- [![Rdoc](http://www.rdocumentation.org/badges/version/AMMIStbp)](http://www.rdocumentation.org/packages/AMMIStbp) -->
<!-- [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.841963.svg)](https://doi.org/10.5281/zenodo.841963) -->
<!-- [![Analytics](https://pro-pulsar-193905.appspot.com/UA-116683292-1/welcome-page)](https://github.com/aravind-j/google-analytics-beacon) -->

-----

## Warning

Package is not stable at present. It is under development. Please use
only after release of stable version.

## Description

Provides functions to compute various germination indices such as
germinability, median germination time, mean germination time, mean
germination rate, speed of germination, Timson’s index, germination
value, coefficient of uniformity of germination, uncertainty of
germination process, synchrony of germination etc. from germination
count data. Includes functions for fitting cumulative seed germination
curves using four-parameter hill function and computation of associated
parameters. See the vignette for more, including full list of citations
for the methods implemented.

## Installation

### Install development version from Github

``` r
devtools::install_github("ajaygpb/AMMIStbP")
```

## Detailed tutorial

For a detailed tutorial on how to used this package type:

``` r
browseVignettes(package = 'AMMIStbP')
```

## What’s new

To know whats new in this version type:

``` r
news(package='AMMIStbP')
```

## Links

<!-- [CRAN page](https://cran.r-project.org/package=AMMIStbP) -->

[Github page](https://github.com/ajaygpb/AMMIStbP)

<!-- [Github website](https://ajaygpb.github.io/AMMIStbP/) -->

<!-- [Zenodo DOI](https://doi.org/10.5281/zenodo.1310011) -->

## Citing `AMMIStbP`

To cite the methods in the package
    use:

``` r
citation("AMMIStbP")
```

    Warning in citation("AMMIStbP"): no date field in DESCRIPTION file of
    package 'AMMIStbP'
    Warning in citation("AMMIStbP"): could not determine year for 'AMMIStbP'
    from package DESCRIPTION file
    
    To cite package 'AMMIStbP' in publications use:
    
      Ajay Basapura Chandrashekar, J. Aravind, R. Abdul Fiyaz and S.
      K. Bera (NA). AMMIStbP: Additive Main Effects and Multiplicative
      Interaction Model Stability Parameters. R package version
      0.0.0.9000. http://github.com/ajaygpb/AMMIStbp
    
    A BibTeX entry for LaTeX users is
    
      @Manual{,
        title = {AMMIStbP: Additive Main Effects and Multiplicative Interaction Model Stability Parameters},
        author = {{Ajay Basapura Chandrashekar} and {J. Aravind} and {R. Abdul Fiyaz} and {S. K. Bera}},
        note = {R package version 0.0.0.9000},
        url = {http://github.com/ajaygpb/AMMIStbp},
      }
