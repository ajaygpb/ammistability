
## `ammistability`: Additive Main Effects and Multiplicative Interaction Model Stability Parameters <img src="https://raw.githubusercontent.com/ajaygpb/ammistability/master/inst/extdata/ammistability.png" align="right" alt="logo" width="173" height = "200" style = "border: none; float: right;">

###### Version : [0.0.0.9000](https://ajaygpb.github.io/ammistability/); License: [GPL-2|GPL-3](https://www.r-project.org/Licenses/)

##### *Ajay, B. C.<sup>1</sup>, Aravind, J.<sup>2</sup> and Abdul Fiyaz, R.<sup>3</sup>*

1.  RRS, ICAR-Directorate of Groundnut Research, Anantapur
2.  ICAR-National Bureau of Plant Genetic Resources, New Delhi
3.  ICAR-Indian Institute of Rice Research, Hyderabad

-----

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.0.2-6666ff.svg)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version-last-release/ammistability)](https://cran.r-project.org/package=ammistability)
<!-- [![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/ammistability?color=green)](https://CRAN.R-project.org/package=ammistability) -->
<!-- [![packageversion](https://img.shields.io/badge/Package%20version-0.2.3.3-orange.svg)](https://github.com/ajaygpb/ammistability) -->
[![develVersion](https://img.shields.io/badge/devel%20version-0.0.0.9000-orange.svg)](https://github.com/ajaygpb/ammistability)
<!-- [![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/ajaygpb/ammistability/total.svg)] -->
[![Project Status:
WIP](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Last-changedate](https://img.shields.io/badge/last%20change-2018--07--30-yellowgreen.svg)](/commits/master)
<!-- [![Rdoc](http://www.rdocumentation.org/badges/version/ammistability)](http://www.rdocumentation.org/packages/ammistability) -->
<!-- [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.841963.svg)](https://doi.org/10.5281/zenodo.841963) -->
[![Analytics](https://pro-pulsar-193905.appspot.com/UA-123032895-2/welcome-page)](https://github.com/aravind-j/google-analytics-beacon)

-----

## Description

Computes various stability parameters from Additive Main Effects and
Multiplicative Interaction (AMMI) analysis results such as Modified AMMI
Stability Value (MASV), Sums of the Absolute Value of the Interaction
Principal Component Scores (SIPC), Sum Across Environments of
Genotype-Environment Interaction Modelled by AMMI (AMGE), Sum Across
Environments of Absolute Value of Genotype-Environment Interaction
Modelled by AMMI (AV\_(AMGE)), AMMI Stability Index (ASI), Modified ASI
(MASI), AMMI Based Stability Parameter (ASTAB), Annicchiarico’s D
Parameter (DA), Zhang’s D Parameter (DZ), Averages of the Squared
Eigenvector Values (EV), Stability Measure Based on Fitted AMMI Model
(FA), Absolute Value of the Relative Contribution of IPCs to the
Interaction (Za). Further calculates the Simultaneous Selection Index
for Yield and Stability from the computed stability parameters. See the
vignette for complete list of citations for the methods implemented.

## Installation

### Install development version from Github

``` r
devtools::install_github("ajaygpb/ammistability")
```

## Detailed tutorial

For a detailed tutorial on how to used this package type:

``` r
browseVignettes(package = 'ammistability')
```

## What’s new

To know whats new in this version type:

``` r
news(package='ammistability')
```

## Links

<!-- [CRAN page](https://cran.r-project.org/package=ammistability) -->

[Github page](https://github.com/ajaygpb/ammistability)

<!-- [Github website](https://ajaygpb.github.io/ammistability/) -->

<!-- [Zenodo DOI](https://doi.org/10.5281/zenodo.1310011) -->

## Citing `ammistability`

To cite the methods in the package use:

``` r
citation("ammistability")
```

``` 

To cite the R package 'ammistability' in publications use:

  Ajay B. C., J. Aravind, and R. Abdul Fiyaz (2018).
  ammistability: Additive Main Effects and Multiplicative
  Interaction Model Stability Parameters. R package version
  0.0.0.9000, https://ajaygpb.github.io/ammistability/.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {ammistability: Additive Main Effects and Multiplicative Interaction Model Stability Parameters},
    author = {{Ajay Basapura Chandrashekar} and {J. Aravind} and {R. Abdul Fiyaz}},
    year = {2018},
    note = {R package version 0.0.0.9000},
    note = {https://ajaygpb.github.io/ammistability/},
  }

This free and open-source software implements academic research by
the authors and co-workers. If you use it, please support the
project by citing the package.
```
