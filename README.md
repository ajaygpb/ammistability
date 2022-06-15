
## `ammistability`: Additive Main Effects and Multiplicative Interaction Model Stability Parameters <img src="https://raw.githubusercontent.com/ajaygpb/ammistability/master/inst/extdata/ammistability.png" align="right" alt="logo" width="173" height = "200" style = "padding: 10px; border: none; float: right;">

###### Version : [0.1.2.9000](https://ajaygpb.github.io/ammistability/); Copyright (C) 2017-2022: [ICAR-DGR](http://www.dgr.org.in/); License: [GPL-2\|GPL-3](https://www.r-project.org/Licenses/)

##### *Ajay, B. C.<sup>1</sup>, Aravind, J.<sup>2</sup> and Abdul Fiyaz, R<sup>3</sup>*

1.  RRS, ICAR-Directorate of Groundnut Research, Anantapur.
2.  ICAR-National Bureau of Plant Genetic Resources, New Delhi.
3.  ICAR-Indian Institute of Rice Research, Hyderabad.

------------------------------------------------------------------------

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.0.2-6666ff.svg?logo=R)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/ammistability)](https://cran.r-project.org/package=ammistability)
[![Dependencies](https://tinyverse.netlify.com/badge/ammistability)](https://cran.r-project.org/package=ammistability)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/ammistability?color=green)](https://CRAN.R-project.org/package=ammistability)
[![develVersion](https://img.shields.io/badge/devel%20version-0.1.2.9000-orange.svg)](https://github.com/ajaygpb/ammistability)
[![Github Code
Size](https://img.shields.io/github/languages/code-size/ajaygpb/ammistability.svg)](https://github.com/ajaygpb/ammistability)
[![R-CMD-check](https://github.com/ajaygpb/ammistability/workflows/R-CMD-check/badge.svg)](https://github.com/ajaygpb/ammistability/actions)
[![Project Status:
Inactive](http://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Last-changedate](https://img.shields.io/badge/last%20change-2022--06--16-yellowgreen.svg)](https://github.com/ajaygpb/ammistability/commits/master)
[![Zenodo
DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1344756.svg)](https://doi.org/10.5281/zenodo.1344756)
[![Website -
pkgdown](https://img.shields.io/website-up-down-green-red/https/ajaygpb.github.io/ammistability.svg)](https://ajaygpb.github.io/ammistability/)
[![.](https://pro-pulsar-193905.appspot.com/UA-123032895-2/welcome-page)](https://github.com/aravind-j/google-analytics-beacon)
<!-- [![Rdoc](https://www.rdocumentation.org/badges/version/ammistability)](https://www.rdocumentation.org/packages/ammistability) -->
<!-- [![packageversion](https://img.shields.io/badge/Package%20version-0.2.3.3-orange.svg)](https://github.com/ajaygpb/ammistability) -->
<!-- [![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/ajaygpb/ammistability/total.svg)] -->

------------------------------------------------------------------------

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

The package can be installed from CRAN as follows:

``` r
# Install from CRAN
install.packages('ammistability', dependencies=TRUE)
```

The development version can be installed from github as follows:

``` r
# Install development version from Github
devtools::install_github("ajaygpb/ammistability")
```

## Detailed tutorial

For a detailed tutorial (vignette) on how to used this package type:

``` r
browseVignettes(package = 'ammistability')
```

The vignette for the latest version is also available
[online](https://ajaygpb.github.io/ammistability/articles/Introduction.html).

## What’s new

To know whats new in this version type:

``` r
news(package='ammistability')
```

## Links

[CRAN page](https://cran.r-project.org/package=ammistability)

[Github page](https://github.com/ajaygpb/ammistability)

[Documentation website](https://ajaygpb.github.io/ammistability/)

[Zenodo DOI](https://doi.org/10.5281/zenodo.1344756)

## CRAN checks

<table class="table table-striped table-hover" style="width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
Flavour
</th>
<th style="text-align:left;">
CRAN check
</th>
</tr>
</thead>
<tbody>
<tr grouplength="6">
<td colspan="2" style="border-bottom: 1px solid;">
<strong>[![Linux](https://shields.io/badge/Linux--9cf?logo=Linux&style=social)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-devel-linux-x86_64-debian-clang
</td>
<td style="text-align:left;">
[![CRAN check -
r-devel-linux-x86_64-debian-clang](https://cranchecks.info/badges/flavor/r-devel-linux-x86_64-debian-clang/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-devel-linux-x86_64-debian-gcc
</td>
<td style="text-align:left;">
[![CRAN check -
r-devel-linux-x86_64-debian-gcc](https://cranchecks.info/badges/flavor/r-devel-linux-x86_64-debian-gcc/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-devel-linux-x86_64-fedora-clang
</td>
<td style="text-align:left;">
[![CRAN check -
r-devel-linux-x86_64-fedora-clang](https://cranchecks.info/badges/flavor/r-devel-linux-x86_64-fedora-clang/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-devel-linux-x86_64-fedora-gcc
</td>
<td style="text-align:left;">
[![CRAN check -
r-devel-linux-x86_64-fedora-gcc](https://cranchecks.info/badges/flavor/r-devel-linux-x86_64-fedora-gcc/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-patched-linux-x86_64
</td>
<td style="text-align:left;">
[![CRAN check -
r-patched-linux-x86_64](https://cranchecks.info/badges/flavor/r-patched-linux-x86_64/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-release-linux-x86_64
</td>
<td style="text-align:left;">
[![CRAN check -
r-release-linux-x86_64](https://cranchecks.info/badges/flavor/r-release-linux-x86_64/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr grouplength="1">
<td colspan="2" style="border-bottom: 1px solid;">
<strong>[![Solaris](https://shields.io/badge/Solaris--9cf?logo=Oracle&style=social)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-patched-solaris-x86
</td>
<td style="text-align:left;">
[![CRAN check -
r-patched-solaris-x86](https://cranchecks.info/badges/flavor/r-patched-solaris-x86/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr grouplength="3">
<td colspan="2" style="border-bottom: 1px solid;">
<strong>[![Windows](https://shields.io/badge/Windows--9cf?logo=Windows&style=social)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-devel-windows-ix86+x86_64
</td>
<td style="text-align:left;">
[![CRAN check -
r-devel-windows-ix86+x86_64](https://cranchecks.info/badges/flavor/r-devel-windows-ix86+x86_64/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-release-windows-ix86+x86_64
</td>
<td style="text-align:left;">
[![CRAN check -
r-release-windows-ix86+x86_64](https://cranchecks.info/badges/flavor/r-release-windows-ix86+x86_64/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-oldrel-windows-ix86+x86_64
</td>
<td style="text-align:left;">
[![CRAN check -
r-oldrel-windows-ix86+x86_64](https://cranchecks.info/badges/flavor/r-oldrel-windows-ix86+x86_64/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr grouplength="2">
<td colspan="2" style="border-bottom: 1px solid;">
<strong>[![MacOS](https://shields.io/badge/MacOS--9cf?logo=Apple&style=social)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-release-macos-x86_64
</td>
<td style="text-align:left;">
[![CRAN check -
r-release-macos-x86_64](https://cranchecks.info/badges/flavor/r-release-macos-x86_64/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
r-oldrel-macos-x86_64
</td>
<td style="text-align:left;">
[![CRAN check -
r-oldrel-macos-x86_64](https://cranchecks.info/badges/flavor/r-oldrel-macos-x86_64/ammistability)](https://cran.r-project.org/web/checks/check_results_ajaygpb_at_yahoo.co.in.html)
</td>
</tr>
</tbody>
</table>

## Citing `ammistability`

To cite the methods in the package use:

``` r
citation("ammistability")
```


    To cite the R package 'ammistability' in publications use:

      Ajay, B. C., Aravind, J., and Abdul Fiyaz, R. (2019). ammistability: R package for ranking genotypes
      based on stability parameters derived from AMMI model. Indian Journal of Genetics and Plant Breeding
      (The), 79(2), 460-466.
      http://www.isgpb.org/article/ammistability-r-package-for-ranking-genotypes-based-on-stability-parameters-derived-from-ammi-model

      Ajay, B. C., Aravind, J., and Abdul Fiyaz, R. (NA).  ammistability: Additive Main Effects and
      Multiplicative Interaction Model Stability Parameters. R package version 0.1.2.9000,
      https://ajaygpb.github.io/ammistability/, https://CRAN.R-project.org/package=ammistability.

    This free and open-source software implements academic research by the authors and co-workers. If you use
    it, please support the project by citing the package.
    To see these entries in BibTeX format, use 'print(<citation>, bibtex=TRUE)', 'toBibtex(.)', or set
    'options(citation.bibtex.max=999)'.
