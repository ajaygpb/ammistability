## `ammistability`: Additive Main Effects and Multiplicative Interaction Model Stability Parameters ![logo](https://raw.githubusercontent.com/ajaygpb/ammistability/master/inst/extdata/ammistability.png)

###### Version : [0.1.4.9000](https://ajaygpb.github.io/ammistability/); Copyright (C) 2017-2025: [ICAR-DGR](https://en.wikipedia.org/wiki/Directorate_of_Groundnut_Research); License: [GPL-2\|GPL-3](https://www.r-project.org/Licenses/)

##### *Ajay, B. C.¹, Aravind, J.² and Abdul Fiyaz, R³*

1.  RRS, ICAR-Directorate of Groundnut Research, Anantapur.
2.  ICAR-National Bureau of Plant Genetic Resources, New Delhi.
3.  ICAR-Indian Institute of Rice Research, Hyderabad.

------------------------------------------------------------------------

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

[![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)](https://cran.r-project.org/web/checks/check_results_ammistability.html)

|                                   |                                                                                                                                                                                                                        |
|:----------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| r-devel-linux-x86_64-debian-clang | [![CRAN check - r-devel-linux-x86_64-debian-clang](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-debian-clang/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html) |
| r-devel-linux-x86_64-debian-gcc   | [![CRAN check - r-devel-linux-x86_64-debian-gcc](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-debian-gcc/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html)     |
| r-devel-linux-x86_64-fedora-clang | [![CRAN check - r-devel-linux-x86_64-fedora-clang](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-fedora-clang/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html) |
| r-devel-linux-x86_64-fedora-gcc   | [![CRAN check - r-devel-linux-x86_64-fedora-gcc](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-fedora-gcc/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html)     |
| r-patched-linux-x86_64            | [![CRAN check - r-patched-linux-x86_64](https://badges.cranchecks.info/flavor/r-patched-linux-x86_64/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html)                       |
| r-release-linux-x86_64            | [![CRAN check - r-release-linux-x86_64](https://badges.cranchecks.info/flavor/r-release-linux-x86_64/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html)                       |

[![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white)](https://cran.r-project.org/web/checks/check_results_ammistability.html)

|                          |                                                                                                                                                                                                      |
|:-------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| r-devel-windows-x86_64   | [![CRAN check - r-devel-windows-x86_64](https://badges.cranchecks.info/flavor/r-devel-windows-x86_64/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html)     |
| r-release-windows-x86_64 | [![CRAN check - r-release-windows-x86_64](https://badges.cranchecks.info/flavor/r-release-windows-x86_64/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html) |
| r-oldrel-windows-x86_64  | [![CRAN check - r-oldrel-windows-x86_64](https://badges.cranchecks.info/flavor/r-oldrel-windows-x86_64/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html)   |

[![MacOS](https://img.shields.io/badge/mac%20os-000000?style=for-the-badge&logo=apple&logoColor=white)](https://cran.r-project.org/web/checks/check_results_ammistability.html)

|                        |                                                                                                                                                                                                  |
|:-----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| r-release-macos-x86_64 | [![CRAN check - r-release-macos-x86_64](https://badges.cranchecks.info/flavor/r-release-macos-x86_64/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html) |
| r-oldrel-macos-x86_64  | [![CRAN check - r-oldrel-macos-x86_64](https://badges.cranchecks.info/flavor/r-oldrel-macos-x86_64/ammistability.svg)](https://cran.r-project.org/web/checks/check_results_ammistability.html)   |

## Citing `ammistability`

To cite the methods in the package use:

``` r
citation("ammistability")
```

``` R
To cite the R package 'ammistability' in publications use:

  Ajay, B. C., Aravind, J., and Abdul Fiyaz, R. (2019). ammistability:
  R package for ranking genotypes based on stability parameters derived
  from AMMI model. Indian Journal of Genetics and Plant Breeding (The),
  79(2), 460-466.
  https://www.isgpb.org/article/ammistability-r-package-for-ranking-genotypes-based-on-stability-parameters-derived-from-ammi-model

  Ajay, B. C., Aravind, J., and Abdul Fiyaz, R. (2025).  ammistability:
  Additive Main Effects and Multiplicative Interaction Model Stability
  Parameters. R package version 0.1.4.9000,
  https://ajaygpb.github.io/ammistability/,
  https://CRAN.R-project.org/package=ammistability.

This free and open-source software implements academic research by the
authors and co-workers. If you use it, please support the project by
citing the package.

To see these entries in BibTeX format, use 'print(<citation>,
bibtex=TRUE)', 'toBibtex(.)', or set
'options(citation.bibtex.max=999)'.
```
