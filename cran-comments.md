 Version 0.1.3 - First submission

* Fixed issue with re-building of vignette outputs and running of examples.
* As the references are too many, instead of in DESCRIPTION all the citations are added in the package vignette with the doi or isbn as specified in the last line of DESCRIPTION.

### Test environments
* local Windows 10 Pro 21H2, R-release (R 4.2.0) & R-devel (R 4.3.0 Pre-release).
* local Ubuntu 20.04, R-release (R 4.2.0) & R-devel (R 4.3.0 Pre-release).
* win-builder, R-release (R 4.2.0) & R-devel (R 4.3.0 Pre-release).
* rhub:macos-highsierra-release-cran - x86_64-apple-darwin17.0 (64-bit), R-release (R 4.2.0).

### R CMD check results
* There were no ERRORs, WARNINGs.
* 1 NOTE - Possibly misspelled words in DESCRIPTION: AMGE, AMMI, ASI, ASTAB, Annicchiarico's, EV, Genotype, IPCs, MASI, MASV, SIPC, Za, Zhang's. These are false positives and may be ignored. They are either abbreviations or names of authors in citations.
