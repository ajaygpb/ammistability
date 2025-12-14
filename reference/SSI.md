# Simultaneous Selection Indices for Yield and Stability

`SSI` computes the Simultaneous Selection Index for Yield and Stability
(SSI) according to the methods specified in the argument `method`.

## Usage

``` r
SSI(y, sp, gen, method = c("farshadfar", "rao"), a = 1)
```

## Arguments

- y:

  A numeric vector of the mean yield/performance of genotypes.

- sp:

  A numeric vector of the stability parameter/index of the genotypes.

- gen:

  A character vector of the names of the genotypes.

- method:

  The method for the computation of simultaneous selection index. Either
  `"farshadfar"` or `"rao"` (See **Details**).

- a:

  The ratio of the weights given to the stability components for
  computation of SSI when `method = "rao"` (See **Details**).

## Value

A data frame with the following columns:

- SP:

  The stability parameter values.

- SSI:

  The computed values of simultaneous selection index for yield and
  stability.

- rSP:

  The ranks of the stability parameter.

- rY:

  The ranks of the mean yield of genotypes.

- means:

  The mean yield of the genotypes.

The names of the genotypes are indicated as the row names of the data
frame.

## Details

The SSI according to Rao and Prabhakaran (2005) (\\I\_{i}\\) is computed
as follows:

\\I\_{i} = \frac{\overline{Y}\_{i}}{\overline{Y}\_{..}} + \alpha
\frac{\frac{1}{SP\_{i}}}{\frac{1}{T}\sum\_{i=1}^{T}\frac{1}{SP\_{i}}}\\

Where \\SP\_{i}\\ is the stability measure of the \\i\\th genotype under
AMMI procedure; \\\overline{Y}\_{i}\\ is mean performance of \\i\\th
genotype; \\\overline{Y}\_{..}\\ is the overall mean; \\T\\ is the
number of genotypes under test and \\\alpha\\ is the ratio of the
weights given to the stability components (\\w\_{2}\\) and yield
(\\w\_{1}\\) with a restriction that \\w\_{1} + w\_{2} = 1\\. The
weights can be specified as required.

|                |                |                |
|----------------|----------------|----------------|
| **\\\alpha\\** | **\\w\_{1}\\** | **\\w\_{2}\\** |
| 1.00           | 0.5            | 0.5            |
| 0.67           | 0.6            | 0.4            |
| 0.43           | 0.7            | 0.3            |
| 0.25           | 0.8            | 0.2            |

The SSI proposed by Farshadfar (2008) is called the Genotype stability
index (\\GSI\\) or Yield stability index (\\YSI\\) (Farshadfar et al.
2011) and is computed by summation of the ranks of the stability
index/parameter and the ranks of the mean yields.

\\GSI = YSI = R\_{SP} + R\_{Y}\\

Where, \\R\_{SP}\\ is the stability parameter/index rank of the genotype
and \\R\_{Y}\\ is the mean yield rank of the genotype.

## References

Farshadfar E (2008). “Incorporation of AMMI stability value and grain
yield in a single non-parametric index (GSI) in bread wheat.” *Pakistan
Journal of biological sciences*, **11**(14), 1791.  
  
Farshadfar E, Mahmodi N, Yaghotipoor A (2011). “AMMI stability value and
simultaneous estimation of yield and yield stability in bread wheat
(*Triticum aestivum* L.).” *Australian Journal of Crop Science*,
**5**(13), 1837–1844.  
  
Rao AR, Prabhakaran VT (2005). “Use of AMMI in simultaneous selection of
genotypes for yield and stability.” *Journal of the Indian Society of
Agricultural Statistics*, **59**, 76–82.

## See also

[`AMGE.AMMI`](https://ajaygpb.github.io/ammistability/reference/AMGE.AMMI.md),
[`ASI.AMMI`](https://ajaygpb.github.io/ammistability/reference/ASI.AMMI.md),
[`ASTAB.AMMI`](https://ajaygpb.github.io/ammistability/reference/ASTAB.AMMI.md),
[`AVAMGE.AMMI`](https://ajaygpb.github.io/ammistability/reference/AVAMGE.AMMI.md),
[`DA.AMMI`](https://ajaygpb.github.io/ammistability/reference/DA.AMMI.md),
[`DZ.AMMI`](https://ajaygpb.github.io/ammistability/reference/DZ.AMMI.md),
[`EV.AMMI`](https://ajaygpb.github.io/ammistability/reference/EV.AMMI.md),
[`FA.AMMI`](https://ajaygpb.github.io/ammistability/reference/FA.AMMI.md),
[`MASV.AMMI`](https://ajaygpb.github.io/ammistability/reference/MASV.AMMI.md),
[`SIPC.AMMI`](https://ajaygpb.github.io/ammistability/reference/SIPC.AMMI.md),
[`ZA.AMMI`](https://ajaygpb.github.io/ammistability/reference/ZA.AMMI.md)

## Examples

``` r
library(agricolae)
data(plrv)
model <- with(plrv, AMMI(Locality, Genotype, Rep, Yield, console=FALSE))

yield <- aggregate(model$means$Yield, by= list(model$means$GEN),
               FUN=mean, na.rm=TRUE)[,2]
stab <- DZ.AMMI(model)$DZ
genotypes <- rownames(DZ.AMMI(model))

# With default ssi.method (farshadfar)
SSI(y = yield, sp = stab, gen = genotypes)
#>                 SP SSI rSP rY    means
#> 102.18  0.26393535  37  14 23 26.31947
#> 104.22  0.22971564  21   8 13 31.28887
#> 121.31  0.32031744  34  19 15 30.10174
#> 141.28  0.39838535  23  22  1 39.75624
#> 157.26  0.53822924  33  28  5 36.95181
#> 163.9   0.26659011  42  15 27 21.41747
#> 221.19  0.19563325  29   3 26 22.98480
#> 233.11  0.25167755  27  10 17 28.66655
#> 235.6   0.46581370  28  24  4 38.63477
#> 241.2   0.21481887  28   6 22 26.34039
#> 255.7   0.30862904  31  17 14 30.58975
#> 314.12  0.22603261  25   7 18 28.17335
#> 317.6   0.20224771  14   5  9 35.32583
#> 319.20  0.50675112  29  26  3 38.75767
#> 320.16  0.23280596  30   9 21 26.34808
#> 342.15  0.25989774  36  12 24 26.01336
#> 346.2   0.37125512  45  20 25 23.84175
#> 351.26  0.43805896  31  23  8 36.11581
#> 364.21  0.07409309  12   2 10 34.05974
#> 402.7   0.02004533  20   1 19 27.47748
#> 405.2   0.26238837  29  13 16 28.98663
#> 406.12  0.28179394  28  16 12 32.68323
#> 427.7   0.20176581  11   4  7 36.19020
#> 450.3   0.25465368  17  11  6 36.19602
#> 506.2   0.30899851  29  18 11 33.26623
#> Canchan 0.37201039  41  21 20 27.00126
#> Desiree 0.52005815  55  27 28 16.15569
#> Unica   0.48083049  27  25  2 39.10400

# With  ssi.method = "rao"
SSI(y = yield, sp = stab, gen = genotypes, method = "rao")
#>                 SP        SSI rSP rY    means
#> 102.18  0.26393535  1.5536988  14 23 26.31947
#> 104.22  0.22971564  1.8193399   8 13 31.28887
#> 121.31  0.32031744  1.5545939  19 15 30.10174
#> 141.28  0.39838535  1.7570779  22  1 39.75624
#> 157.26  0.53822924  1.5459114  28  5 36.95181
#> 163.9   0.26659011  1.3869397  15 27 21.41747
#> 221.19  0.19563325  1.6878048   3 26 22.98480
#> 233.11  0.25167755  1.6641025  10 17 28.66655
#> 235.6   0.46581370  1.6538090  24  4 38.63477
#> 241.2   0.21481887  1.7134093   6 22 26.34039
#> 255.7   0.30862904  1.5922105  17 14 30.58975
#> 314.12  0.22603261  1.7307783   7 18 28.17335
#> 317.6   0.20224771  2.0595024   5  9 35.32583
#> 319.20  0.50675112  1.6259792  26  3 38.75767
#> 320.16  0.23280596  1.6476346   9 21 26.34808
#> 342.15  0.25989774  1.5545233  12 24 26.01336
#> 346.2   0.37125512  1.2718506  20 25 23.84175
#> 351.26  0.43805896  1.5966462  23  8 36.11581
#> 364.21  0.07409309  3.5881882   2 10 34.05974
#> 402.7   0.02004533 10.0539968   1 19 27.47748
#> 405.2   0.26238837  1.6447637  13 16 28.98663
#> 406.12  0.28179394  1.7171135  16 12 32.68323
#> 427.7   0.20176581  2.0898536   4  7 36.19020
#> 450.3   0.25465368  1.9010808  11  6 36.19602
#> 506.2   0.30899851  1.6787677  18 11 33.26623
#> Canchan 0.37201039  1.3738642  21 20 27.00126
#> Desiree 0.52005815  0.8797586  27 28 16.15569
#> Unica   0.48083049  1.6568004  25  2 39.10400

# Changing the ratio of weights for Rao's SSI
SSI(y = yield, sp = stab, gen = genotypes, method = "rao", a = 0.43)
#>                 SP       SSI rSP rY    means
#> 102.18  0.26393535 1.1572429  14 23 26.31947
#> 104.22  0.22971564 1.3638258   8 13 31.28887
#> 121.31  0.32031744 1.2279220  19 15 30.10174
#> 141.28  0.39838535 1.4944208  22  1 39.75624
#> 157.26  0.53822924 1.3514985  28  5 36.95181
#> 163.9   0.26659011 0.9944318  15 27 21.41747
#> 221.19  0.19563325 1.1529329   3 26 22.98480
#> 233.11  0.25167755 1.2483375  10 17 28.66655
#> 235.6   0.46581370 1.4291726  24  4 38.63477
#> 241.2   0.21481887 1.2263072   6 22 26.34039
#> 255.7   0.30862904 1.2531668  17 14 30.58975
#> 314.12  0.22603261 1.2678419   7 18 28.17335
#> 317.6   0.20224771 1.5421234   5  9 35.32583
#> 319.20  0.50675112 1.4194898  26  3 38.75767
#> 320.16  0.23280596 1.1981670   9 21 26.34808
#> 342.15  0.25989774 1.1519083  12 24 26.01336
#> 346.2   0.37125512 0.9899993  20 25 23.84175
#> 351.26  0.43805896 1.3577771  23  8 36.11581
#> 364.21  0.07409309 2.1759278   2 10 34.05974
#> 402.7   0.02004533 4.8338929   1 19 27.47748
#> 405.2   0.26238837 1.2459704  13 16 28.98663
#> 406.12  0.28179394 1.3457828  16 12 32.68323
#> 427.7   0.20176581 1.5712389   4  7 36.19020
#> 450.3   0.25465368 1.4901748  11  6 36.19602
#> 506.2   0.30899851 1.3401295  18 11 33.26623
#> Canchan 0.37201039 1.0925852  21 20 27.00126
#> Desiree 0.52005815 0.6785528  27 28 16.15569
#> Unica   0.48083049 1.4391795  25  2 39.10400
```
