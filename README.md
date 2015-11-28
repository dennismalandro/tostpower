<!-- README.md is generated from README.Rmd. Please edit that file -->
The tostpower package is a power/sample-size calculator for experiments designed to detect either that two population means are *equivalent*, or that one population mean is *not inferior to* second a second population mean, for either paired or unpaired observations.

The power calculation uses on the excellent `pmvt` function in the [mvtnorm](http://mvtnorm.R-forge.R-project.org) package of Genz, et al.

It is the computational engine underlying the shiny app of the same name, but it can be used as a standalone package.

``` r
install_github(username = 'dennismalandro', 'tostpower')

library(tostpower)
tost_power(20)
```
