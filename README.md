This is an R package used for jump detection (or changepoint detection) for high dimensional non-stationary time series. The corresponding paper this method is based on can be found at https://arxiv.org/abs/2410.23706. Within the documentation of the package, any references to sections/equations in the paper reflect those of the version v2 submitted to arXiv in November 2024. To install the package from Github, run:

```R
library(devtools)
devtools::install_github("daveveitch/ajdn")
```

A useful walkthrough of using the ajdn package with real and simulated data is provided via a vignette. To install the package with the vignette, and then view the vignette, use:

```R
library(devtools)
devtools::install_github("daveveitch/ajdn",build_vignettes=TRUE)
vignette('ajdn_walkthrough')
```
