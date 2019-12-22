# tidyqpcr - Quantitative PCR analysis with the tidyverse.

tidyqpcr is an R package for quantitative PCR analysis that aims to:
* be FAIR - [Findable, Accessible, Interoperable, and Reusable](https://www.force11.org/group/fairgroup/fairprinciples)
* use modern data science tools from [the tidyverse](https://www.tidyverse.org/)
* be accessible to novice programmers, e.g. those who have taken [data carpentry in R](https://datacarpentry.org/R-ecology-lesson/) 
* be compliant with [MIQE guidelines](1373/clinchem.2008.112797)

Current capabilities include:

* lay out and display 96/384-well plates for easy experimental setup
* read-in Ct and raw data from Roche LightCycler machines with single-channel fluorescence
* calibration of primer sets including estimating efficiencies and visualization of curves
* visualization of amplification and melt curves
* normalization of Ct data to one or more reference probe sets by delta delta count method
* flexible assignment of metadata to samples for visualisation with [ggplot2](https://ggplot2.tidyverse.org/)

Future priorities include:

* including primer efficiencies in quantification
* an open-source and tested Ct/Cq calculation algorithm
* extend to 1536-well plates 
* files for automatic plate loading with [Opentrons](https://opentrons.com/) and [Labcyte Echo](https://www.labcyte.com/products/liquid-handling/echo-liquid-handlers) liquid handlers.

As of Dec 2019, this software is in development. Edward Wallace wrote basic functions and documentation needed to do qPCR analysis in [the Wallace lab](https://ewallace.github.io/), and is making them freely available. We would be delighted to work with you to answer questions, add features, and fix problems. Please file an issue or email Edward dot Wallace at his University email address, (ed.ac.uk). 

## Install

```
library(devtools)
devtools::install_github("ewallace/tidyqpcr",build_vignettes = TRUE)
```

## Use
The best place to start is the vignettes, see `vignette(package="tidyqpcr")` from your R session.

* calibration_vignette.Rmd - qPCR primer calibration
* multifactor_vignette.Rmd - RT-qPCR of gene expression in multifactorial design

Individual R functions are also documented, use R's standard help system after loading the package.
