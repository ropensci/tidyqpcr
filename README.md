# tidyqpcr - tidy Quantitative PCR analysis.

Empowering scientists to conduct reproducible, flexible, and MIQE best-practice compliant quantitative PCR analysis.

### Empowering

* free software
* open-source software so users know *exactly* what happened during every step of the analysis 
* accessible to beginner programmers, e.g. those who have taken [data carpentry in R](https://datacarpentry.org/R-ecology-lesson/) 
* FAIR - [Findable, Accessible, Interoperable, and Reusable](https://www.force11.org/group/fairgroup/fairprinciples)

### Reproducible

* enable completely scripted qPCR analysis, from instrument data to figures
* repeat the same analysis on a different computer
* easily extend to new datasets

### Flexible

* allow user flexibility to extend analyses with generic tools
* use modern data science tools from [the tidyverse](https://www.tidyverse.org/)

### Best-practice compliant

* aims to be compliant with [MIQE guidelines](1373/clinchem.2008.112797)
* aims to make it easier to design experiments that are MIQE-compliant, e.g. with multiple reference genes
* we still have work to do to add MIQE-recommended quantification features, see "Future priorities"

## Status

As of Feb 2020, this software is in development. Edward Wallace wrote basic functions and documentation needed to do qPCR analysis in [the Wallace lab](https://ewallace.github.io/), and is making them freely available. We would be delighted to work with you to answer questions, add features, and fix problems. Please file an issue or email Edward dot Wallace at his University email address, (ed.ac.uk). 


# Features 

## Current features include:

* lay out and display 96/384-well plates for easy experimental setup
* read-in Ct and raw data from Roche LightCycler machines with single-channel fluorescence
* calibration of primer sets including estimating efficiencies and visualization of curves
* visualization of amplification and melt curves
* normalization of Ct data to one or more reference probe sets by delta delta count method
* flexible assignment of metadata to samples for visualisation with [ggplot2](https://ggplot2.tidyverse.org/)

## Future priorities include:

* including primer efficiencies in quantification
* an open-source and tested Ct/Cq calculation algorithm
* extend to 1536-well plates 
* files for automatic plate loading with [Opentrons](https://opentrons.com/) and [Labcyte Echo](https://www.labcyte.com/products/liquid-handling/echo-liquid-handlers) liquid handlers.


# Getting started

## Install instructions

#### For Windows users

Before anything else, that you  need a working installation of [Rtools]().

Jeffrey Leek made [slides on installation and testing of Rtools](http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/).

### For all R users

Install the devtools R package, see [devtools installation instructions](https://www.r-project.org/nosvn/pandoc/devtools.html). 

```
library(devtools)
devtools::install_github("ewallace/tidyqpcr",build_vignettes = TRUE)
```

Then load tidyqpcr as a standard package:

```
library(tidyqpcr)
```

## Use

The best place to start is the vignettes, see `vignette(package="tidyqpcr")` from your R session.

* calibration_vignette.Rmd - qPCR primer calibration
* multifactor_vignette.Rmd - RT-qPCR of gene expression in multifactorial design

Individual R functions are also documented, use R's standard help system after loading the package, e.g. `?create_blank_plate`.
