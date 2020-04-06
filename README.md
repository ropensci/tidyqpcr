# tidyqpcr - Quantitative PCR analysis in the tidyverse.

*Empowering scientists to conduct reproducible, flexible, and MIQE best-practice compliant quantitative PCR analysis.*

### Empowering

tidyqpcr combines a free, open-source qPCR analysis R package with online teaching materials. 

We want our users to be able to know and understand what happens at every step of their analysis. Users are able to know what occurs at each step as all tidyqpcr tools are open source and follow the FAIR principles - [Findable, Accessible, Interoperable, and Reusable](https://www.force11.org/group/fairgroup/fairprinciples). Users should also find each step understandable as we aim to produce educational resources as extensions of [data carpentry](https://datacarpentry.org/) workshops, such as [Data Analysis and Visualization in R for Ecologists](https://datacarpentry.org/R-ecology-lesson/), accessible to beginner programmers. 

### Reproducible

tidyqpcr scripts produce paper-ready figures straight from raw data with identical results across computers.

We want to promote reproducible research so collaborators, reviewers or students can easily confirm and extend results and conclusions. tidyqpcr analysis will repeat exactly on different computers, enabling scientists to share raw data and analysis scripts rather than just  processed figures. An R or R markdown script using tidyqpcr to analyse a set of qPCR data could be directly uploaded to a repository such as [figshare](https://figshare.com/), as encouraged by many journal publishers.

### Flexible

tidyqpcr follows the 'tidy' data paradigm to ensure scalability and adaptability.

We want to create a tool that is flexible enough to analyse high or low throughput experimental data whilst integrating easily into multi-omic data analyses. tidyqpcr uses powerful generic data science tools from [the tidyverse](https://www.tidyverse.org/) R package, lightly overlaid with qPCR-specific scripts. Manipulating and plotting qPCR data without creating bespoke data structures allows tidyqpcr scripts to be easily integrated and scaled according to the needs of your experiments.

### Best-practice compliant

tidyqpcr encourages standardised, reliable experimental design by prioritising MIQE-compliant best practices.

We want to make it easier for scientists to produce reliable and interpretable results. The final version of tidyqpcr will, by default, request the relevant experimental conditions and assay characteristics, as described in the [MIQE guidelines](https://academic.oup.com/clinchem/article/55/4/611/5631762), to allow reviewers/readers to rigorously assess the validity of a result. See "Future Priorities" below to get updates on tidyqpcr's MIQE compliant features.

## Status

As of Feb 2020, this software is in development. Edward Wallace wrote basic functions and documentation needed to do qPCR analysis in [the Wallace lab](https://ewallace.github.io/), and is making them freely available. Sam Haynes is helping develop as part of the [eLife Open Innovation Leaders programme](https://elifesciences.org/labs/fdcb6588/innovation-leaders-2020-introducing-the-cohort). 

## Contribute

We would be delighted to work with you to answer questions, add features, and fix problems. Please [file an issue](https://github.com/ewallace/tidyqpcr/issues) or email Edward dot Wallace at his University email address, (ed.ac.uk).

We will be following the [code of conduct from the tidyverse](https://dplyr.tidyverse.org/CODE_OF_CONDUCT).

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

## Installing tidyqpcr

First install [R](https://www.r-project.org/). 

#### For Windows users

Next, you need a working installation of [Rtools]().

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

## Using tidyqpcr

The best place to start is the vignettes, which offer tutorials and example data analyses including figures. Currently there are 3 vignettes:

* [IntroDesignPlatesetup](vignettes/platesetup_vignette.Rmd) - Introduction to designing an experiment and setting up a plate plan in tidyqpcr.
* [MultifactorialExample](vignettes/multifactor_vignette.Rmd) - Example design and analysis of a (real) multifactorial qPCR experiment.
* [ProbeCalibration](vignettes/platesetup_vignette.Rmd) - Example design and analysis of calibrating qPCR primer sets from a (real) experimental test

To find these from your R session, enter `browseVignettes(package="tidyqpcr")`. 

Individual R functions are also documented, use R's standard help system after loading the package, e.g. `?create_blank_plate`. To see a list of all the functions and links to their help pages use `help(package="tidyqpcr")`.
