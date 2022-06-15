<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test coverage](https://codecov.io/gh/ropensci/tidyqpcr/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/tidyqpcr/branch/main)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/5287/badge)](https://bestpractices.coreinfrastructure.org/projects/5287)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tidyqpcr)](https://cran.r-project.org/package=tidyqpcr)
[![R-CMD-check](https://github.com/ropensci/tidyqpcr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/tidyqpcr/actions)
<!-- badges: end -->

# tidyqpcr - Quantitative PCR analysis in the tidyverse.

*Empowering scientists to conduct reproducible, flexible, and MIQE best-practice compliant quantitative PCR analysis.*

# Contents

* [Motivation](#Motivation)
	* [Empowering](#Empowering)
	* [Reproducible](#Reproducible)
	* [Flexible](#Flexible)
	* [Best-practice compliant](#Best-practice-compliant)
* [Getting started](#Getting-started)
	* [Installing tidyqpcr](#Installing-tidyqpcr)
	* [Using tidyqpcr](#Using-tidyqpcr)
* [Status](#Status)
    * [News - see NEWS.md](NEWS.md)
* [Features](#Features)
	* [Current features](#Current-features)
	* [Future priorities](#Future-priorities)
* [Contribute](#Contribute)


# Motivation

Quantitative Polymerase Chain Reaction (qPCR) is a highly adaptable experimental technique used across biology and medicine to measure the amounts of nucleic acids (DNA or RNA).
tidyqpcr is a software package for qPCR data analysis that builds on the tidyverse collection of data science tools in the R programming language.

### Empowering

tidyqpcr combines a free, open-source qPCR analysis R package with online teaching materials. 

We want our users to be able to know and understand what happens at every step of their analysis.
Users are able to know what occurs at each step as all tidyqpcr tools are open source and follow the FAIR principles - [Findable, Accessible, Interoperable, and Reusable](https://force11.org/info/the-fair-data-principles/). 
Users should also find each step understandable as we aim to produce educational resources as extensions of [data carpentry](https://datacarpentry.org/) workshops, such as [Data Analysis and Visualization in R for Ecologists](https://datacarpentry.org/R-ecology-lesson/), accessible to beginner programmers. 

### Reproducible

tidyqpcr scripts produce paper-ready figures straight from raw data with identical results across computers.

We want to promote reproducible research so collaborators, reviewers or students can easily confirm and extend results and conclusions. tidyqpcr analysis will repeat exactly on different computers, enabling scientists to share raw data and analysis scripts rather than just  processed figures. An R or R markdown script using tidyqpcr to analyse a set of qPCR data could be directly uploaded to a repository such as [figshare](https://figshare.com/), as encouraged by many journal publishers.

### Flexible

tidyqpcr follows the 'tidy' data paradigm to ensure scalability and adaptability.

We want to create a tool that is flexible enough to analyse high or low throughput experimental data whilst integrating easily into other data analyses. tidyqpcr uses powerful generic data science tools from [the tidyverse](https://www.tidyverse.org/) R package, lightly overlaid with qPCR-specific scripts. As far as possible, every object in tidyqpcr is stored as a generic tibble / data frame. Manipulating and plotting qPCR data without creating bespoke data structures allows tidyqpcr scripts to be easily integrated and scaled according to the needs of your experiments.

### Best-practice compliant

tidyqpcr encourages standardised, reliable experimental design by following the Minimum Information for Publication of Quantitative Real-Time PCR Experiments (MIQE) best practices.

We want to make it easier for scientists to produce reliable and interpretable results. The MIQE best practices are a framework to facilitate the full disclosure of all reagents, sequences, and analysis methods necessary to enable other investigators to reproduce results. The final version of tidyqpcr will, by default, request the relevant experimental conditions and assay characteristics, as described in the [MIQE guidelines](https://academic.oup.com/clinchem/article/55/4/611/5631762), to allow reviewers/readers to rigorously assess the validity of a result. See "Future Priorities" below to get updates on tidyqpcr's MIQE compliant features.


# Getting started

## Installing tidyqpcr

First install [R](https://www.r-project.org/). 

#### For Windows users

Next, you need a working installation of [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html).

Jeffrey Leek made [slides on installation and testing of Rtools](http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/).

### Install via devtools (for all R users)

Install the devtools R package, see [devtools installation instructions](https://www.r-project.org/nosvn/pandoc/devtools.html). 

```
library(devtools)
devtools::install_github("ropensci/tidyqpcr", build_vignettes = TRUE)

## Alternatively, install without building the vignettes 
## (Not recommended as vignettes contain the tutorials on using tidyqpcr)
devtools::install_github("ropensci/tidyqpcr")
```

**Note**
older versions of the remotes package automatically convert warnings to errors during installation.
Please update your remotes package to >2.3.0 in order to remove this default.

Then load tidyqpcr as a standard package:

```
library(tidyqpcr)
```

**Note**
tidyqpcr automatically imports and loads several external packages for basic functionality, including; tidy, dplyr and ggplot2.
This allows tidyqpcr to be used immediately but may cause NAMESPACE clashes if the user already has many other package libraries loaded.
Restarting the R session and loading tidyqpcr separately may solve such issues.

## Using tidyqpcr

The best place to start is by viewing the articles on the [tidyqpcr website](https://docs.ropensci.org/tidyqpcr/index.html).
Here you will find the vignettes, which offer tutorials and example data analyses including figures.
Currently there are 4 vignettes:

* [IntroDesignPlatesetup](https://docs.ropensci.org/tidyqpcr/articles/platesetup_vignette.html) - Introduction to designing an experiment and setting up a plate plan in tidyqpcr.
* [DeltaCq96wellExample](https://docs.ropensci.org/tidyqpcr/articles/deltacq_96well_vignette.html) - Example analysis of 96-well RT-qPCR data including relative quantification with delta Cq, from a real experiment.
* [MultifactorialExample](https://docs.ropensci.org/tidyqpcr/articles/multifactor_vignette.html) - Example design and analysis of a (real) multifactorial RT-qPCR experiment.
* [PrimerCalibration](https://docs.ropensci.org/tidyqpcr/articles/calibration_vignette.html) - Example design and analysis of calibrating qPCR primer sets from a (real) experimental test

To find these from your R session, enter `browseVignettes(package="tidyqpcr")`. 

Individual R functions are also documented, use R's standard help system after loading the package, e.g. `?create_blank_plate`. To see a list of all the functions and links to their help pages use `help(package="tidyqpcr")`.

A basic use case for designing a 12 well plate is given below, see [IntroDesignPlatesetup](https://docs.ropensci.org/tidyqpcr/articles/platesetup_vignette.html) for more details.

```
rowkey4 <- tibble(
  well_row = LETTERS[1:4],
  target_id = c("ACT1", "BFG2", "CDC19", "DED1")
)


colkey3 <- tibble(
  well_col = 1:3,
  sample_id = c("rep1", "rep2", "rep3"),
  prep_type = "+RT"
)

create_blank_plate(well_row = LETTERS[1:4], well_col = 1:3)

plate_plan12 <- label_plate_rowcol(
  plate = create_blank_plate(well_row = LETTERS[1:4], well_col = 1:3),
  rowkey = rowkey4,
  colkey = colkey3
)

display_plate_qpcr(plate_plan12)
```

# Status

As of June 2022, this software is fully useable, and under active development development.
It is particularly good at designing qPCR experiments in microwell plates (96-well and 384-well), and at relative quantification by the delta Cq method.

[Edward Wallace](https://github.com/ewallace) wrote basic functions and documentation needed to do qPCR analysis in [the Wallace lab](https://ewallace.github.io/), then started building them into an R package.
[Sam Haynes](https://github.com/dimmestp) is actively developing, initially as part of the [eLife Open Innovation Leaders programme 2020](https://elifesciences.org/labs/fdcb6588/innovation-leaders-2020-introducing-the-cohort).

If there is a feature that you need for your work, please ask us! 

## News - see [NEWS.md](NEWS.md).

# Features 

tidyqpcr can be used to analyse qPCR data from any nucleic acid source - DNA for qPCR or ChIP-qPCR, RNA for RT-qPCR.

Currently tidyqpcr has functions that support relative quantification by the delta Cq method, but not yet absolute quantification.


## Current features

* every object is a tibble / data frame, no special data classes to learn
* lay out and display 96/384-well plates for easy experimental setup (`label_plate_rowcol`, `create_blank_plate`, ...).
* consistently describe samples and target amplicons with reserved variable names (`sample_id`, `target_id`).
* flexibly assign metadata to samples for visualisation with [ggplot2](https://ggplot2.tidyverse.org/) (see vignettes).
* read in quantification cycle (Cq) and raw data from Roche LightCycler machines with single-channel fluorescence (`read_lightcycler_1colour_cq`, `read_lightcycler_1colour_raw`).
* calibration of primer sets including estimating efficiencies and visualization of curves (`calculate_efficiency`, and see vignettes)
* visualization of amplification and melt curves (`calculate_drdt_plate`, and see vignettes)
* delta Cq: normalization/ relative quantification of Cq data to one or more reference targets by delta count method (`calculate_normcq`, `calculate_deltacq_bysampleid`)
* delta delta Cq: normalization of delta Cq data across multiple samples (`calculate_deltadeltacq_bytargetid`)

## Future priorities

* including primer efficiencies in quantification
* an open-source and tested Cq calculation function, from amplification curves
* multi-colour (hydrolysis probe) detection
* extend to 1536-well plates 
* metadata handling compatible with RDML format
* files for automatic plate loading with [Opentrons](https://opentrons.com/) and [Labcyte Echo](https://www.labcyte.com/products/liquid-handling/echo-liquid-handlers) liquid handlers.

# Comparison of qPCR R packages with respect to the MIQE guidelines

Table of package features corresponding to the essential information on qPCR validation and data analysis, that are outlined by the MIQE guidelines for publication of qPCR results. 

| MIQE Guidelines | tidyqpcr | HTqPCR | NormqPCR | qpcR | pcr |
--- | --- | --- | --- | --- | ---
Version | 0.5.0 | 1.48.0 | 1.40.0 | 1.4.1  | 1.2.2
For SYBR Green I, Cq of the NTC | Yes + Docs | Yes + Docs | Yes + Docs | Yes | Yes
Calibration curves with slope and y intercept | Slope | No | No | Yes + Docs | Yes + Docs
PCR efficiency calculated from slope | Yes + Doc | No | Yes + Doc | Yes + Doc | Yes
r2 of calibration curve | Yes + Doc | No | No | Yes | Yes
Linear dynamic range ‡ | No | No | No | No | No
Cq variation at LOD ‡ | No | No | No | No | No
Evidence for LOD ‡ | No | No | No | No | No
If multiplex, efficiency and LOD of each assay ‡ | No | No | No | No | No
Method of Cq determination | N/A | N/A | Sigmoidal model selection | Sigmoidal model selection | N/A
Outlier identification and disposition | No | Yes | Yes | Yes | No
Results for NTCs | Yes + Doc | Yes + Docs | Yes + Docs | Yes | Yes
Justification of number and choice of reference genes | User defined (vignettes encourage 3) | User defined (vignettes encourage 2) | Automatic Selection (vignettes encourage 2) | User defined | One 
Description of normalization method | Relative | Relative | Relative | Relative or absolute | Relative
Number and stage (reverse transcription or qPCR) of technical replicates | User defined (vignettes encourage 3) | User defined (vignettes encourage 3) | User defined (vignettes encourage 2) | User defined | User defined (vignettes encourage 6)
Repeatability | No† | Yes + Docs | Yes + Docs | Yes + Docs | No
Statistical methods for results significance | No† | Yes + Docs | No | Yes + Docs | Yes

Note: 
- Yes means the package includes the functionality to complete this analysis.
- Yes + Docs means this step is explicitly shown in either the function documentation or a vignette.
- No† means that the package lacks explicit functionality, but generic R capabilities for statistical testing can be applied to the data.
- ‡ Linear dynamic range and limit of detection (LOD) calculations would be enabled by these packages from additional short scripted analyses from a well-designed experiment, but the functionality is not specifically documented.


# Contribute

We would be delighted to work with you to answer questions, add features, and fix problems. Please [file an issue](https://github.com/ropensci/tidyqpcr/issues) or email Edward dot Wallace at his University email address, (ed.ac.uk).

## Code of conduct

This package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.

## How to contribute code: style, checking, development cycle

If you want to fix bugs or add features yourself, that's great. tidyqpcr development aims to follow best practices which we have outlined in the [CONTRIBUTING.md](.github/CONTRIBUTING.md) file in the .github folder.

## Thank you

Many thanks to everyone who has helped with tidyqpcr

Users and interviewees: Jamie Auxillos, Rosey Bayne, Liz Hughes, Rachael Murray, Elliott Chapman, Laura Tuck, Amy Newell, David Barrass, Christopher Katanski,  Magnus Gwynne and Stuart McKeller.
Reviewers: @seaaan, @kelshmo and @jooolia
