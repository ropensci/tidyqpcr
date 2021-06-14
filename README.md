[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

# tidyqpcr - Quantitative PCR analysis in the tidyverse.

*Empowering scientists to conduct reproducible, flexible, and MIQE best-practice compliant quantitative PCR analysis.*

# Contents

* [Motivation](#Motivation)
	* [Empowering](#Empowering)
	* [Reproducible](#Reproducible)
	* [Flexible](#Flexible)
	* [Best-practice compliant](#Best-practice-compliant)
* [Status](#Status)
* [Features](#Features)
	* [Current features](#Current-features)
	* [Future priorities](#Future-priorities)
* [Getting started](#Getting-started)
	* [Installing tidyqpcr](#Installing-tidyqpcr)
	* [Using tidyqpcr](#Using-tidyqpcr)
* [Contribute](#Contribute)

# Motivation

### Empowering

tidyqpcr combines a free, open-source qPCR analysis R package with online teaching materials. 

We want our users to be able to know and understand what happens at every step of their analysis. Users are able to know what occurs at each step as all tidyqpcr tools are open source and follow the FAIR principles - [Findable, Accessible, Interoperable, and Reusable](https://www.force11.org/group/fairgroup/fairprinciples). Users should also find each step understandable as we aim to produce educational resources as extensions of [data carpentry](https://datacarpentry.org/) workshops, such as [Data Analysis and Visualization in R for Ecologists](https://datacarpentry.org/R-ecology-lesson/), accessible to beginner programmers. 

### Reproducible

tidyqpcr scripts produce paper-ready figures straight from raw data with identical results across computers.

We want to promote reproducible research so collaborators, reviewers or students can easily confirm and extend results and conclusions. tidyqpcr analysis will repeat exactly on different computers, enabling scientists to share raw data and analysis scripts rather than just  processed figures. An R or R markdown script using tidyqpcr to analyse a set of qPCR data could be directly uploaded to a repository such as [figshare](https://figshare.com/), as encouraged by many journal publishers.

### Flexible

tidyqpcr follows the 'tidy' data paradigm to ensure scalability and adaptability.

We want to create a tool that is flexible enough to analyse high or low throughput experimental data whilst integrating easily into other data analyses. tidyqpcr uses powerful generic data science tools from [the tidyverse](https://www.tidyverse.org/) R package, lightly overlaid with qPCR-specific scripts. As far as possible, every object in tidyqpcr is stored as a generic tibble / data frame. Manipulating and plotting qPCR data without creating bespoke data structures allows tidyqpcr scripts to be easily integrated and scaled according to the needs of your experiments.

### Best-practice compliant

tidyqpcr encourages standardised, reliable experimental design by prioritising MIQE-compliant best practices.

We want to make it easier for scientists to produce reliable and interpretable results. The final version of tidyqpcr will, by default, request the relevant experimental conditions and assay characteristics, as described in the [MIQE guidelines](https://academic.oup.com/clinchem/article/55/4/611/5631762), to allow reviewers/readers to rigorously assess the validity of a result. See "Future Priorities" below to get updates on tidyqpcr's MIQE compliant features.

# Status

As of August 2020, this software is in development. [Edward Wallace](https://github.com/ewallace) wrote basic functions and documentation needed to do qPCR analysis in [the Wallace lab](https://ewallace.github.io/), and is making them freely available. [Sam Haynes](https://github.com/dimmestp) is helping develop as part of the [eLife Open Innovation Leaders programme](https://elifesciences.org/labs/fdcb6588/innovation-leaders-2020-introducing-the-cohort). 

## News

* August 2020, relative quantification (delta delta Cq) added with function `calculate_deltadeltacq_bytargetid`, and a vignette illustrationg this with small data from a 96-well plate.
* June 2020, upgrades that break previous code. All function and variable names have been changed to snake case, i.e. lower case with underscore. Commits up to #ee6d192 change variable and function names. tidyqpcr now uses `sample_id` for nucleic acid sample (replaces Sample or SampleID), `target_id` for primer set/ probe (replaces TargetID or Probe), `prep_type` for nucleic acid preparation type (replaces Type), and `cq` for quantification cycle (replaces Cq or Ct). 
It should be possible to upgrade old analysis code by (case-sensitive) search and replace. 

Alternatively, pre-April 2020 analysis code should run from release v0.1-alpha, see [releases](https://github.com/ewallace/tidyqpcr/releases).

# Features 

tidyqpcr can be used to analyse qPCR data from any nucleic acid source - DNA for qPCR or ChIP-qPCR, RNA for RT-qPCR.

Currently tidyqpcr has functions that support relative quantification, but not yet absolute quantification.


## Current features

* every object is a tibble / data frame, no special data classes to learn
* lay out and display 96/384-well plates for easy experimental setup (`label_plate_rowcol`, `create_blank_plate`, ...)
* flexible assignment of metadata to samples for visualisation with [ggplot2](https://ggplot2.tidyverse.org/) (see vignettes)
* read-in Cq and raw data from Roche LightCycler machines with single-channel fluorescence (`read_lightcycler_1colour_cq`, `read_lightcycler_1colour_raw`)
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
devtools::install_github("ewallace/tidyqpcr",build_vignettes = TRUE) ## Vignettes require cowplot package

## Alternatively, install without building the vignetttes to remove cowplot dependency (Not recommended as vignettes contain the tutorials on using tidyqpcr)
devtools::install_github("ewallace/tidyqpcr")
```
**Note**
older versions of the remote package automatically convert warnings to errors during installation. Please update your remotes package to >2.3.0 in order to remove this default.

Then load tidyqpcr as a standard package:

```
library(tidyqpcr)
```
**Note**
tidyqpcr automatically imports and loads several external packages for basic functionality, including; tidy, dplyr and ggplot2. This allows tidyqpcr to be used immediately but may cause NAMESPACE clashes if the user already has many other package libraries loaded. Restarting the R session and loading tidyqpcr separately may solve such issues.

## Using tidyqpcr

The best place to start is the vignettes, which offer tutorials and example data analyses including figures. Currently there are 3 vignettes:

* [IntroDesignPlatesetup](vignettes/platesetup_vignette.Rmd) - Introduction to designing an experiment and setting up a plate plan in tidyqpcr.
* [DeltaCq96wellExample](vignettes/deltacq_96well_vignette.Rmd) - Example analysis of 96-well RT-qPCR data including relative quantification with delta Cq, from a real experiment.
* [MultifactorialExample](vignettes/multifactor_vignette.Rmd) - Example design and analysis of a (real) multifactorial RT-qPCR experiment.
* [PrimerCalibration](vignettes/calibration_vignette.Rmd) - Example design and analysis of calibrating qPCR primer sets from a (real) experimental test

To find these from your R session, enter `browseVignettes(package="tidyqpcr")`. 


Individual R functions are also documented, use R's standard help system after loading the package, e.g. `?create_blank_plate`. To see a list of all the functions and links to their help pages use `help(package="tidyqpcr")`.

# Contribute

We would be delighted to work with you to answer questions, add features, and fix problems. Please [file an issue](https://github.com/ewallace/tidyqpcr/issues) or email Edward dot Wallace at his University email address, (ed.ac.uk).

## Code of conduct

We will be following the [code of conduct from the tidyverse](https://dplyr.tidyverse.org/CODE_OF_CONDUCT).

## How to contribute code: style, checking, development cycle

If you want to fix bugs or add features yourself, that's great. tidyqpcr development aims to follow best practices:

* first file an issue on the tidyqpcr github repository, then create a branch (in a fork) of the repository that references the issue number, e.g. `add-exciting-feature-55` if it were issue 55.
* add code in the branch, continue to ask questions on the issue ticket.
* follow the [tidyverse style guide](https://style.tidyverse.org/).
* document functions with [roxygen2](https://roxygen2.r-lib.org/), as described in [the R packages book](http://r-pkgs.had.co.nz/man.html).
* check the package with `R CMD check` / `devtools::check()`, as explained in [the R packages book](http://r-pkgs.had.co.nz/check.html).
* including, check that all the vignettes run.
* put in a pull request to the main repository, we will review, then we will accept or suggest changes.
