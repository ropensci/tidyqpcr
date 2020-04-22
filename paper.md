---
title: 'tidyqpcr: Quantitative PCR analysis in the tidyverse.'
tags:
  - quantitative PCR
  - qPCR
  - tidyverse
  - R
  - MIQE
authors:
  - name: Edward W. J. Wallace
    orcid: 0000-0001-8025-6361
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Samuel J. Haynes
    affiliation: 1
affiliations:
 - name: Institute for Cell Biology, School of Biological Sciences, The University of Edinburgh,
   index: 1
date: 22 April 2020
bibliography: paper.bib # does not yet exist

---

# Summary

Quantitative PCR is a fundamental technique in molecular biology to detect and quantify DNA and RNA. Although technologies for qPCR measurement are mature, the software landscape for qPCR analysis is still in flux. Here we present the tidyqpcr software package for user-friendly qPCR analysis using the tidyverse suite of R packages. tidyqpcr empowers scientists to conduct reproducible, flexible, and MIQE best-practice compliant quantitative PCR analysis.

## Background

Quantitative PCR.
Key steps in data analysis.
MIQE guidelines.

Existing software landscape.
Manufacturer software is closed-source and platform-dependent.
Open-source packages from Spiess lab, e.g. qpcR.

Impact of the tidyverse.
Principles of consistent data formats.

This created an opportunity for an open-source qPCR analysis package that exploits the strengths of the tidyverse.


## Motivation for tidyqpcr

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

# Key functions of tidyqpcr

# Acknowledgements

We thank Emmy Tsang and the eLife Innovation Leaders 2020 program for all their help developing this, in particular our mentor Aidan Budd.
Edward Wallace is a Sir Henry Dale Fellows, jointly funded by the Wellcome Trust and the Royal Society (Grant Number 208779/Z/17/Z).
Samuel Haynes is funded by the EASTBIO UKRI-BBSRC DTP.

# References