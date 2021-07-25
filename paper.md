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
    orcid: 0000-0002-3366-1812
    affiliation: 1
affiliations:
 - name: Institute for Cell Biology, School of Biological Sciences, The University of Edinburgh,
   index: 1
date: 25 July 2021
bibliography: paper.bib # does not yet exist

---

# Summary

Quantitative polymerase chain reaction (qPCR) is a fundamental technique in molecular biology to detect and quantify DNA and RNA.
Here we present the tidyqpcr software package for user-friendly qPCR analysis using the tidyverse suite of R packages. 
tidyqpcr offers a consistent user interface and structure for qPCR analysis, within the tidyverse paradigm of spreadsheet-like rectangular data frames and generic functions that build up complex analyses in a series of simple steps.
tidyqpcr focuses on experimental design in microwell plates, and relative quantification using changes in quantification cycle ($\Delta Cq$).
Overall, tidyqpcr empowers scientists to conduct reproducible, flexible, and best-practice compliant quantitative PCR analysis. 


# Statement of need

Quantitative PCR is among the most common techniques in biological and biomedical research, used for quantification of DNA and RNA.
There is a critical need for rigorous analysis and reporting of qPCR experiments, codified  in the [minimum information on a qPCR experiment (MIQE) guidelines](https://academic.oup.com/clinchem/article/55/4/611/5631762).
Yet it is still common for qPCR to be analyzed by closed-source software supplied with machines, or unreliable home-written scripts.
Some open-source libraries for qpcr analysis are available, notably qpcR (cite) which offers many features within an object-oriented approach using S4 classes.
Recently, the tidyverse suite of generic data-science tools democratised a user-friednly paradigm of tidy data (spreadsheet-like rectangular data frames) and generic functions that build up complex analyses in a series of simple steps.
This created a need for qPCR analysis package that integrates easily with the tidyverse, including data visualization with ggplot2.

Our package, tidyqpcr, addresses the need for a user-friendly qPCR analysis in the tidyverse paradigm. 
tidyqpcr aims to be:

* Empowering: tidyqpcr combines a free, open-source qPCR analysis R package with online teaching materials. 
* Reproducible: tidyqpcr scripts produce paper-ready figures straight from raw data with identical results across computers.
* Flexible: tidyqpcr follows the 'tidy' data paradigm to ensure scalability and adaptability.
* Best-practice compliant: tidyqpcr encourages standardised, reliable experimental design by prioritising MIQE-compliant best practices.

tidyqpcr can be used to analyse qPCR data from any nucleic acid source - DNA for qPCR or ChIP-qPCR, RNA for RT-qPCR.
Currently tidyqpcr has functions that support relative quantification, but not yet absolute quantification.

tidyqpcr's current features allow users to:

* use a single data type for analysis as every object is a tibble / data frame.
* lay out and display 96/384-well plates for easy experimental setup (`label_plate_rowcol`, `create_blank_plate`, ...).
* flexibly assign metadata to samples for visualisation with [ggplot2](https://ggplot2.tidyverse.org/) (see vignettes).
* read in quantification cycle (Cq) and raw data from Roche LightCycler machines with single-channel fluorescence (`read_lightcycler_1colour_cq`, `read_lightcycler_1colour_raw`).
* calibrate primer sets including estimating efficiencies and visualization of curves (`calculate_efficiency`).
* visualize of amplification and melt curves (`calculate_drdt_plate`)
* perform normalization and relative quantification to one or more reference targets by the $\Delta Cq$ method (`calculate_normcq`, `calculate_deltacq_bysampleid`).
* delta delta Cq: normalization of delta Cq data across multiple samples (`calculate_deltadeltacq_bytargetid`).
* accelerate downstream analysis and visualization by writing tidy data frames summarized target id (amplicon or gene), sample id, or by any other user-supplied metadata variable.

We have conducted a series of user interviews to improve tidyqpcr's capabilities and documentation.
The ease-of-use and documentation of tidyqpcr will enable efficient best-practice analysis of qPCR data by both novice and experienced programmers.


# Acknowledgements

We thank everyone in the eLife Innovation Leaders 2020 program for all their help developing tidyqpcr, in particular program leader Ammy Tsang and our mentor Aidan Budd.
We thank Stefanie Butland and Sean Hughes for feedback and encouragement.
We thank all those who have agreed to user interviews, including Flic Anderson, Jamie Auxillos, David Barrass, Rosey Bayne, Elliott Chapman, Magnus Gwynne, Liz Hughes, Chris Katanski, .
Edward Wallace is a Sir Henry Dale Fellow, jointly funded by the Wellcome Trust and the Royal Society (Grant Number 208779/Z/17/Z).
Samuel Haynes is funded by the EASTBIO UKRI-BBSRC DTP.


# References