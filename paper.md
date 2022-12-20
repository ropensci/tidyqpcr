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
 - name: Institute for Cell Biology and SynthSys, School of Biological Sciences, The University of Edinburgh,
   index: 1
date: 30 May 2022
bibliography: paper.bib
---

# Summary

Quantitative polymerase chain reaction (qPCR) is a fundamental technique in molecular biology to detect and quantify DNA and RNA.
Here we present the tidyqpcr software package for user-friendly qPCR analysis using the tidyverse suite of R packages. 
tidyqpcr offers a consistent user interface and structure for qPCR analysis, within the tidyverse paradigm of spreadsheet-like rectangular data frames and generic functions that build up complex analyses in a series of simple steps.
tidyqpcr focuses on experimental design in microwell plates, and relative quantification using changes in quantification cycle ($\Delta Cq$).
Overall, tidyqpcr empowers scientists to conduct reproducible, flexible, and best-practice compliant quantitative PCR analysis. 


# Statement of need

Quantitative PCR is among the most common techniques in biological and biomedical research, used for the quantification of DNA and RNA.
There is a critical need for rigorous analysis and reporting of qPCR experiments, codified  in the [minimum information for publication of quantitative real-time PCR experiments (MIQE) guidelines](https://academic.oup.com/clinchem/article/55/4/611/5631762) [@Bustin:2009].
Yet it is common for qPCR to be analysed either by closed-source software supplied by the manufacturers of PCR machines, or by highly variable, in-house analysis scripts that have not been peer-reviewed.

Our package, tidyqpcr, addresses the need for a qPCR analysis package that fully integrates with the user-friendly tidyverse, encourages the use of MIQE best-practice compliant experimental design, and provides detailed example analysis pipelines as R vignettes.
Following the tidy data paradigm integrates tidyqpcr into the wider collection of data analysis packages provided by the tidyverse, while being accessible to novice users of R.
Other open-source libraries for qPCR analysis are available with distinct aims. 
HTqPCR [@Dvinge:2009], ReadqPCR/NormqPCR [@Perkins:2012], and qpcR [@Spiess:2018] have similar and in some respects greater functionality than tidyqpcr, but follow object oriented approaches with specialised data objects.
By contrast, pcr [@Ahmed:2018] aligns most closely with tidyqpcr but its function inputs are not tidy data frames.
Although alternatives are available, tidyqpcr's aims and approach are distinct: to improve the quality of qPCR experiments from plate design to analysis, by exposing all data in a consistent tidy format that integrates with the tidyverse.

tidyqpcr aims to be:

* Empowering: tidyqpcr combines a free, open-source qPCR analysis R package with online teaching materials. 
* Reproducible: tidyqpcr scripts produce paper-ready figures straight from raw data with identical results across computers.
* Flexible: tidyqpcr follows the 'tidy' data paradigm to ensure scalability and adaptability.
* Best-practice compliant: tidyqpcr encourages standardised, reliable experimental design by prioritising MIQE-compliant best practices.

tidyqpcr can be used to analyse qPCR data from any nucleic acid source - DNA for qPCR or ChIP-qPCR, RNA for RT-qPCR.
Currently tidyqpcr has functions that explicitly support relative quantification by the $\Delta Cq$ method, but not yet absolute quantification.

tidyqpcr's current features allow users to:

* use a single data type for analysis as every object is a tibble / data frame.
* lay out and display 96/384-well plates for easy experimental setup (`label_plate_rowcol`, `create_blank_plate`, ...).
* consistently describe samples and target amplicons with reserved variable names (`sample_id`, `target_id`).
* flexibly assign metadata to samples for visualisation with [ggplot2](https://ggplot2.tidyverse.org/) (see vignettes).
* read in quantification cycle (Cq) and raw data from Roche LightCycler machines with single-channel fluorescence (`read_lightcycler_1colour_cq`, `read_lightcycler_1colour_raw`).
* calibrate primer sets including estimating efficiencies and visualization of curves (`calculate_efficiency`).
* visualize amplification and melt curves (`calculate_drdt_plate`)
* perform normalisation and relative quantification to one or more reference targets by the $\Delta Cq$ method (`calculate_normcq`, `calculate_deltacq_bysampleid`).
* estimate differential expression across multiple samples by the $\Delta \Delta Cq$ method (`calculate_deltadeltacq_bytargetid`).
* accelerate further downstream analysis and visualization by writing tidy data frames that are fully compatible with the tidyverse suite.

We have conducted a series of user interviews to improve tidyqpcr's capabilities and documentation.
The ease-of-use and documentation of tidyqpcr will enable efficient best-practice analysis of qPCR data by both novice and experienced programmers.


# Acknowledgements

We thank everyone in the eLife Innovation Leaders 2020 program for all their help developing tidyqpcr, in particular program leader Emmy Tsang and our mentor Aidan Budd.
We thank Sander Granneman, Stefanie Butland and Sean Hughes for feedback and encouragement.
We thank rOpenSci, Julia Gustavsen, and Kelsey Montgomery for constructive reviews.
We thank all those who have participated in interviews, including; Flic Anderson, Jamie Auxillos, David Barrass, Rosey Bayne, Elliott Chapman, Magnus Gwynne, Liz Hughes, Chris Katanski and Stuart McKellar.
Edward Wallace is a Sir Henry Dale Fellow, jointly funded by Wellcome and the Royal Society [208779/Z/17/Z].
Samuel Haynes is funded by the EASTBIO UKRI-BBSRC DTP [BB/M010996/1].


# References
