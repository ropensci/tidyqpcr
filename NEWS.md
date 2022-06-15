## News

* June 2022, removed plot helper functions `scale_..._nice` and `scale_loglog` from tidyqpcr, because those capabilities are now available in the [scales package](https://scales.r-lib.org) using `label_log` and similar functions. Older code may need to change `scale_y_log10nice` to `scale_y_log10(labels = scales::label_log())`, for example.
* May 2022, Improvements in documentation and testing. Reorganized `display_plate` function to be more flexible, so older code will need to use `display_plate_qpcr` to ensure that `sample_id` and `target_id` info displays. Updated to v0.5.
* January 2022, Improvements in documentation and argument-checking for v0.4.
* October 2021, Unit tests now cover over 75% of tidyqpcr code.
* June 2021, [tidyqpcr blogpost in eLife labs](https://elifesciences.org/labs/f23e268f/tidyqpcr-quantitative-pcr-analysis-in-the-tidyverse)
* August 2020, relative quantification (delta delta Cq) added with function `calculate_deltadeltacq_bytargetid`, and a vignette illustrating this with a small data set from a 96-well plate.
* June 2020, upgrades that break previous code. All function and variable names have been changed to snake case, i.e. lower case with underscore. Commits up to #ee6d192 change variable and function names. tidyqpcr now uses `sample_id` for nucleic acid sample (replaces Sample or SampleID), `target_id` for primer set/ probe (replaces TargetID or Probe), `prep_type` for nucleic acid preparation type (replaces Type), and `cq` for quantification cycle (replaces Cq or Ct). It should be possible to upgrade old analysis code by (case-sensitive) search and replace. 
