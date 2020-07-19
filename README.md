# Metapopulation dynamics in *Tetranychus urticae*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3951429.svg)](https://doi.org/10.5281/zenodo.3951429)

This repo contains data and code needed to re-do the analyses and figures in our manuscript:

"Individual heterogeneity: a neglected driver of metapopulation dynamics"
(by Stefano Masier, Maxime Dahirel, Frederik Mortier, Dries Bonte)



data in `.csv` format are in the `data` folder, R script in`Rmd` and `html` format (including detailed information about the analysis) in the `R` folder.

This folder is a RStudio project folder, and the script uses the `here` package (see also here). This means all file paths are relative, and the analysis should work on your computer no questions asked, whether you use the R project or not, no code line to change as long as you download the entire repo (you just need to install all the needed packages first, of course).

If you run the script for the first time, models and some other time-consuming outputs will be saved in the `R_output` folder so you don't have to re-run them everytime.
