# Metapopulation dynamics in *Tetranychus urticae*

(doi link to Zenodo archive redacted for anonymised peer-review)

This repo contains data and code needed to re-do the analyses and figures in our manuscript:

"Individual heterogeneity: a neglected driver of metapopulation dynamics"
(authors names redacted for anonymised peer-review)



data in `.csv` format are in the `data` folder, R script in`Rmd` format (including detailed information about the analysis) in the `R` folder.

This folder is a RStudio project folder, and the script uses the `here` package (see also [here](https://github.com/jennybc/here_here)). This means all file paths are relative, and the analysis should work on your computer no questions asked, whether you use the R project or not, no code line to change as long as you download the entire repo (you just need to install all the needed packages first, of course).

If you run the script for the first time, models and some other time-consuming outputs will be saved in the `R_output` folder so you don't have to re-run them everytime.
