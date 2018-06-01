#!/usr/bin/env Rscript
# command line script to call Histogram_model.Rmd
# first argument should be the full path to the histo.tsv file
# ! Update path the Histogram_model.Rmd if running on a different computer.

library(tidyverse)
library(stringr)

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop(paste("Supply the full path to the histo file.\n",
             "Usage: Rscript --vanilla run_histogram_model.R <full input file path>"),
       call. = FALSE)
}

histo_file = args[1]

out_name = histo_file %>% str_replace("tsv", "fits")

rmarkdown::render("/scratch/users/mfincker/cell_cycle_model/Histogram_model.Rmd",
                  output_file = out_file,
                  params = list(histo_file = histo_file))