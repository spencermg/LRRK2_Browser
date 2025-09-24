#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Required packages
pkgs <- c(
  "shiny",
  "shinydashboard",
  "shinydashboardPlus",
  "DT",
  "plotly",
  "shinyWidgets",
  "data.table",
  "bslib"
)

# Load packages, installing any that are missing
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) {
  install.packages(to_install)
}
invisible(lapply(pkgs, library, character.only = TRUE))
