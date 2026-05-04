#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

library(bslib)
library(data.table)
library(DT)
library(plotly)
library(shiny)
library(shinyalert)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(shinyWidgets)

# Source local files
source("util.R")
mod_files <- list.files("modules", pattern = "\\.[Rr]$", full.names = TRUE)
invisible(
    lapply(mod_files, function(f) {
        sys.source(f, envir = globalenv())
    })
)
