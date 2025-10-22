#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load dependencies
if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
}
pacman::p_load(
    shiny,
    shinydashboard,
    shinydashboardPlus,
    shinyalert,
    DT,
    plotly,
    shinyWidgets,
    data.table,
    bslib
)

# Source local files
source ("util.R")
mod_files <- list.files("modules", pattern = "\\.[Rr]$", full.names = TRUE)
invisible(
    lapply(mod_files, function(f) {
        sys.source(f, envir = globalenv())
    })
)
