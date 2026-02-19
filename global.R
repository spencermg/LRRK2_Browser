#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load dependencies
if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
}
pacman::p_load(
    bslib,
    data.table,
    DT,
    plotly,
    shiny,
    shinyalert,
    shinydashboard,
    shinydashboardPlus,
    shinyjs,
    shinymanager,
    shinyWidgets
)

# Source local files
source ("util.R")
mod_files <- list.files("modules", pattern = "\\.[Rr]$", full.names = TRUE)
invisible(
    lapply(mod_files, function(f) {
        sys.source(f, envir = globalenv())
    })
)
