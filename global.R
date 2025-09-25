#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  shiny,
  shinydashboard,
  shinydashboardPlus,
  DT,
  plotly,
  shinyWidgets,
  data.table,
  bslib
)

source("modules/gene_overview.R", local = TRUE)
source("modules/annotation_summary_table.R", local = TRUE)
source("modules/protein_diagram.R", local = TRUE)
source("modules/cdna_diagram.R", local = TRUE)
source("modules/gene_var_table.R", local = TRUE)
