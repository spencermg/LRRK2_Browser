#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("shiny", quietly = TRUE)) {
    install.packages("shiny")
}
if (!requireNamespace("shinydashboardPlus", quietly = TRUE)) {
    install.packages("shinydashboardPlus")
}
if (!requireNamespace("DT", quietly = TRUE)) {
    install.packages("DT")
}
if (!requireNamespace("plotly", quietly = TRUE)) {
    install.packages("plotly")
}
if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
    install.packages("shinyWidgets")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}

library(bslib)
library(data.table) 
library(DT) 
library(plotly)
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)

runApp("/Users/grantsm/Desktop/projects/lrrk2_browser")
