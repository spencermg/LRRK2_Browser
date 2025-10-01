#!/usr/bin/env Rscript

shiny::runApp(
    "./",
    host = "0.0.0.0", 
    port = 8080, 
    launch.browser = TRUE
)
