#!/usr/bin/env Rscript

shiny::runApp(
    "/Users/grantsm/Desktop/projects/lrrk2_browser",
    host = "0.0.0.0", 
    port = 8080, 
    launch.browser = TRUE
)
