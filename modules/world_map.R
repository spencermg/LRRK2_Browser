#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

worldMapUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            width = 12,
            tags$p(
                paste(
                    "World map with proportions of LRRK2 disease-related variants per ancestry illustrated",
                    "in donut charts. Numbers in brackets indicate the total number of individuals carrying",
                    "LRRK2 disease-related variants per ancestry, and colors reflect protein domains in which",
                    "variants are located, according to the protein diagram above."
                ),
                style = "font-size: 14px; color: #555; margin-top: 8px; text-align: center;"
            )
        ),
        column(
            width = 12,
            tags$p(
                paste0(
                    "*Note: Charts are positioned approximately according to the geographic origin of the ",
                    "respective ancestries and cohorts and are intended for illustrative purposes only; their ",
                    "placement should not be interpreted as precise geographic localization."
                ),
                style = "font-size: 11px; font-style: italic; color: #555; margin-top: 8px; text-align: center;"
            )
        ),
        uiOutput(ns("worldMap"))
    )
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

worldMapServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        output$worldMap <- renderUI({
            tags$div(
                tags$img(
                    src = "world_map.png",
                    style = "width: 75%; height: auto; display: block; margin: 0 auto;"
                )
            )
        })
    })
}
