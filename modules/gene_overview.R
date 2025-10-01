#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

geneOverviewUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("geneOverview"))
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

geneOverviewServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        output$geneOverview <- renderUI({
            fluidRow(
                # Chromosome number
                column(
                    width = 3,
                    align = "center",
                    descriptionBlock(
                        header = tags$b("CHR"),
                        text = 12,
                        rightBorder = TRUE,
                        marginBottom = FALSE
                    )
                ),
                # GrCh38 coordinates
                column(
                    width = 3,
                    align = "center",
                    descriptionBlock(
                        header = tags$b("BP (GrCh38)"),
                        text = "40224997 - 40369285",
                        rightBorder = TRUE,
                        marginBottom = FALSE
                    )
                ),
                # T2T coordinates
                column(
                    width = 3,
                    align = "center",
                    descriptionBlock(
                        header = tags$b("BP (T2T)"),
                        text = "40177355 - 40321422",
                        rightBorder = TRUE,
                        marginBottom = FALSE
                    )
                ),
                # GrCh37 coordinates
                column(
                    width = 3,
                    align = "center",
                    descriptionBlock(
                        header = tags$b("BP (GrCh37)"),
                        text = "40618799 - 40763087",
                        rightBorder = TRUE,
                        marginBottom = FALSE
                    )
                )
            )
        })
    })
}