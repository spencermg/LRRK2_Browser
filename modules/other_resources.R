#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

otherResourcesUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("otherResources"))
}
# =========================================================================
# SERVER FUNCTION
# =========================================================================

otherResourcesServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        output$otherResources <- renderUI({
            fluidRow(
                column(
                    width = 3,
                    align = "center",
                    tags$a(
                        href = "https://www.ncbi.nlm.nih.gov/gene/120892",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "NCBI (Gene Summary)"
                    ),
                ),
                column(
                    width = 3,
                    align = "center",
                    tags$a(
                        href = "https://www.uniprot.org/uniprotkb/Q5S007/entry",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "UniProt (Protein Function & Domains)"
                    ),
                ),
                column(
                    width = 3,
                    align = "center",
                    tags$a(
                        href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=LRRK2",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "GeneCards (Integrative database)"
                    ),
                ),
                column(
                    width = 3,
                    align = "center",
                    tags$a(
                        href = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000188906",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "Ensembl (Functional Annotations)"
                    ),
                )
            )
        })
    })
}
