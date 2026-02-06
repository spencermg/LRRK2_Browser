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
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://www.ncbi.nlm.nih.gov/clinvar/?term=%22LRRK2%22%5BGENE%5D&redir=gene",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "Clinvar"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000188906",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "Ensembl"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=LRRK2",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "GeneCards"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://gnomad.broadinstitute.org/gene/ENSG00000188906",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "gnomAD"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://pdgenetics.shinyapps.io/GP2Browser/",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "GP2 GWAS Browser"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://gp2.broadinstitute.org/gene/ENSG00000188906",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "GP2 Variant Browser"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://ldlink.nih.gov/ldpair",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "LD Pair"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://www.mdsgene.org/genes/PARK-LRRK2",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "MDSGene"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://www.ncbi.nlm.nih.gov/gene/120892",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "NCBI"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://www.omim.org/entry/609007",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "OMIM"
                    ),
                ),
                column(
                    width = 2,
                    align = "center",
                    tags$a(
                        href = "https://www.uniprot.org/uniprotkb/Q5S007/entry",
                        target = "_blank",
                        style = "font-size:18px; font-weight:bold; color:#0C8DC3; text-decoration:none;",
                        "UniProt"
                    ),
                ),
            )
        })
    })
}
