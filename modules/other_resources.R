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

height <- 90
horizontal_spacing <- 30

otherResourcesServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        output$otherResources <- renderUI({
            tags$div(
                style = paste0("display: flex; flex-wrap: wrap; gap: ", horizontal_spacing, "px; justify-content: center;"),
                tags$a(
                    href = "https://www.ncbi.nlm.nih.gov/clinvar/?term=%22LRRK2%22%5BGENE%5D&redir=gene",
                    target = "_blank",
                    tags$img(
                        src = "icons/clinvar.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000188906",
                    target = "_blank",
                    tags$img(
                        src = "icons/ensembl.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=LRRK2",
                    target = "_blank",
                    tags$img(
                        src = "icons/genecards.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://gnomad.broadinstitute.org/gene/ENSG00000188906",
                    target = "_blank",
                    tags$img(
                        src = "icons/gnomad.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://pdgenetics.shinyapps.io/GP2Browser/",
                    target = "_blank",
                    tags$img(
                        src = "icons/gwas_browser.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://gp2.broadinstitute.org/gene/ENSG00000188906",
                    target = "_blank",
                    tags$img(
                        src = "icons/gwas_browser.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://ldlink.nih.gov/ldpair",
                    target = "_blank",
                    tags$img(
                        src = "icons/ldlink.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://www.mdsgene.org/genes/PARK-LRRK2",
                    target = "_blank",
                    tags$img(
                        src = "icons/mdsgene.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://www.ncbi.nlm.nih.gov/gene/120892",
                    target = "_blank",
                    tags$img(
                        src = "icons/ncbi.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://www.omim.org/entry/609007",
                    target = "_blank",
                    tags$img(
                        src = "icons/omim.png",
                        style = paste0("height:", height, "px;")
                    )
                ),
                tags$a(
                    href = "https://www.uniprot.org/uniprotkb/Q5S007/entry",
                    target = "_blank",
                    tags$img(
                        src = "icons/uniprot.png",
                        style = paste0("height:", height, "px;")
                    )
                )
            )
        })
    })
}
