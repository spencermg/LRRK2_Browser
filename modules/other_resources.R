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

height <- 75
horizontal_spacing <- 20
vertical_spacing <- 30

otherResourcesServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        output$otherResources <- renderUI({
            tagList(
                tags$div(
                    style = paste0(
                        "display: flex; flex-wrap: wrap; ",
                        "column-gap: ", horizontal_spacing, "px; ",
                        "row-gap: ", vertical_spacing, "px; ",
                        "justify-content: center; margin-bottom: 20px;"
                    ),
                    tags$a(
                        href = "https://pdgenetics.shinyapps.io/GP2Browser/",
                        target = "_blank",
                        title = "GP2 GWAS Browser -- characterization of underlying molecular mechanisms of LRRK2 that contribute to PD risk and pathogenesis",
                        tags$img(
                            src = "gp2_gwas_browser.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://gp2.broadinstitute.org/gene/ENSG00000188906",
                        target = "_blank",
                        title = "GP2 Genome Browser -- variant-level data across LRRK2 from individuals of diverse ancestries",
                        tags$img(
                            src = "gp2_genome_browser.png",
                            style = paste0("height:", height, "px;")
                        )
                    )
                ),
                tags$div(
                    style = paste0(
                        "display: flex; flex-wrap: wrap; ",
                        "column-gap: ", horizontal_spacing, "px; ",
                        "row-gap: ", vertical_spacing, "px; ",
                        "justify-content: center;"
                    ),
                    tags$a(
                        href = "https://www.ncbi.nlm.nih.gov/clinvar/?term=%22LRRK2%22%5BGENE%5D&redir=gene",
                        target = "_blank",
                        title = "ClinVar -- relationships between LRRK2 variations and their associated health outcomes",
                        tags$img(
                            src = "clinvar.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000188906",
                        target = "_blank",
                        title = "Ensembl -- gene annotations, variation, regulation, and comparative genomics for LRRK2 across a wide range of species",
                        tags$img(
                            src = "ensembl.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=LRRK2",
                        target = "_blank",
                        title = "GeneCards -- integrated genomic, proteomic, transcriptomic, ad functional information for LRRK2",
                        tags$img(
                            src = "genecards.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://gnomad.broadinstitute.org/gene/ENSG00000188906",
                        target = "_blank",
                        title = "Genome Aggregation Database (gnomAD) -- harmonized exome and genome data with summary statistics and population-level frequencies for LRRK2 variants",
                        tags$img(
                            src = "gnomad.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://ldlink.nih.gov/ldpair",
                        target = "_blank",
                        title = "LDLink -- suite of web applications for exploring linkage disequilibrium for LRRK2 variants",
                        tags$img(
                            src = "ldlink.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://www.mdsgene.org/genes/PARK-LRRK2",
                        target = "_blank",
                        title = "Movement Disorder Society Genetic Mutation Database (MDSGene) -- comprehensive overview of data on movement disorder patients with causative LRRK2 variants",
                        tags$img(
                            src = "mdsgene.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://www.ncbi.nlm.nih.gov/gene/120892",
                        target = "_blank",
                        title = "National Center for Biotechnology Information (NCBI) -- comprehensive biomedical information about LRRK2",
                        tags$img(
                            src = "ncbi.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://www.omim.org/entry/609007",
                        target = "_blank",
                        title = "Online Mendelian Inheritance in Man (OMIM) -- database for LRRK2-phenotype relationships",
                        tags$img(
                            src = "omim.png",
                            style = paste0("height:", height, "px;")
                        )
                    ),
                    tags$a(
                        href = "https://www.uniprot.org/uniprotkb/Q5S007/entry",
                        target = "_blank",
                        title = "UniProt -- protein sequence and functional information for LRRK2 variants",
                        tags$img(
                            src = "uniprot.png",
                            style = paste0("height:", height, "px;")
                        )
                    )
                )
            )
        })
    })
}
