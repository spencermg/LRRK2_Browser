#!/usr/bin/env Rscript

# =========================================================================
# DEFINE VARIABLES
# =========================================================================
header_tooltips <- c(
    "Variant (GrCh38)" = "Genomic coordinates for variants (chromosome, position, reference allele, alternate allele)",
    "rsID" = "",
    "cDNA change" = "Based on reference transcript ENST00000298910.12",
    "Exon #" = "Based on reference transcript ENST00000298910.12",
    "AA change" = "",
    "Protein domain" = "",
    "Region" = "",
    "Functional consequence" = "",
    "CADD" = "",
    "Clinvar Pathogenic" = "",
    "Conservation score" = "",
    "Kinase activity (mean pRAB10/RAB10)" = "",
    "PD frequency (WGS)" = "",
    "Control frequency (WGS)" = "",
    "PD frequency (Imputed)" = "",
    "Control frequency (Imputed)" = "",
    "PD frequency (Raw genotyping)" = "",
    "Control frequency (Raw genotyping)" = "",
    "PD frequency (Clinical exome)" = "",
    "gnomAD allele frequency" = "Using whole-genome data from the gnomAD v4 dataset (76,215 genomes across all ancestries)"
)


# =========================================================================
# UI FUNCTION
# =========================================================================

geneVarTableUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 12,
                tags$p(
                    "*Some ancestries do not have reported allele frequencies through gnomAD, and are thus left empty.",
                    style = "font-size: 11px; font-style: italic; color: #555; margin-top: 8px; text-align: center;"
                ),
                tags$p(
                    "**Individual ancestry assignments are not present for clinical exome data, and are thus excluded from ancestry-specific tables.",
                    style = "font-size: 11px; font-style: italic; color: #555; margin-top: 8px; text-align: center;"
                )
            ),
            # Add dropdown for ancestry selection
            column(
                width = 3,
                selectInput(
                    inputId  = ns("dataset"),
                    label    = "Choose ancestry:",
                    choices  = NULL,
                    width    = "100%"
                )
            ),
            # Add buttons to vilter variants
            column(
                width = 6,
                shinyWidgets::checkboxGroupButtons(
                    inputId      = ns("filters"),
                    label        = "Filters:",
                    choiceNames  = list(
                        tags$span(title = "Variants in exonic regions", "Exonic"),
                        tags$span(title = "GP2-determined high-confidence variants (according to at least two of Clinvar, HGMD, and MDSGene)", "Pathogenic"),
                        tags$span(title = "Predicted deleterious variants (CADD > 20)", "Deleterious"),
                        tags$span(title = "Evolutionarily conserved variants (conservation score > 5)", "Conserved"),
                        tags$span(title = "Kinase-active variants (mean pRAB10/RAB10 > 1.40)", "Kinase active")
                    ),
                    choiceValues = c("Exonic", "Pathogenic", "Deleterious", "Conserved", "Kinase active"),
                    selected     = NULL,
                    status       = "primary",
                    justified    = FALSE,
                    width        = "100%" ,
                    checkIcon    = list(
                        yes = icon("ok", lib = "glyphicon"),
                        no  = icon("remove", lib = "glyphicon")
                    )
                )
            ),
            column(
                width = 3,
                textInput(
                    inputId = ns("search"),
                    label   = "Search:",
                    placeholder = "e.g., 'G2019S' or 'Pathogenic'"
                )
            )
        ),
        
        fluidRow(
            column(
                width = 12,
                DT::DTOutput(ns("table"))
            )
        )
    )
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

geneVarTableServer <- function(id, all_tables_cleaned, clicked_variant = NULL) {
    moduleServer(id, function(input, output, session) {
        # Make cols_to_keep reactive based on ancestry selection
        cols_to_keep <- reactive({
            base_cols <- c(
                "Variant (GrCh38)",
                "rsID",
                "cDNA change",
                "Exon #",
                "AA change",
                "Protein domain",
                "Region",
                "Functional consequence",
                "CADD",
                "Clinvar Pathogenic",
                "Conservation score",
                "Kinase activity (mean pRAB10/RAB10)",
                "PD frequency (WGS)",
                "Control frequency (WGS)",
                "PD frequency (Imputed)",
                "Control frequency (Imputed)",
                "PD frequency (Raw genotyping)",
                "Control frequency (Raw genotyping)"
            )
            
            # Only include exome column for Combined ancestry
            if (input$dataset == "Combined") {
                base_cols <- c(base_cols, "PD frequency (Clinical exome)")
            }
            
            # Add remaining columns
            base_cols <- c(
                base_cols,
                "gnomAD allele frequency"
            )
            
            base_cols
        })

        ns <- session$ns

        # Initialize dataset selector
        updateSelectInput(
            session, "dataset",
            choices  = names(all_tables_cleaned),
            selected = if ("Combined" %in% names(all_tables_cleaned)) "Combined" else names(all_tables_cleaned)[1]
        )

        # Apply selected filters
        dat_rx <- reactive({
            dat <- all_tables_cleaned[[ req(input$dataset) ]]

            freq_cols <- c(
                "PD frequency (WGS)",
                "Control frequency (WGS)",
                "PD frequency (Imputed)",
                "Control frequency (Imputed)",
                "PD frequency (Raw genotyping)",
                "Control frequency (Raw genotyping)"
            )
            
            # Only include exome in frequency check for Combined
            if (input$dataset == "Combined") {
                freq_cols <- c(freq_cols, "PD frequency (Clinical exome)")
            }
            
            # Filter to only available columns
            freq_cols <- freq_cols[freq_cols %in% names(dat)]
            
            keep <- apply(dat[, ..freq_cols], 1, function(row) {
                any(!is.na(row) & row != 0)
            })
            dat <- dat[keep, ]

            sel <- if (is.null(input$filters)) character(0) else input$filters

            if ("Exonic" %in% sel) {
                region <- dat[["Region"]]
                keep <- !is.na(region) & region == "exonic"
                dat <- dat[ keep, , drop = FALSE]
            }
            if ("Pathogenic" %in% sel) {
                pathogenic <- dat[["Pathogenic"]]
                keep <- !is.na(pathogenic) & pathogenic == 1
                dat <- dat[ keep, , drop = FALSE]
            }
            if ("Deleterious" %in% sel) {
                deleterious <- dat[["Deleterious?"]]
                keep <- !is.na(deleterious) & deleterious == "Yes"
                dat <- dat[ keep, , drop = FALSE]
            }
            if ("Conserved" %in% sel) {
                conserved <- dat[["Conserved?"]]
                keep <- !is.na(conserved) & conserved == "Yes"
                dat <- dat[ keep, , drop = FALSE]
            }
            if ("Kinase active" %in% sel) {
                kinase_active <- dat[["Kinase active?"]]
                keep <- !is.na(kinase_active) & kinase_active == "Yes"
                dat <- dat[ keep, , drop = FALSE]
            }
            search_term <- trimws(input$search)
            if (nzchar(search_term)) {
                matches <- apply(dat, 1, function(row) {
                    any(grepl(search_term, row, ignore.case = TRUE, fixed = FALSE))
                })
                dat <- dat[matches, ]
            }

            # Only keep columns that exist in the data
            cols <- cols_to_keep()
            cols <- cols[cols %in% names(dat)]
            dat <- dat[, ..cols]
            dat
        })

        output$table <- DT::renderDT({
            dat <- dat_rx()

            # Add tooltips to column headers
            make_tooltip_headers <- function(cols, tooltips) {
                sapply(cols, function(col) {
                    if (!is.null(tooltips[[col]])) {
                        sprintf(
                            '<span title="%s">%s</span>',
                            htmltools::htmlEscape(tooltips[[col]]),
                            htmltools::htmlEscape(col)
                        )
                    } else {
                        htmltools::htmlEscape(col)
                    }
                }, USE.NAMES = FALSE)
            }
            headers <- make_tooltip_headers(colnames(dat), header_tooltips)

            # Make first column (variant ID) clickable
            if (nrow(dat) > 0) {
                dat[[1]] <- paste0(
                    '<a href="#" class="variant-link" style="color: #0C8DC3; text-decoration: underline;">',
                    dat[[1]],
                    '</a>'
                )
            }

            # Columns to show in scientific notation
            sci_cols <- c(
                "PD frequency (WGS)", 
                "Control frequency (WGS)", 
                "PD frequency (Imputed)", 
                "Control frequency (Imputed)", 
                "PD frequency (Raw genotyping)",
                "Control frequency (Raw genotyping)",
                "PD frequency (Clinical exome)",
                "gnomAD allele frequency"
            )
            sci_cols <- intersect(sci_cols, colnames(dat))
            sci_targets <- match(sci_cols, colnames(dat)) - 1L

            DT::datatable(
                dat,
                colnames = headers,
                extensions = "Buttons",
                rownames   = FALSE,
                escape     = FALSE,
                selection  = "none",
                options    = list(
                    dom          = "Blrtip",
                    buttons      = c("copy", "csv", "excel", "pdf", "print"),
                    paging       = TRUE,
                    pageLength   = 10,
                    lengthChange = TRUE,
                    lengthMenu   = list(c(10,25,50,100,500,-1), c("10","25","50","100","500","All")),
                    scrollX      = TRUE,
                    deferRender  = TRUE,
                    columnDefs   = if (length(sci_targets)) list(
                        list(
                            targets   = sci_targets,
                            className = "dt-right",
                            render    = DT::JS(
                                "function(data, type, row, meta){",
                                "  if (data === null || data === undefined || data === '') {",
                                "       return (type === 'display') ? '' : null;",
                                "  }",
                                "  var num = Number(data);",
                                "  if (!isFinite(num)) {",
                                "       return (type === 'display') ? data : null;",
                                "  }",
                                "  if (num === 0) {",
                                "       return (type === 'display') ? '0' : 0;",
                                "  }",
                                "  if (type === 'display') {",
                                "       if (Math.abs(num) < 0.001) {",
                                "            return num.toExponential(3);",
                                "       } else {",
                                "            return Math.round(num * 10000) / 10000;",
                                "       }",
                                "  }",
                                "  return num;",
                                "}"
                            )
                        )
                    ) else NULL
                )
            ) |>
            DT::formatStyle(
                columns         = colnames(dat),
                backgroundColor = "#FFFFFF",
                color           = "black"
            )
        })

        # Capture clicks on variant links
        observeEvent(input$table_cell_clicked, {
            info <- input$table_cell_clicked
            dat <- dat_rx()

            # Only trigger on first column (variant ID)
            if (!is.null(info$col) && info$col == 0 && !is.null(clicked_variant)) {
                variant_id <- dat[[1]][info$row]
                ancestry <- input$dataset
                clicked_variant(list(
                    variant_id = variant_id,
                    ancestry = ancestry,
                    counter = Sys.time()
                ))
            }
        })
    })
}
