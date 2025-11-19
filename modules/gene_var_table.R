#!/usr/bin/env Rscript

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
                    "*Some ancestries do not have reported allele frequencies through gnomAD, and are thus indicated as 'N/A'.",
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
                width = 4,
                shinyWidgets::checkboxGroupButtons(
                    inputId   = ns("filters"),
                    label     = "Filters:",
                    choices   = c("Deleterious", "Conserved", "Kinase active"),
                    selected  = NULL,
                    status    = "primary",
                    justified = FALSE,
                    width     = "100%" ,
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon"),
                        no  = icon("remove", lib = "glyphicon")
                    )
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
        cols_to_keep <- c(
            "Variant (GrCh38)",
            "PD frequency",
            "Control frequency",
            "gnomAD allele frequency",
            "Region",
            "Functional consequence",
            "CADD",
            "rsID",
            "Clinical significance",
            "Clinical disease name",
            "Exon #",
            "cDNA change",
            "AA change",
            "Protein domain",
            "Conservation score",
            "Kinase activity (mean pRAB10/RAB10)"
        )

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
            dat <- dat[, ..cols_to_keep]
            sel <- input$filters %||% character(0)

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

            dat
        })

        output$table <- DT::renderDT({
            dat <- dat_rx()

            # Make first column (variant ID) clickable
            if (nrow(dat) > 0) {
                dat[[1]] <- paste0(
                    '<a href="#" class="variant-link" style="color: #0C8DC3; text-decoration: underline;">',
                    dat[[1]],
                    '</a>'
                )
            }

            # Columns to show in scientific notation
            sci_cols    <- intersect(c("PD frequency", "Control frequency", "gnomAD allele frequency"), colnames(dat))
            sci_targets <- match(sci_cols, colnames(dat)) - 1L

            DT::datatable(
                dat,
                extensions = "Buttons",
                rownames   = FALSE,
                escape     = FALSE,
                selection  = "none",
                options    = list(
                    dom          = "Blfrtip",
                    buttons      = c("copy", "csv", "excel", "pdf", "print"),
                    paging       = TRUE,
                    pageLength   = 25,
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
                                "       return num.toExponential(3);",
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
