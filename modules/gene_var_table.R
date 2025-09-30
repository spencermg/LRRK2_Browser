#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

geneVarTableUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            div(
                # Dataset selector
                selectInput(
                    inputId  = ns("dataset"),
                    label    = "Choose ancestry:",
                    choices  = NULL,
                    width    = "250px"
                ),

                # Filter buttons
                shinyWidgets::checkboxGroupButtons(
                    inputId   = ns("filters"),
                    label     = "Filters:",
                    choices   = c("Deleterious", "Conserved", "Kinase active"),
                    selected  = NULL,
                    status    = "primary",
                    justified = FALSE,
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon"),
                        no  = icon("remove", lib = "glyphicon")
                    )
                ),

                # Reset button
                actionButton(
                    inputId = ns("reset_filters"),
                    label   = "Reset filters",
                    class   = "btn btn-secondary",
                    icon    = icon("undo")
                ),

                # Table
                DT::DTOutput(ns("table")),

                style = "margin: 12px 50px 50px 12px;"
            )
        )
    )
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

geneVarTableServer <- function(id, all_tables_cleaned, kinase_activation_threshold) {
    moduleServer(id, function(input, output, session) {

        # Initialize dataset selector
        updateSelectInput(
            session, "dataset",
            choices  = names(all_tables_cleaned),
            selected = if ("Combined" %in% names(all_tables_cleaned)) "Combined" else names(all_tables_cleaned)[1]
        )

        # Add reset button
        observeEvent(input$reset_filters, {
            shinyWidgets::updateCheckboxGroupButtons(session, "filters", selected = character(0))
        })

        # Apply selected filters
        dat_rx <- reactive({
            dat <- all_tables_cleaned[[ req(input$dataset) ]]
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
            # Grab dataset with current filters applied
            dat <- dat_rx()

            # Columns to show in scientific notation
            sci_cols    <- intersect(c("PD frequency", "Control frequency", "Gnomad allele frequency"), colnames(dat))
            sci_targets <- match(sci_cols, colnames(dat)) - 1L

            # Populate data table
            DT::datatable(
                dat,
                extensions = "Buttons",
                rownames   = FALSE,
                escape     = FALSE,
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
    })
}
