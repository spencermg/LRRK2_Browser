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

                # Multi-select button filters
                shinyWidgets::checkboxGroupButtons(
                    inputId  = ns("filters"),
                    label    = "Filters:",
                    choices  = c("CADD > 20", "Kinase active"),
                    selected = NULL,                 # none selected by default
                    status   = "primary",
                    justified = FALSE,
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon"),
                        no  = icon("remove", lib = "glyphicon")
                    )
                ),

                # Reset button (clears selected filter buttons)
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

geneVarTableServer <- function(id, all_tables_cleaned) {
    moduleServer(id, function(input, output, session) {

        # Initialize dataset selector
        updateSelectInput(
            session, "dataset",
            choices  = names(all_tables_cleaned),
            selected = if ("Combined" %in% names(all_tables_cleaned)) "Combined" else names(all_tables_cleaned)[1]
        )

        # Reset button clears the filter buttons
        observeEvent(input$reset_filters, {
            shinyWidgets::updateCheckboxGroupButtons(session, "filters", selected = character(0))
        })

        # Base data (selected dataset)
        dat_base <- reactive({
            all_tables_cleaned[[ req(input$dataset) ]]
        })

        # Apply selected filters (can be both, one, or none)
        dat_rx <- reactive({
            dat <- dat_base()
            sel <- input$filters %||% character(0)

            # --- CADD > 20 ---
            if ("CADD > 20" %in% sel) {
                cadd <- dat[["CADD"]]
                # robust coercion to numeric; keep NA as FALSE in filter
                cadd_num <- suppressWarnings(as.numeric(cadd))
                keep <- !is.na(cadd_num) & cadd_num > 20
                dat <- dat[ keep, , drop = FALSE]
            }

            # --- Kinase active ---
            if ("Kinase active" %in% sel) {
                ka <- dat[["Kinase activity (mean pRAB10/RAB10)"]]
                # robust coercion to numeric; keep NA as FALSE in filter
                kinase_active_num <- suppressWarnings(as.numeric(ka))
                keep <- !is.na(kinase_active_num) & kinase_active_num > 1.40
                dat <- dat[ keep, , drop = FALSE]
            }

            dat
        })

        output$table <- DT::renderDT({
            dat <- dat_rx()

            # columns to show in scientific notation
            sci_cols    <- intersect(c("PD frequency", "Control frequency", "Gnomad allele frequency"), colnames(dat))
            sci_targets <- match(sci_cols, colnames(dat)) - 1L  # DataTables is 0-indexed

            DT::datatable(
                dat,
                extensions = "Buttons",
                options = list(
                    dom          = "Blfrtip",
                    buttons      = c("copy", "csv", "excel", "pdf", "print"),
                    paging       = TRUE,
                    pageLength   = 25,
                    lengthChange = TRUE,
                    lengthMenu   = list(c(10,25,50,100,500,-1), c("10","25","50","100","500","All")),
                    scrollX      = TRUE,
                    deferRender  = TRUE,
                    columnDefs = if (length(sci_targets)) list(
                        list(
                            targets   = sci_targets,
                            className = "dt-right",
                            render    = DT::JS(
                                "function(data, type, row, meta){",
                                "  if (data === null || data === undefined || data === '') {",
                                "    return (type === 'display') ? '' : null;",
                                "  }",
                                "  var num = Number(data);",
                                "  if (!isFinite(num)) {",
                                "    return (type === 'display') ? data : null;",
                                "  }",
                                "  if (type === 'display') {",
                                "    return num.toExponential(3);",
                                "  }",
                                "  return num;",
                                "}"
                            )
                        )
                    ) else NULL
                ),
                rownames = FALSE,
                escape   = FALSE
            ) |>
            DT::formatStyle(
                columns = colnames(dat),
                backgroundColor = "#FFFFFF",
                color = "black"
            )
        })
    })
}
