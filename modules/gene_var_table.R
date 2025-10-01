#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

geneVarTableUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Disclaimer message
    div(
      "DISCLAIMER: Data displayed here should be used for research purposes clinical guidelines.",
      style = "color: red; font-style: italic; margin: 5px 0 15px 0;"
    ),
    
    # Row with ancestry dropdown + filter buttons
    fluidRow(
      column(
        width = 3,
        selectInput(
          inputId  = ns("dataset"),
          label    = "Choose ancestry:",
          choices  = NULL,
          width    = "100%"
        )
      ),
      column(
        width = 3,
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
    
    # Table
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

geneVarTableServer <- function(id, all_tables_cleaned) {
    moduleServer(id, function(input, output, session) {

        # Initialize dataset selector
        updateSelectInput(
            session, "dataset",
            choices  = names(all_tables_cleaned),
            selected = if ("Combined" %in% names(all_tables_cleaned)) "Combined" else names(all_tables_cleaned)[1]
        )

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
