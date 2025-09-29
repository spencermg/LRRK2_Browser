#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

annotationSummaryTableUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("annotationSummaryTable"))
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

annotationSummaryTableServer <- function(id, all_tables_cleaned) {
    variantTable <- all_tables_cleaned$Combined

    # Define exonic variant annotation summary table to display in the UI
    variantTable_exonic <- variantTable[variantTable$Region == "exonic"]
    variantTable_exonic.counts <- variantTable_exonic[, .N, by = `Functional consequence`]
    colnames(variantTable_exonic.counts) <- c(
        "Exonic variant functional consequence",
        paste0("Count (Total: ", dim(variantTable_exonic)[1], ")")
    )

    # Define noncoding variant annotation summary table to display in the UI
    variantTable_noncoding <- variantTable[variantTable$Region != "exonic"]
    variantTable_noncoding.counts <- variantTable_noncoding[, .N, by = Region]
    colnames(variantTable_noncoding.counts) <- c(
        "Noncoding variant category",
        paste0("Count (Total: ", dim(variantTable_noncoding)[1], ")")
    )

    moduleServer(id, function(input, output, session) {
        output$annotationSummaryTable <- renderUI({
            total_count <- dim(variantTable_exonic)[1]
            tagList(
                div(
                    # Left table for exonic variants
                    column(
                        width = 6,
                        renderDT({
                            datatable(
                                variantTable_exonic.counts,
                                rownames = F,
                                selection = 'none',
                                options = list(
                                    paging = F,
                                    dom = 't'
                                )
                            ) %>% formatStyle(
                                columns = "Count",
                                valueColumns = "Exonic variant functional consequence",
                                target = 'cell',
                                color = "black"
                            ) %>% formatStyle(columns="Exonic variant functional consequence", backgroundColor = "white", color = "black")
                        })
                    ),
                    # Right table for noncoding variants
                    column(
                        width = 6,
                        renderDT({
                            datatable(
                                variantTable_noncoding.counts,
                                rownames = F,
                                selection = 'none',
                                options = list(
                                    paging = F,
                                    dom = 't'
                                )
                            ) %>% formatStyle(
                                columns = "Count",
                                valueColumns = "Noncoding variant category",
                                target = 'cell',
                                color = "black"
                            ) %>% formatStyle(columns="Noncoding variant category", backgroundColor = "white", color = "black")
                        })
                    )
                )
            )
        })
    })
}