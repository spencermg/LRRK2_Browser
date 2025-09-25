#!/usr/bin/env Rscript

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
                shinyWidgets::radioGroupButtons(
                    inputId = ns("filter"),
                    label = div(
                        "Filter variants:",
                        title = "For more options, please use the \"Search\" box on the right!"
                    ),
                    choices  = c("No filter", "Nonsynonymous", "Frameshift", "Stop gain/loss"),
                    selected = "No filter",
                    status   = "primary",
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon"),
                        no  = icon("remove", lib = "glyphicon")
                    )
                ),

                # DT output target for renderDT()
                DT::DTOutput(ns("table")),

                style = "margin: 12px 50px 50px 12px;"
            )
        )
    )
}

geneVarTableServer <- function(id, all_tables_cleaned) {
    moduleServer(id, function(input, output, session) {
        # Populate the selector after module is created
        updateSelectInput(
            session, 
            "dataset",
            choices  = names(all_tables_cleaned),
            selected = if ("Combined" %in% names(all_tables_cleaned)) "Combined" else names(all_tables_cleaned)[1]
        )

        # Reactive table with filter applied
        dat_rx <- reactive({
            dat <- all_tables_cleaned[[ req(input$dataset) ]]
            fc <- dat[["Functional consequence"]]
            if (is.null(input$filter) || input$filter == "No filter") return(dat)
            switch(
                input$filter,
                "Nonsynonymous"  = dat[ !(fc %in% c(".", "synonymous SNV")), , drop = FALSE],
                "Frameshift"     = dat[  fc %in% c("frameshift deletion", "frameshift insertion", "frameshift block substitution"), , drop = FALSE],
                "Stop gain/loss" = dat[  fc %in% c("stopgain", "stoploss"), , drop = FALSE],
                dat
            )
        })

        # Render the DT into the DTOutput above
        output$table <- DT::renderDT({
            dat <- dat_rx()
            DT::datatable(
                dat,
                extensions = "Buttons",
                options = list(
                    dom = "Blfrtip",
                    buttons = c("copy", "csv", "excel", "pdf", "print"),
                    paging = TRUE,
                    pageLength = 25,
                    lengthChange = TRUE,
                    lengthMenu = list(c(10,25,50,100,500,-1), c("10","25","50","100","500","All")),
                    scrollX = TRUE,
                    deferRender = TRUE
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
