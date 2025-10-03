#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

barChartUI <- function(id) {
    ns <- NS(id)
    tagList(
        div(
            style = "display: inline-block; vertical-align: middle; margin-left: 10px;",
            "Sort by:",
            selectInput(
                ns("sort_order"), 
                label = NULL,
                choices = c("Genomic coordinate" = "coord", "Kinase activity" = "activity"),
                selected = "coord",
                width = "200px"
            )
        ),
        div(
            style = "overflow-x: auto; overflow-y: hidden; width: 100%;",
            plotlyOutput(ns("bar_chart"), height = "500px", width = "2000px")
        )
    )
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

barChartServer <- function(id, variant_data) {
    moduleServer(id, function(input, output, session) {
        output$bar_chart <- renderPlotly({
            req(variant_data)

            # Subset only variants with kinase activity
            dat <- variant_data[
                !is.na(variant_data[["Kinase activity (mean pRAB10/RAB10)"]]),
                c("Variant (GrCh38)", "Kinase activity (mean pRAB10/RAB10)", "cDNA change", "AA change"),
                with = FALSE
            ]

            # Sort based on dropdown selection
            dat <- switch(input$sort_order,
                coord    = dat,  # leave as-is (assumes already in genomic order)
                activity = dat[order(-dat[["Kinase activity (mean pRAB10/RAB10)"]]), ]
            )

            dat[, `Variant (GrCh38)` := factor(`Variant (GrCh38)`, levels = dat$`Variant (GrCh38)`)]
            
            # Convert to data frame for plotly
            dat_df <- as.data.frame(dat)

            custom_data <- lapply(1:nrow(dat_df), function(i) {
                list(dat_df$`cDNA change`[i], dat_df$`AA change`[i])
            })

            # Make the bar plot
            plot_ly(
                dat_df,
                x = ~`Variant (GrCh38)`,
                y = ~`Kinase activity (mean pRAB10/RAB10)`,
                type = "bar",
                marker = list(color = "#0C8DC3"),
                customdata = custom_data,
                hovertemplate = "%{x}<br>cDNA change: %{customdata[0]}<br>AA change: %{customdata[1]}<br>Kinase activity: %{y:.3f}<extra></extra>"
            ) %>%
            layout(
                title = NULL,
                xaxis = list(
                    title = "Variant",
                    tickangle = 45,
                    fixedrange = TRUE
                ),
                yaxis = list(
                    title = "Kinase activity",
                    fixedrange = TRUE
                ),
                margin = list(l = 100),
                showlegend = FALSE
            ) %>%
            config(
                displayModeBar = FALSE,
                scrollZoom = FALSE
            )
        })
    })
}