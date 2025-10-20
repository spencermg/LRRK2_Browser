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

barChartServer <- function(id, variant_data, kinase_activation_threshold, kinase_inactivation_threshold, domain_colors) {
    moduleServer(id, function(input, output, session) {
        output$bar_chart <- renderPlotly({
            req(variant_data)

            # Subset only variants with kinase activity
            dat <- variant_data[
                !is.na(variant_data[["Kinase activity (mean pRAB10/RAB10)"]]),
                c("AA change", "Kinase activity (mean pRAB10/RAB10)", "cDNA change", "AA change", "Exon #"),
                with = FALSE
            ]

            # Sort based on dropdown selection
            dat <- switch(input$sort_order,
                coord    = dat,  # leave as-is (assumes already in genomic order)
                activity = dat[order(-dat[["Kinase activity (mean pRAB10/RAB10)"]]), ]
            )

            dat[, `AA change` := factor(`AA change`, levels = dat$`AA change`)]
            
            # Convert to data frame for plotly
            dat_df <- as.data.frame(dat)

            custom_data <- lapply(1:nrow(dat_df), function(i) {
                list(dat_df$`cDNA change`[i], dat_df$`AA change`[i])
            })

            if (input$sort_order == "coord") {
                color_palette <- domain_colors[dat_df$`Exon #`]
            } else {
                activity_vals <- dat_df$`Kinase activity (mean pRAB10/RAB10)`
                color_palette <- ifelse(
                    activity_vals > kinase_activation_threshold, "#8c4e9f",
                    ifelse(activity_vals < kinase_inactivation_threshold, "#34a270", "#e3e4e3")
                )
            }

            # Make the bar plot
            p <- plot_ly(
                dat_df,
                x = ~`AA change`,
                y = ~`Kinase activity (mean pRAB10/RAB10)`,
                type = "bar",
                marker = list(color = color_palette),
                customdata = custom_data,
                hovertemplate = "%{x}<br>cDNA change: %{customdata[0]}<br>AA change: %{customdata[1]}<br>Kinase activity: %{y:.3f}<extra></extra>"
            )

            if (input$sort_order == "coord") {
                line_shape <- list(
                    type = "line",
                    xref = "paper",
                    x0 = 0,
                    x1 = 1,
                    yref = "y",
                    y0 = 1.4,
                    y1 = 1.4,
                    line = list(color = "#8c4e9f", width = 2, dash = "dash")
                )

                line_annotation <- list(
                    xref = "paper",
                    x = 1.00,
                    y = 1.4,
                    text = "Kinase activating",
                    showarrow = FALSE,
                    font = list(color = "#8c4e9f", size = 14),
                    xanchor = "left",
                    yanchor = "middle"
                )

                exon_domain_map <- data.frame(
                    Exon = 1:sum(c(17, 2, 9, 3, 3, 4, 5, 8)),
                    Domain = rep(c("ARM", "ANK", "LRR", "ROC", "COR-A", "COR-B", "KIN", "WD-40"),
                                times = c(17, 2, 9, 3, 3, 4, 5, 8))
                )

                dat_df <- merge(dat_df, exon_domain_map, by.x = "Exon #", by.y = "Exon", all.x = TRUE, sort = FALSE)

                dat_df$`AA change` <- factor(dat_df$`AA change`, levels = dat_df$`AA change`)

                domain_annotations <- lapply(unique(dat_df$Domain), function(domain_name) {
                    domain_variants <- which(dat_df$Domain == domain_name)
                    if (length(domain_variants) == 0) return(NULL)
                    midpoint <- mean(domain_variants) - 1
                    list(
                        x = midpoint,
                        y = max(dat_df$`Kinase activity (mean pRAB10/RAB10)`) * 1.05,
                        text = domain_name,
                        showarrow = FALSE,
                        font = list(size = 14, color = "black"),
                        xref = "x",
                        yref = "y",
                        xanchor = "center",
                        yanchor = "bottom"
                    )
                })
            } else {
                line_shape <- NULL
                line_annotation <- NULL
                domain_annotations <- NULL
            }

            p <- p %>%
                layout(
                    title = NULL,
                    xaxis = list(
                        title = "Variant",
                        tickangle = 45,
                        fixedrange = TRUE
                    ),
                    yaxis = list(
                        title = "Kinase activity (mean pRAB10/RAB10)",
                        fixedrange = TRUE
                    ),
                    margin = list(l = 100, r = 120, t = 20),
                    showlegend = FALSE,
                    shapes = line_shape,
                    annotations = c(list(line_annotation), domain_annotations)
                ) %>%
                config(
                    displayModeBar = FALSE,
                    scrollZoom = FALSE
                )
            p
        })
    })
}