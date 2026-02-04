#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

barChartUI <- function(id) {
    ns <- NS(id)
    tagList(
        # Sort order dropdown by either genomic coordinate or kinase activity
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
        # Bar chart output
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

            # Subset only variants with kinase activity, keeping only necessary columns
            dat <- variant_data[
                !is.na(variant_data[["Kinase activity (mean pRAB10/RAB10)"]]),
                c("Kinase activity (mean pRAB10/RAB10)", "cDNA change", "AA change", "Exon #", "Variant (GrCh38)"),
                with = FALSE
            ]

            print(variant_data[variant_data$`AA change` == "p.R1728H", ])

            # Aggregate by AA change, combining DNA/cDNA variants
            dat_aggregated <- dat[, .(
                `Kinase activity (mean pRAB10/RAB10)` = mean(`Kinase activity (mean pRAB10/RAB10)`, na.rm = TRUE),
                `cDNA change` = paste(unique(`cDNA change`), collapse = "<br>                         "),
                `Variant (GrCh38)` = paste(unique(`Variant (GrCh38)`), collapse = "<br>                         "),
                `Exon #` = first(`Exon #`)
            ), by = `AA change`]

            # Sort based on dropdown selection
            dat_aggregated <- switch(
                input$sort_order,
                coord    = dat_aggregated,
                activity = dat_aggregated[order(-dat_aggregated[["Kinase activity (mean pRAB10/RAB10)"]]), ]
            )

            # Convert AA change column to factor for plotting on horizontal axis
            dat_aggregated[, `AA change` := factor(`AA change`, levels = unique(dat_aggregated$`AA change`))]
            
            # Convert to data frame for plotly
            dat_aggregated <- as.data.frame(dat_aggregated)

            # Define color palette based on sort order
            if (input$sort_order == "coord") {
                color_palette <- domain_colors[dat_aggregated$`Exon #`]
            } else {
                activity_vals <- dat_aggregated$`Kinase activity (mean pRAB10/RAB10)`
                color_palette <- ifelse(
                    activity_vals > kinase_activation_threshold, "#8C4E9F",
                    ifelse(activity_vals < kinase_inactivation_threshold, "#34A270", "#E3E3E3")
                )
            }

            # Store genomic and cDNA changes to include in tooltip
            custom_data <- lapply(1:nrow(dat_aggregated), function(i) {
                list(dat_aggregated$`Variant (GrCh38)`[i], dat_aggregated$`cDNA change`[i])
            })

            # Make the bar plot
            p <- plot_ly(
                dat_aggregated,
                x = ~`AA change`,
                y = ~`Kinase activity (mean pRAB10/RAB10)`,
                type = "bar",
                marker = list(color = color_palette),
                customdata = custom_data,
                hovertemplate = "DNA change:    %{customdata[0]}<br>cDNA change:  %{customdata[1]}<br>AA change:       %{x}<br>Kinase activity: %{y:.3f}<extra></extra>"
            )

            # Add kinase activation threshold line and domain annotations if sorting by genomic coordinate
            if (input$sort_order == "coord") {
                line_shape <- list(
                    type = "line",
                    xref = "paper",
                    x0 = 0,
                    x1 = 1,
                    yref = "y",
                    y0 = 1.4,
                    y1 = 1.4,
                    line = list(color = "#8C4E9F", width = 2, dash = "dash")
                )
                line_annotation <- list(
                    xref = "paper",
                    x = 1.00,
                    y = 1.4,
                    text = "Kinase activating",
                    showarrow = FALSE,
                    font = list(color = "#8C4E9F", size = 14),
                    xanchor = "left",
                    yanchor = "middle"
                )

                exon_domain_map <- data.frame(
                    Exon = 1:sum(c(17, 2, 9, 3, 3, 4, 5, 8)),
                    Domain = rep(
                        c("ARM", "ANK", "LRR", "ROC", "COR-A", "COR-B", "KIN", "WD-40"),
                        times = c(17, 2, 9, 3, 3, 4, 5, 8)
                    )
                )
                dat_aggregated <- merge(dat_aggregated, exon_domain_map, by.x = "Exon #", by.y = "Exon", all.x = TRUE, sort = FALSE)
                dat_aggregated$`AA change` <- factor(dat_aggregated$`AA change`, levels = unique(dat_aggregated$`AA change`))
                domain_annotations <- lapply(unique(dat_aggregated$Domain), function(domain_name) {
                    domain_variants <- which(dat_aggregated$Domain == domain_name)
                    if (length(domain_variants) == 0) return(NULL)
                    midpoint <- mean(domain_variants) - 1
                    list(
                        x = midpoint,
                        y = max(dat_aggregated$`Kinase activity (mean pRAB10/RAB10)`) * 1.05,
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

            # Update plot layout
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
