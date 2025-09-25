#!/usr/bin/env Rscript

cdnaDiagramUI <- function(id) {
  ns <- NS(id)
  plotlyOutput(ns("cdnaDiagram"))
}

cdnaDiagramServer <- function(id, exon_colors, exon_positions, subdomain_gap) {
    # Number of exons
    num_exons <- length(exon_positions) - 1
    # Adjusted start positions of each exon, accounting for line segments 
    exon_positions_adj <- exon_positions + (0:(num_exons)) * subdomain_gap
    # Metadata used for the boxes for each protein exon
    exon_names <- paste("", 1:num_exons)
    # Add tooltips for each box
    exon_info <- paste("Exon", 1:num_exons)
    # Create data tables with metadata for each subdomain and exon
    exons <- data.frame(
        start = exon_positions_adj[-(num_exons+1)],
        end   = exon_positions_adj[-1],
        label = exon_names,
        color = exon_colors,
        info  = exon_info
    )

    moduleServer(id, function(input, output, session) {
        output$cdnaDiagram <- renderPlotly({
            # Start and end plot coordinate positions for each exon
            exons$start_adj <- exons$start + subdomain_gap/2
            exons$end_adj <- exons$end - subdomain_gap/2

            # Initialize plotly plot
            p <- plot_ly(source = "protein") %>%
                layout(
                    hovermode = "closest",
                    xaxis = list(showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE,
                                range = c(min(exons$start) - 50, max(exons$end) + 50)),
                    yaxis = list(showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE, range = c(-1, 1))
                )

            # Draw line segments connecting exons
            line_x <- c(
                exons$start_adj[1] - subdomain_gap,
                as.vector(rbind(exons$start_adj, exons$end_adj)),
                exons$end_adj[num_exons] + subdomain_gap
            )
            line_y <- rep(0, length(line_x))
            
            # Add line segments to the plot
            p <- add_trace(
                p,
                x = line_x,
                y = line_y,
                type = "scatter",
                mode = "lines",
                line = list(color = "black", width = 2),
                showlegend = FALSE,
                hoverinfo = "none"
            )

            for (i in 1:nrow(exons)) {
                # Rectangle coordinates
                x_rect <- c(exons$start_adj[i], exons$end_adj[i], exons$end_adj[i], exons$start_adj[i], exons$start_adj[i])
                y_rect <- c(-0.2, -0.2, 0.2, 0.2, -0.2)
                
                # Draw rectangle with tooltip
                p <- add_trace(
                    p,
                    x = x_rect,
                    y = y_rect,
                    type = "scatter",
                    mode = "lines",
                    fill = "toself",
                    fillcolor = exons$color[i],
                    line = list(color = "black"),
                    text = exons$info[i],
                    hoverinfo = "text",
                    hoveron = "fills",
                    showlegend = FALSE
                )
            }

            # Add starting cDNA position under first connector segment
            p <- add_annotations(
                p,
                x = 1,
                y = -0.3,
                text = "1",
                showarrow = FALSE,
                font = list(size = 14)
            )

            for (i in seq(1000, max(exon_positions) - 1, by = 1000)) {
                # Midpoint of connector segment
                x <- i + subdomain_gap * sum(exon_positions < i)
                
                # Add cDNA position label under the given connector segment
                label <- paste0(i)
                p <- add_annotations(
                    p,
                    x = x,
                    y = -0.3,
                    text = label,
                    showarrow = FALSE,
                    font = list(size = 14)
                )
            }
            p
        })
    })
}