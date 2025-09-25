#!/usr/bin/env Rscript

proteinDiagramUI <- function(id) {
  ns <- NS(id)
  plotlyOutput(ns("proteinDiagram"))
}

proteinDiagramServer <- function(id, protein_domain_colors, protein_domain_positions, subdomain_gap) {
    # Number of protein subdomains
    num_protein_domains <- length(protein_domain_positions) - 1
    # Adjusted start positions of each protein subdomain and exon, accounting for line segments 
    protein_domain_positions_adj <- protein_domain_positions + (0:(num_protein_domains)) * subdomain_gap
    # Metadata used for the boxes for each protein subdomain
    protein_domain_names <- c("ARM", "ANK", "LRR", "ROC", "COR-A", "COR-B", "KIN", "WD-40")
    # Add tooltips for each box
    ### TODO: Replace info tooltips with something meaningful
    protein_domain_info <- rep("Protein subdomain", num_protein_domains)
    # Create data tables with metadata for each subdomain and exon
    protein_domains <- data.frame(
        start = protein_domain_positions_adj[-(num_protein_domains+1)],
        end   = protein_domain_positions_adj[-1],
        label = protein_domain_names,
        color = protein_domain_colors,
        info  = protein_domain_info
    )

    moduleServer(id, function(input, output, session) {
        output$proteinDiagram <- renderPlotly({
            # Start and end plot coordinate positions for each subdomain
            protein_domains$start_adj <- protein_domains$start + subdomain_gap/2
            protein_domains$end_adj <- protein_domains$end - subdomain_gap/2

            # Initialize plotly plot
            p <- plot_ly(source = "protein") %>%
                layout(
                    hovermode = "closest",
                    xaxis = list(showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE,
                                range = c(min(protein_domains$start) - 50, max(protein_domains$end) + 50)),
                    yaxis = list(showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE, range = c(-1, 1))
                )

            # Draw line segments connecting subdomains
            line_x <- c(
                protein_domains$start_adj[1] - subdomain_gap,
                as.vector(rbind(protein_domains$start_adj, protein_domains$end_adj)),
                protein_domains$end_adj[num_protein_domains] + subdomain_gap
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

            for (i in 1:nrow(protein_domains)) {
                # Rectangle coordinates
                x_rect <- c(protein_domains$start_adj[i], protein_domains$end_adj[i], protein_domains$end_adj[i], protein_domains$start_adj[i], protein_domains$start_adj[i])
                y_rect <- c(-0.2, -0.2, 0.2, 0.2, -0.2)
                
                # Draw rectangle with tooltip
                p <- add_trace(
                    p,
                    x = x_rect,
                    y = y_rect,
                    type = "scatter",
                    mode = "lines",
                    fill = "toself",
                    fillcolor = protein_domains$color[i],
                    line = list(color = "black"),
                    text = protein_domains$info[i],
                    hoverinfo = "text",
                    hoveron = "fills",
                    showlegend = FALSE
                )
                
                # Add subdomain label
                p <- add_annotations(
                    p,
                    x = (protein_domains$start_adj[i] + protein_domains$end_adj[i])/2,
                    y = 0,
                    text = protein_domains$label[i],
                    showarrow = FALSE,
                    font = list(size = 14)
                )
            }

            # Add starting amino acid position under first connector segment
            p <- add_annotations(
                p,
                x = 1,
                y = -0.3,
                text = "1",
                showarrow = FALSE,
                font = list(size = 14)
            )

            for (i in 1:nrow(protein_domains)) {
                # Midpoint of connector segment
                x_mid <- protein_domains$end_adj[i] + subdomain_gap / 2
                
                # Add amino acid position label under the given connector segment
                label <- paste0(protein_domain_positions[i+1])
                p <- add_annotations(
                    p,
                    x = x_mid,
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