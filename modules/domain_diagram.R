#!/usr/bin/env Rscript

# =========================================================================
# AUXILIARY FUNCTIONS
# =========================================================================

get_text_color <- function(hex) {
    rgb <- col2rgb(hex) / 255
    brightness <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
    if (brightness > 0.5) "black" else "white"
}


# =========================================================================
# UI FUNCTION
# =========================================================================

diagramUI <- function(id) {
  ns <- NS(id)
  plotlyOutput(ns("diagram"))
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

diagramServer <- function(
    id,
    variant_data,
    domains,
    domain_positions,
    subdomain_gap,
    variants,
    mode,
    top_n,
    variant_bus,
    max_stack_size = 6,     # Max num labels to stack on each other
    min_label_sep_px = 50,  # Num pixel gap below which labels stack vertically
    label_lift = 0.07,      # Vertical lift per stacked label
    label_offset = 1        # Additional vertical offset above lollipop for first label
) {
    # Map position to x-value accounting for connecting line segments
    map_pos_to_x <- function(pos_vec) {
        subdomain_idx <- findInterval(pos_vec, domain_positions, left.open = TRUE)
        pos_vec + (subdomain_idx - 0.5) * subdomain_gap
    }

    # Find variant(s) associated with a given AA/cDNA change
    find_variant_ids <- function(change, type) {
        if (type == "cDNA") {
            match_row <- variant_data[`cDNA change` == change]
        } else {
            match_row <- variant_data[`AA change` == change]
        }
        if (nrow(match_row) == 0) return(NULL)
        match_row$`Variant (GrCh38)`
    }

    moduleServer(id, function(input, output, session) {
        # Create new columns converting cDNA/AA IDs to positions
        if (!"aa_pos" %in% names(variant_data)) {
            variant_data[, aa_pos := as.integer(gsub("[^0-9]", "", `AA change`))]
        }
        if (!"cdna_pos" %in% names(variant_data)) {
            variant_data[, cdna_pos := as.integer(gsub("[^0-9]", "", `cDNA change`))]
        }

        # Full x-range and interactive window
        x_full   <- c(
            min(domains$start), 
            max(domains$end)
        )
        x_window <- reactiveVal(x_full)

        # Find all variants within the diagram and precompute their x/y coords
        v_all <- local({
            v <- variants[
                is.finite(variants$pos) &
                variants$pos >= min(domain_positions) &
                variants$pos <= max(domain_positions),
                , drop = FALSE
            ]
            if (nrow(v) == 0) return(NULL)

            v$x <- map_pos_to_x(v$pos)
            v$y_top <- 0.6
            v
        })

        # Recompute Top N variants dynamically based on current x-window
        visible_variants <- reactive({
            if (is.null(v_all)) return(NULL)
            x_range <- x_window()
            v <- v_all[v_all$x >= x_range[1] & v_all$x <= x_range[2], , drop = FALSE]
            if (nrow(v) == 0) return(v)
            v <- v[order(-ifelse(is.finite(v$value), v$value, -Inf), v$pos), , drop = FALSE]
            head(v, top_n)
        })

        # Track zoom/pan/autoscale/double-click to update x-window
        observeEvent(plotly::event_data("plotly_relayout", source = id), {
            ev <- plotly::event_data("plotly_relayout", source = id)

            # Autorange gives full spectrum of variants
            if (isTRUE(ev[["xaxis.autorange"]])) {
                x_window(x_full)
                return()
            }

            # Update x values for zooming and/or panning
            x0 <- ev[["xaxis.range[0]"]] 
            x1 <- ev[["xaxis.range[1]"]]
            if (!is.null(x0) && !is.null(x1)) {
                x_window(c(x0, x1))
                return()
            }

            # If a 2-length vector is returned in xaxis.range
            xr <- ev[["xaxis.range"]]
            if (!is.null(xr) && length(xr) == 2) {
                x_window(c(xr[1], xr[2]))
                return()
            }
        }, ignoreInit = TRUE)
        observeEvent(plotly::event_data("plotly_doubleclick", source = id), {
            x_window(x_full)
        }, ignoreInit = TRUE)

        # Respond to users clicking on a variant
        observeEvent(
            plotly::event_data("plotly_click", source = id, priority = "event"),
            ignoreInit = TRUE,
            {
                click_data <- plotly::event_data("plotly_click", source = id)
                req(click_data$customdata)

                # Find all genomic variants at the position clicked by the user
                clicked_pos <- as.numeric(unlist(click_data$customdata)[1])
                rows <- if (mode == "protein") {
                    variant_data[aa_pos == clicked_pos]
                } else {
                    variant_data[cdna_pos == clicked_pos]
                }
                ids <- rows$`Variant (GrCh38)`

                # Prompt user to select which variant to show in the popup if more than one option
                if (length(ids) == 1) {
                    variant_bus$publish(list(variant_id = ids[1]))
                } else {
                    showModal(
                        modalDialog(
                            title = "Multiple genomic variants found",
                            selectInput(
                                session$ns("variant_choice"),
                                "Select variant:",
                                choices = ids
                            ),
                            footer = tagList(
                                modalButton("Cancel"),
                                actionButton(session$ns("confirm_variant"), "Show Variant")
                            )
                        )
                    )
                }
            }
        )
        observeEvent(input$confirm_variant, {
            req(input$variant_choice)
            removeModal()
            variant_bus$publish(
                list(variant_id = input$variant_choice)
            )
        })

        output$diagram <- renderPlotly({
            # Start/end for the boxes, adjusted for connecting line segments
            domains$start_adj <- domains$start + subdomain_gap/2
            domains$end_adj   <- domains$end   - subdomain_gap/2

            # Find height of y-axis to fit all labels
            y_max <- max(1, 0.6 + (label_offset + max(1, max_stack_size)) * label_lift)

            # Set up base plot
            p <- plotly::plot_ly(source = id) %>%
                plotly::layout(
                    uirevision = "topN",
                    hovermode  = "closest",
                    showlegend = FALSE,
                    dragmode   = "pan",
                    xaxis = list(
                        showticklabels = FALSE, 
                        zeroline = FALSE, 
                        showgrid = FALSE,
                        range = x_window(),
                        fixedrange = FALSE
                    ),
                    yaxis = list(
                        showticklabels = FALSE, 
                        zeroline = FALSE, 
                        showgrid = FALSE,
                        range = c(-1, y_max), 
                        fixedrange = TRUE
                    ),
                    margin = list(t = 30)
                ) %>%
                plotly::config(
                    displaylogo = FALSE,
                    scrollZoom  = FALSE,
                    modeBarButtonsToRemove = c("zoom2d", "lasso2d", "select2d")
                )

            # Add connector line segments
            line_x <- c(
                domains$start_adj[1] - subdomain_gap,
                as.vector(rbind(domains$start_adj, domains$end_adj)),
                domains$end_adj[length(domains$label)] + subdomain_gap
            )
            p <- plotly::add_trace(
                p,
                x = line_x, 
                y = rep(0, length(line_x)),
                type = "scatter", 
                mode = "lines",
                line = list(color = "black", width = 2),
                hoverinfo = "none", 
                showlegend = FALSE
            )

            # Domain rectangles + labels
            for (i in seq_len(nrow(domains))) {
                # Add boxes for each domain
                p <- plotly::add_trace(
                    p,
                    x = c(
                        domains$start_adj[i], domains$end_adj[i],
                        domains$end_adj[i], domains$start_adj[i],
                        domains$start_adj[i]
                    ),
                    y = c(-0.2, -0.2, 0.2, 0.2, -0.2),
                    type = "scatter", 
                    mode = "lines",
                    fill = "toself", 
                    fillcolor = domains$color[i],
                    line = list(color = "black"),
                    text = domains$info[i], 
                    hoverinfo = "text", 
                    hoveron = "fills",
                    showlegend = FALSE
                )

                # Add protein domain labels
                p <- plotly::add_annotations(
                    p,
                    x = (domains$start_adj[i] + domains$end_adj[i]) / 2,
                    y = 0, 
                    text = if (mode == "protein") domains$label[i] else "",
                    showarrow = FALSE, 
                    font = list(size = 14, color = get_text_color(domains$color[i]))
                )
            }

            # Add positional labels along the horizontal axis
            if (mode == "protein") {
                p <- plotly::add_annotations(
                    p, 
                    x = 1, 
                    y = -0.3, 
                    text = "1",
                    showarrow = FALSE, 
                    font = list(size = 14)
                )
                for (i in seq_len(nrow(domains))) {
                    p <- plotly::add_annotations(
                        p, 
                        x = domains$end_adj[i] + subdomain_gap/2, 
                        y = -0.3,
                        text = paste0(domain_positions[i + 1]),
                        showarrow = FALSE, 
                        font = list(size = 14)
                    )
                }
            }
            else if (mode == "cDNA") {
                # Add cDNA position labels every 1000 bases
                for (i in seq(1000, max(domain_positions) - 1, by = 1000)) {
                    # Midpoint of connector segment
                    x <- i + subdomain_gap * sum(domain_positions < i)
                    
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
            }

            # Add variant lollipops (grouped by position)
            v_visible <- visible_variants()

            if (!is.null(v_visible) && nrow(v_visible) > 0) {
                # Aggregate visible variants by position
                v_grouped <- v_visible[, .(
                    count = .N
                ), by = .(pos, x)]
                v_grouped$y_top <- 0.6

                # Define tooltip
                tooltip <- vapply(v_grouped$pos, function(p) {
                    if (mode == "protein"){
                        rows <- variant_data[aa_pos == p]
                    } else {
                        rows <- variant_data[cdna_pos == p]
                    }
                    paste0(
                        "<b>", nrow(rows), " variant", if (nrow(rows) != 1) "s" else "", " at position ", p, "</b><br><br>",
                        paste(
                            apply(rows, 1, function(r) {
                                paste0(
                                    "<b>DNA:</b> ", r[["Variant (GrCh38)"]], "<br>",
                                    "<b>cDNA:</b> ", r[["cDNA change"]], "<br>",
                                    "<b>AA:</b> ", r[["AA change"]]
                                )
                            }),
                            collapse = "<br><br>"
                        )
                    )
                }, FUN.VALUE = character(1))

                # Add line segments for the lollipop stems
                p <- plotly::add_segments(
                    p,
                    x = v_grouped$x, xend = v_grouped$x,
                    y = 0.2, yend = v_grouped$y_top,
                    line = list(color = "#444444", width = 1),
                    hoverinfo = "none",
                    showlegend = FALSE
                )

                # Include both position and tooltip into customdata
                custom_data <- Map(function(pos, txt) {
                    list(pos, txt)
                }, v_grouped$pos, tooltip)

                # Add markers at the top of the lollipops
                p <- plotly::add_trace(
                    p,
                    x = v_grouped$x,
                    y = v_grouped$y_top,
                    type = "scatter",
                    mode = "markers",
                    marker = list(
                        size = 9 + 2 * log1p(v_grouped$count),
                        color = "#444444",
                        line = list(color = "black", width = 0.5)
                    ),
                    customdata = custom_data,
                    hovertemplate = "%{customdata[1]}<extra></extra>",
                    showlegend = FALSE
                )

                # Build labels grouped by position
                label_rows <- v_grouped[, {
                    labs <- unique(v_visible[pos == .BY$pos]$label)
                    .(label = labs)
                }, by = .(pos, x)][order(x)]

                # Find minimum soacing between labels in coordinate units
                out_id <- session$ns("diagram")
                plot_width_px <- session$clientData[[paste0("output_", out_id, "_width")]]
                if (is.null(plot_width_px) || plot_width_px <= 0) plot_width_px <- 800
                x_range <- x_window()
                x_spacing <- (x_range[2] - x_range[1]) * (min_label_sep_px / plot_width_px)

                # Assign each label to the "level" it will be at vertically
                max_levels <- max_stack_size
                last_x <- rep(-Inf, max_levels)
                label_level <- rep(NA_integer_, nrow(label_rows))
                for (i in seq_len(nrow(label_rows))) {
                    for (L in seq_len(max_levels)) {
                        if (label_rows$x[i] - last_x[L] >= x_spacing) {
                            label_level[i] <- L - 1
                            last_x[L] <- label_rows$x[i]
                            break
                        }
                    }
                }

                # Remove labels that couldn't be placed
                label_rows <- label_rows[!is.na(label_level)]

                # Compute y positions
                label_rows$y <- 0.6 + (label_level + label_offset) * label_lift

                # Set maximum number of variant labels to stack vertically
                max_labels_per_pos <- max(v_grouped$count)
                y_max <<- max(1, 0.6 + (label_offset + max_labels_per_pos) * label_lift)

                # Prevent overflow
                label_rows$y <- pmin(y_max - 0.02, label_rows$y)

                # Add labels
                p <- plotly::add_trace(
                    p,
                    x = label_rows$x,
                    y = label_rows$y,
                    type = "scatter",
                    mode = "text",
                    text = label_rows$label,
                    textposition = "middle center",
                    textfont = list(size = 12),
                    hoverinfo = "skip",
                    showlegend = FALSE,
                    cliponaxis = FALSE
                )
            }

            p
        })
    })
}
