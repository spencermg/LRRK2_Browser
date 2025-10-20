#!/usr/bin/env Rscript

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
    domains,
    domain_positions,
    subdomain_gap,
    variants,
    mode,
    top_n,
    max_labels = 6,          # Max num labels to stack on each other
    min_label_sep_px = 50,   # Num pixel gap below which labels stack vertically
    label_lift = 0.07,
    label_offset = 1,
    y_padding = 0.07
) {
    # Map position to x-value accounting for connecting line segments
    map_pos_to_x <- function(pos_vec) {
        subdomain_idx <- findInterval(pos_vec, domain_positions, left.open = TRUE)
        pos_vec + (subdomain_idx - 0.5) * subdomain_gap
    }

    moduleServer(id, function(input, output, session) {
        # Full x-range and reactive current window
        x_full   <- c(
            min(domains$start), 
            max(domains$end)
        )
        x_window <- reactiveVal(x_full)

        v_all <- local({
            # Keep only the variants that fit in the range of the diagram
            v <- variants[
                is.finite(variants$pos) &
                variants$pos >= min(domain_positions) &
                variants$pos <= max(domain_positions),
                , drop = FALSE
            ]
            if (nrow(v) == 0) return(NULL)

            # Precompute variant coords/heights upfront
            v$x <- map_pos_to_x(v$pos)
            v$y_top <- 0.6
            v
        })

        # Recompute visible Top N variants whenever the window changes
        v_visible <- reactive({
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

        output$diagram <- renderPlotly({
            # Adjusted start/end for the boxes
            domains$start_adj <- domains$start + subdomain_gap/2
            domains$end_adj   <- domains$end   - subdomain_gap/2

            # Find height of y-axis to fit all labels
            y_max <- max(1, 0.6 + (label_offset + max(0, max_labels - 1)) * label_lift + y_padding)

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

            get_text_color <- function(hex) {
                rgb <- col2rgb(hex) / 255
                brightness <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
                if (brightness > 0.5) "black" else "white"
            }

            # Domain rectangles + labels
            for (i in seq_len(nrow(domains))) {
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
                # Add domain boundary tick labels
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


            # Lollipops: only Top N within current window (with label collision avoidance)
            vv <- v_visible()
            if (!is.null(vv) && nrow(vv) > 0) {
                # stems
                x_seg <- as.numeric(t(cbind(vv$x, vv$x, NA)))
                y_seg <- as.numeric(t(cbind(rep(0.2, nrow(vv)), vv$y_top, NA)))
                p <- plotly::add_trace(
                    p, 
                    x = x_seg, 
                    y = y_seg,
                    type = "scatter", 
                    mode = "lines",
                    line = list(color = "#444444", width = 1),
                    hoverinfo = "none", 
                    showlegend = FALSE
                )
                # markers
                p <- plotly::add_trace(
                    p,
                    x = vv$x, 
                    y = vv$y_top,
                    type = "scatter", 
                    mode = "markers",
                    marker = list(size = 9, color = vv$color, line = list(color = "black", width = 0.5)),
                    hovertemplate = "%{text}<br>Pos: %{customdata}<extra></extra>",
                    text = vv$label, 
                    customdata = vv$pos,
                    showlegend = FALSE
                )

                # ---- Stacked labels (uses tunables above) ----
                out_id <- session$ns("diagram")
                w_px <- session$clientData[[paste0("output_", out_id, "_width")]]
                if (is.null(w_px) || w_px <= 0) w_px <- 800
                xr <- x_window()
                thr_x <- (xr[2] - xr[1]) * (min_label_sep_px / w_px)

                ord <- order(vv$x, vv$y_top)
                last_x <- rep(-Inf, max_labels)
                labels  <- rep(NA_integer_, nrow(vv))

                for (i in ord) {
                    for (L in seq_len(max_labels)) {
                        if (vv$x[i] - last_x[L] >= thr_x) {
                            labels[i] <- L - 1
                            last_x[L] <- vv$x[i]
                            break
                        }
                    }
                    # if not placed, label omitted
                }

                lab <- vv[!is.na(labels), , drop = FALSE]
                if (nrow(lab) > 0) {
                    lab_y <- lab$y_top + (labels[!is.na(labels)] + label_offset) * label_lift
                    lab_y <- pmin(y_max - y_padding/2, lab_y)  # clamp to top

                    p <- plotly::add_trace(
                        p,
                        x = lab$x, 
                        y = lab_y,
                        type = "scatter",
                        mode = "text",
                        text = lab$label, 
                        textposition = "middle center",
                        textfont = list(size = 12),
                        hoverinfo = "skip",
                        showlegend = FALSE,
                        cliponaxis = FALSE
                    )
                }
            }

            p
        })
    })
}
