#!/usr/bin/env Rscript

proteinDiagramUI <- function(id) {
  ns <- NS(id)
  plotlyOutput(ns("proteinDiagram"))
}

proteinDiagramServer <- function(
    id,
    protein_domain_colors,
    protein_domain_positions,
    subdomain_gap,
    variants,
    top_n,
    max_label_lanes = 6,
    min_label_sep_px = 50,
    lane_lift = 0.07,
    lane_offset = 1,
    y_padding = 0.07         # extra headroom at top of plot
) {
    num_protein_domains <- length(protein_domain_positions) - 1
    protein_domain_positions_adj <- protein_domain_positions + (0:num_protein_domains) * subdomain_gap
    protein_domain_names <- c("ARM","ANK","LRR","ROC","COR-A","COR-B","KIN","WD-40")
    protein_domain_tooltips <- rep("Protein subdomain", num_protein_domains)

    protein_domains <- data.frame(
        start = protein_domain_positions_adj[-(num_protein_domains + 1)],
        end   = protein_domain_positions_adj[-1],
        label = protein_domain_names,
        color = protein_domain_colors,
        info  = protein_domain_tooltips
    )

    # Map AA position -> x on gapped axis
    map_aa_to_x <- function(pos_vec) {
        subdomain_idx <- findInterval(pos_vec, protein_domain_positions, left.open = TRUE)
        pos_vec + subdomain_idx * subdomain_gap
    }

    moduleServer(id, function(input, output, session) {
        # Full x-range and reactive current window (updates on zoom/pan)
        x_full <- c(min(protein_domains$start) - 50, max(protein_domains$end) + 50)
        xwin   <- reactiveVal(x_full)

        # Precompute variant coords/heights once; we’ll subset reactively
        v_all <- local({
            if (is.null(variants) || nrow(variants) == 0) return(NULL)
            v <- variants[
                is.finite(variants$aa_pos) &
                variants$aa_pos >= min(protein_domain_positions) &
                variants$aa_pos <= max(protein_domain_positions),
                , drop = FALSE
            ]

            if (nrow(v) == 0) return(NULL)

            v$x <- map_aa_to_x(v$aa_pos)
            if (!("value" %in% names(v))) v$value <- 1

            rng <- range(v$value[is.finite(v$value)])
            if (!all(is.finite(rng)) || diff(rng) == 0) {
                v$y_top <- 0.6
            } else {
                v$y_top <- 0.5 + 0.25 * ((v$value - rng[1]) / diff(rng))  # 0.5–0.75
            }

            if (!("label" %in% names(v))) v$label <- paste0("AA ", v$aa_pos)
            if (!("color" %in% names(v))) v$color <- "red"
            v
        })

        # Recompute visible Top N variants whenever the window changes
        vis_v <- reactive({
            if (is.null(v_all)) return(NULL)
            rng <- xwin()
            v <- v_all[v_all$x >= rng[1] & v_all$x <= rng[2], , drop = FALSE]
            if (nrow(v) == 0) return(v)
            v <- v[order(-ifelse(is.finite(v$value), v$value, -Inf), v$aa_pos), , drop = FALSE]
            head(v, top_n)
        })

        # Track zoom/pan/autoscale to update x-window
        observeEvent(plotly::event_data("plotly_relayout", source = "protein"), {
            ev <- plotly::event_data("plotly_relayout", source = "protein")

            # 1) Autoscale / reset: autorange -> full view
            if (isTRUE(ev[["xaxis.autorange"]])) {
                xwin(x_full)
                return()
            }

            # 2) Explicit numeric ranges (zoom/pan via +/- or drag)
            x0 <- ev[["xaxis.range[0]"]]; x1 <- ev[["xaxis.range[1]"]]
            if (!is.null(x0) && !is.null(x1)) {
                xwin(c(x0, x1))
                return()
            }

            # 3) Some events send a 2-length vector instead
            xr <- ev[["xaxis.range"]]
            if (!is.null(xr) && length(xr) == 2) {
                xwin(c(xr[1], xr[2]))
                return()
            }
        }, ignoreInit = TRUE)

        # Double-click autoscale also restores full view
        observeEvent(plotly::event_data("plotly_doubleclick", source = "protein"), {
            xwin(x_full)
        }, ignoreInit = TRUE)

        output$proteinDiagram <- renderPlotly({
            # Adjusted start/end for the boxes
            protein_domains$start_adj <- protein_domains$start + subdomain_gap/2
            protein_domains$end_adj   <- protein_domains$end   - subdomain_gap/2

            # --- NEW: compute how tall the y-axis must be to fit all lanes ---
            y_top_max <- if (!is.null(v_all)) max(v_all$y_top, na.rm = TRUE) else 0.6
            # head height + (lane offset + max lanes - 1)*lift + padding; at least 1
            y_max <- max(1, y_top_max + (lane_offset + max(0, max_label_lanes - 1)) * lane_lift + y_padding)

            # Base plot (legend removed, vertical locked, no scroll zoom; only +/- buttons)
            p <- plotly::plot_ly(source = "protein") %>%
                plotly::layout(
                    uirevision = "protein-topN",       # preserve UI state across redraws
                    hovermode  = "closest",
                    showlegend = FALSE,
                    dragmode   = "pan",
                    xaxis = list(
                        showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE,
                        range = xwin(),                     # keep current window
                        fixedrange = FALSE
                    ),
                    yaxis = list(
                        showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE,
                        range = c(-1, y_max), fixedrange = TRUE     # <-- more vertical headroom
                    ),
                    margin = list(t = 30)                           # small top margin
                ) %>%
                plotly::config(
                    displaylogo = FALSE,
                    scrollZoom  = FALSE,                  # touchpad/wheel does nothing
                    modeBarButtonsToRemove = c("zoom2d", "lasso2d", "select2d")
                )

            # Connector line
            line_x <- c(
                protein_domains$start_adj[1] - subdomain_gap,
                as.vector(rbind(protein_domains$start_adj, protein_domains$end_adj)),
                protein_domains$end_adj[num_protein_domains] + subdomain_gap
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
            for (i in seq_len(nrow(protein_domains))) {
                p <- plotly::add_trace(
                    p,
                    x = c(
                        protein_domains$start_adj[i], protein_domains$end_adj[i],
                        protein_domains$end_adj[i], protein_domains$start_adj[i],
                        protein_domains$start_adj[i]
                    ),
                    y = c(-0.2, -0.2, 0.2, 0.2, -0.2),
                    type = "scatter", mode = "lines",
                    fill = "toself", fillcolor = protein_domains$color[i],
                    line = list(color = "black"),
                    text = protein_domains$info[i], hoverinfo = "text", hoveron = "fills",
                    showlegend = FALSE
                )
                p <- plotly::add_annotations(
                    p,
                    x = (protein_domains$start_adj[i] + protein_domains$end_adj[i]) / 2,
                    y = 0, 
                    text = protein_domains$label[i],
                    showarrow = FALSE, 
                    font = list(size = 14)
                )
            }

            # Boundary tick labels
            p <- plotly::add_annotations(
                p, 
                x = 1, 
                y = -0.3, 
                text = "1",
                showarrow = FALSE, 
                font = list(size = 14)
            )
            for (i in seq_len(nrow(protein_domains))) {
                p <- plotly::add_annotations(
                    p, 
                    x = protein_domains$end_adj[i] + subdomain_gap/2, 
                    y = -0.3,
                    text = paste0(protein_domain_positions[i + 1]),
                    showarrow = FALSE, 
                    font = list(size = 14)
                )
            }

            # Lollipops: only Top N within current window (with label collision avoidance)
            vv <- vis_v()
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
                    hovertemplate = "%{text}<br>AA pos: %{customdata}<extra></extra>",
                    text = vv$label, 
                    customdata = vv$aa_pos,
                    showlegend = FALSE
                )

                # ---- Lane-stacked labels (uses tunables above) ----
                out_id <- session$ns("proteinDiagram")
                w_px <- session$clientData[[paste0("output_", out_id, "_width")]]
                if (is.null(w_px) || w_px <= 0) w_px <- 800
                xr <- xwin()
                thr_x <- (xr[2] - xr[1]) * (min_label_sep_px / w_px)

                ord <- order(vv$x, vv$y_top)
                last_x <- rep(-Inf, max_label_lanes)
                lanes  <- rep(NA_integer_, nrow(vv))

                for (i in ord) {
                    for (L in seq_len(max_label_lanes)) {
                        if (vv$x[i] - last_x[L] >= thr_x) {
                            lanes[i] <- L - 1
                            last_x[L] <- vv$x[i]
                            break
                        }
                    }
                    # if not placed in any lane, label omitted
                }

                lab <- vv[!is.na(lanes), , drop = FALSE]
                if (nrow(lab) > 0) {
                    lab_y <- lab$y_top + (lanes[!is.na(lanes)] + lane_offset) * lane_lift
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
