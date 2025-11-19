#!/usr/bin/env Rscript

# =========================================================================
# UI FUNCTION
# =========================================================================

variantDetailUI <- function(id) {
    ns <- NS(id)
}


# =========================================================================
# SERVER FUNCTION
# =========================================================================

variantDetailServer <- function(id, variant_data, all_tables_cleaned) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Observe when a variant is clicked
        observeEvent(variant_data(), {
            req(variant_data())
            
            variant_id <- variant_data()$variant_id

            # Create one data table with data from each ancestry and combined
            variant_details <- list()
            for (anc in c(setdiff(names(all_tables_cleaned), "Combined"), "Combined")) {
                dat <- all_tables_cleaned[[anc]]
                dat <- dat[`Variant (GrCh38)` == variant_id]
                if (nrow(dat) > 0) {
                    dat[, Ancestry := anc]
                    setcolorder(dat, c("Ancestry", setdiff(names(dat), "Ancestry")))
                    variant_details[[anc]] <- dat
                } else {
                    placeholder <- copy(all_tables_cleaned$Combined[`Variant (GrCh38)` == variant_id])
                    placeholder$Ancestry <- anc
                    placeholder$`PD frequency` <- 0
                    placeholder$`Control frequency` <- 0
                    variant_details[[anc]] <- placeholder
                }
            }
            variant_details <- rbindlist(variant_details, fill = TRUE)

            # Indicate columns to display in the table and convert frequencies to scientific notation
            variant_display <- variant_details[, c("Ancestry", "PD frequency", "Control frequency", "gnomAD allele frequency"), drop = FALSE]
            for (col in c("PD frequency", "Control frequency", "gnomAD allele frequency")) {
                variant_display[[col]] <- ifelse(
                    is.na(variant_display[[col]]),
                    "N/A",
                    format(variant_display[[col]], scientific = TRUE, digits = 5)
                )
                variant_display[[col]] <- gsub("^0e[\\+\\-]0+$", "0", variant_display[[col]])
            }

            # Build the popup content
            popup_content <- if (nrow(variant_details) > 0) {
                tagList(
                    # Variant Annotations
                    tags$h3(
                        "Overview", 
                        style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                    ),
                    tags$div(
                        style = "margin-bottom: 20px;",
                        if ("rsID" %in% colnames(variant_details) && !is.na(variant_details$rsID[1])) {
                            tags$p(tags$strong("rsID: "), variant_details$rsID[1])
                        },
                        if ("Functional consequence" %in% colnames(variant_details)) {
                            tags$p(
                                tags$strong("Functional Consequence: "), 
                                ifelse(
                                    !is.na(variant_details$`Functional consequence`[1]), 
                                    variant_details$`Functional consequence`[1], 
                                    "N/A"
                                )
                            )
                        },
                        if ("Region" %in% colnames(variant_details)) {
                            tags$p(tags$strong("Region: "), variant_details$Region[1])
                        },
                        if ("CADD" %in% colnames(variant_details)) {
                            tags$p(
                                tags$strong("CADD Score: "), 
                                ifelse(
                                    !is.na(variant_details$CADD[1]), 
                                    round(variant_details$CADD[1], 2), 
                                    "N/A"
                                )
                            )
                        },
                        if ("Conservation score" %in% colnames(variant_details)) {
                            tags$p(
                                tags$strong("Conservation Score: "), 
                                ifelse(
                                    !is.na(variant_details$`Conservation score`[1]), 
                                    round(variant_details$`Conservation score`[1], 2), 
                                    "N/A"
                                )
                            )
                        },
                        if ("Kinase activity (mean pRAB10/RAB10)" %in% colnames(variant_details)) {
                            tags$p(
                                tags$strong("Kinase Activity (mean pRAB10/RAB10): "), 
                                ifelse(
                                    !is.na(variant_details$`Kinase activity (mean pRAB10/RAB10)`[1]),
                                    round(variant_details$`Kinase activity (mean pRAB10/RAB10)`[1], 3), 
                                    "N/A"
                                )
                            )
                        },
                        if ("Clinical significance" %in% colnames(variant_details)) {
                            tags$p(
                                tags$strong("Clinical Significance: "), 
                                ifelse(
                                    !is.na(variant_details$`Clinical significance`[1]), 
                                    variant_details$`Clinical significance`[1], 
                                    "N/A"
                                )
                            )
                        },
                        if ("Clinical disease name" %in% colnames(variant_details)) {
                            tags$p(
                                tags$strong("Clinical Disease Name: "), 
                                ifelse(
                                    !is.na(variant_details$`Clinical disease name`[1]), 
                                    variant_details$`Clinical disease name`[1], 
                                    "N/A"
                                )
                            )
                        },
                        if ("cDNA change" %in% colnames(variant_details) && variant_details$Region[1] == "exonic") {
                            tags$p(tags$strong("cDNA Change: "), variant_details$`cDNA change`[1])
                        },
                        if ("AA change" %in% colnames(variant_details) && variant_details$Region[1] == "exonic") {
                            tags$p(tags$strong("Amino Acid Change: "), variant_details$`AA change`[1])
                        },
                        if ("Exon #" %in% colnames(variant_details) && variant_details$Region[1] == "exonic") {
                            tags$p(tags$strong("Exon #: "), variant_details$`Exon #`[1])
                        },
                        if ("Protein domain" %in% colnames(variant_details) && variant_details$Region[1] == "exonic") {
                            tags$p(tags$strong("Protein Domain: "), variant_details$`Protein domain`[1])
                        }
                    ),
                    
                    # Allele frequencies across each ancestry
                    tags$h3(
                        "Variant Frequencies by Ancestry", 
                        style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                    ),
                    tags$div(
                        style = "overflow-x: auto; margin-top: 10px;",
                        tags$table(
                            class = "table table-striped table-hover",
                            style = "width: 100%; font-size: 12px;",
                            tags$thead(
                                tags$tr(
                                    lapply(colnames(variant_display), function(col) {
                                        tags$th(col, style = "background-color: #0C8DC3; color: white; padding: 8px;")
                                    })
                                )
                            ),
                            tags$tbody(
                                lapply(1:nrow(variant_display), function(i) {
                                    tags$tr(
                                        lapply(variant_display[i, ], function(val) {
                                            tags$td(
                                                ifelse(is.na(val), "N/A", as.character(val)),
                                                style = "padding: 8px;"
                                            )
                                        })
                                    )
                                })
                            ),
                            tags$p(
                                "*Some ancestries do not have reported allele frequencies through gnomAD, and are thus indicated as 'N/A'.",
                                style = "font-size: 11px; font-style: italic; color: #555; margin-top: 5px; text-align: center;"
                            )
                        )
                    ),

                    # Family History Pie Charts
                    tags$h3(
                        "Family History of Parkinson's Disease Among Carriers",
                        style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                    ),
                    fluidRow(
                        column(6, plotlyOutput(ns("pd_pie"), height = "300px")),
                        column(6, plotlyOutput(ns("control_pie"), height = "300px"))
                    ),

                    # AAO Histogram
                    tags$h3(
                        "Age at Onset (AAO) Distribution Among Carriers",
                        style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                    ),
                    plotlyOutput(ns("aao_hist"), height = "300px")
                )
            } else {
                tags$p("No details available for this variant.")
            }

            # Prepare family history data for pie charts
            pd_fh <- colSums(variant_details[, c("No family history (PD)", "Family history (PD)", "Unknown family history (PD)"), with = FALSE], na.rm = TRUE)
            control_fh <- colSums(variant_details[, c("No family history (Control)", "Family history (Control)", "Unknown family history (Control)"), with = FALSE], na.rm = TRUE)
            fh_labels <- c("No", "Yes", "Unknown")
            fh_colors <- c("#8C4E9F", "#4FC190", "#D3D3D3")

            # Keep zeros as NA to preserve order but hide from display
            pd_fh_display <- ifelse(pd_fh == 0, NA, pd_fh)
            control_fh_display <- ifelse(control_fh == 0, NA, control_fh)

            # Render PD pie chart
            output$pd_pie <- renderPlotly({
                # If there are no carriers, show a message instead of a chart
                if (variant_details[Ancestry == "Combined", `PD frequency`][1] == 0) {
                    return(plotly_empty(type = "scatter", mode = "text") %>%
                        layout(
                            annotations = list(
                                text = "No PD carrier data \navailable for this variant",
                                x = 0.5, 
                                y = 0.5, 
                                showarrow = FALSE,
                                font = list(size = 14, color = "#888")
                            )
                        )
                    )
                }
                plot_ly(
                    labels = fh_labels,
                    values = pd_fh_display,
                    type = "pie",
                    sort = FALSE,
                    direction = "counterclockwise",
                    marker = list(colors = fh_colors),
                    name = "PD",
                    hoverinfo = "label+percent"
                ) %>% layout(
                    title = "PD Cases",
                    showlegend = TRUE
                )
            })

            # Render Control pie chart
            output$control_pie <- renderPlotly({
                # If there are no carriers, show a message instead of a chart
                if (variant_details[Ancestry == "Combined", `Control frequency`][1] == 0) {
                    return(plotly_empty(type = "scatter", mode = "text") %>%
                        layout(
                            annotations = list(
                                text = "No healthy control carrier data \navailable for this variant",
                                x = 0.5, 
                                y = 0.5, 
                                showarrow = FALSE,
                                font = list(size = 14, color = "#888")
                            )
                        )
                    )
                }
                plot_ly(
                    labels = fh_labels,
                    values = control_fh_display,
                    type = "pie",
                    sort = FALSE,
                    direction = "counterclockwise",
                    marker = list(colors = fh_colors),
                    name = "Control",
                    hoverinfo = "label+percent"
                ) %>% layout(
                    title = "Healthy Controls",
                    showlegend = TRUE
                )
            })


            # Render AAO histogram
            output$aao_hist <- renderPlotly({
                # Extract the counts for each AAO bin
                aao_cols <- c(
                    "AAO (1-10)", "AAO (11-20)", "AAO (21-30)", "AAO (31-40)",
                    "AAO (41-50)", "AAO (51-60)", "AAO (61-70)", "AAO (71-80)",
                    "AAO (81-90)", "AAO (91-100)"
                )
                aao_counts <- colSums(variant_details[, ..aao_cols], na.rm = TRUE)

                # If there are no carriers, show a message instead of a chart
                if (variant_details[Ancestry == "Combined", `PD frequency`][1] == 0) {
                    return(plotly_empty(type = "scatter", mode = "text") %>%
                        layout(
                            annotations = list(
                                text = "No PD carrier data \navailable for this variant",
                                x = 0.5, 
                                y = 0.5, 
                                showarrow = FALSE,
                                font = list(size = 14, color = "#888")
                            )
                        )
                    )
                }

                # Prepare labels and data
                aao_labels <- gsub("AAO \\(|\\)", "", aao_cols)
                df_aao <- data.frame(
                    Range = factor(aao_labels, levels = aao_labels),
                    Count = as.numeric(aao_counts)
                )

                # Get min/median/max across all ancestries
                dat_combined <- variant_details[Ancestry == "Combined"]
                min_aao <- suppressWarnings(as.numeric(dat_combined$`Minimum AAO`[1]))
                med_aao <- suppressWarnings(as.numeric(dat_combined$`Median AAO`[1]))
                max_aao <- suppressWarnings(as.numeric(dat_combined$`Maximum AAO`[1]))

                # Build the histogram
                plot_ly(
                    data = df_aao,
                    x = ~Range,
                    y = ~Count,
                    type = "bar",
                    marker = list(color = "#0C8DC3"),
                    hovertemplate = "Count: %{y}<extra></extra>"
                ) %>%
                    layout(
                        xaxis = list(
                            title = paste0(
                                "Age at Onset Range (years)<br>",
                                "Min: ", ifelse(is.na(min_aao), "N/A", min_aao),
                                ", Median: ", ifelse(is.na(med_aao), "N/A", med_aao),
                                ", Max: ", ifelse(is.na(max_aao), "N/A", max_aao)
                            )
                        ),
                        yaxis = list(title = "Count"),
                        bargap = 0.2
                    )
            })

            # Show popup
            showModal(modalDialog(
                title = tags$div(
                    style = "font-size: 24px; font-weight: bold; text-align: center;",
                    paste("Variant Details:", variant_id)
                ),
                size = "l",
                easyClose = TRUE,
                footer = modalButton("Close"),
                fluidRow(column(12, tags$div(style = "padding: 20px;", popup_content)))
            ))
        })
    })
}
