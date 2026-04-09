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

variantDetailServer <- function(id, all_tables_cleaned, variant_bus) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        variant_data <- reactive({
            variant_bus$variant()
        })

        # Observer to show charts when warning is acknowledged
        observeEvent(input$acknowledge_warning, {
            shinyjs::hide(id = "charts_warning_section")
            shinyjs::show(id = "charts_section")
        })
        
        # Observe when a variant is clicked
        observeEvent(variant_bus$variant(), {
            variant <- variant_data()
            req(variant)
            variant_id <- variant$variant_id

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
                    placeholder$`PD allele frequency (Imputed)` <- 0
                    placeholder$`Control allele frequency (Imputed)` <- 0
                    placeholder$`PD allele frequency (WGS)` <- 0
                    placeholder$`Control allele frequency (WGS)` <- 0
                    placeholder$`PD allele frequency (Raw genotyping)` <- 0
                    placeholder$`Control allele frequency (Raw genotyping)` <- 0
                    
                    # Only set exome column if it exists (for Combined ancestry)
                    if ("PD allele frequency (CES)" %in% names(placeholder)) {
                        placeholder$`PD allele frequency (CES)` <- 0
                    }
                    
                    variant_details[[anc]] <- placeholder
                }
            }
            variant_details <- rbindlist(variant_details, fill = TRUE)

            # Indicate columns to display in the table and convert frequencies to scientific notation
            display_cols <- c(
                "Ancestry", 
                "PD allele frequency (WGS)", 
                "Control allele frequency (WGS)", 
                "PD allele frequency (Imputed)", 
                "Control allele frequency (Imputed)", 
                "PD allele frequency (Raw genotyping)", 
                "Control allele frequency (Raw genotyping)"
            )

            # Only include exome column if it exists
            if ("PD allele frequency (CES)" %in% names(variant_details)) {
                display_cols <- c(display_cols, "PD allele frequency (CES)")
            }

            display_cols <- c(display_cols, "gnomAD allele frequency")

            # Only select columns that actually exist
            display_cols <- display_cols[display_cols %in% names(variant_details)]
            variant_display <- variant_details[, ..display_cols]

            variant_display <- variant_display[, .SD, .SDcols = !sapply(variant_display, function(col) all(is.na(col)))]

            # Convert frequencies to scientific notation
            freq_cols <- c(
                "PD allele frequency (WGS)", 
                "Control allele frequency (WGS)", 
                "PD allele frequency (Imputed)", 
                "Control allele frequency (Imputed)", 
                "PD allele frequency (Raw genotyping)", 
                "Control allele frequency (Raw genotyping)", 
                "PD allele frequency (CES)", 
                "gnomAD allele frequency"
            )

            for (col in intersect(freq_cols, names(variant_display))) {
                variant_display[[col]] <- ifelse(
                    is.na(variant_display[[col]]),
                    "N/A",
                    ifelse(
                        variant_display[[col]] == 0,
                        "0",
                        ifelse(
                            abs(variant_display[[col]]) < 0.001,
                            format(variant_display[[col]], scientific = TRUE, digits = 3),
                            format(variant_display[[col]], scientific = FALSE, digits = 5)
                        )
                    )
                )
                variant_display[[col]] <- gsub("^0e[\\+\\-]0+$", "0", variant_display[[col]])
            }

            # Build the popup content
            popup_content <- if (nrow(variant_details) > 0) {
                is_pathogenic <- variant_details$Pathogenic[1]
                show_warning <- is.na(is_pathogenic) || is_pathogenic == 0
                
                tagList(
                    # Variant Annotations
                    tags$h3(
                        "Overview", 
                        style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                    ),
                    tags$div(
                        style = "margin-bottom: 20px;",
                        tags$p(
                            style = "text-align: center;",
                            tags$strong(
                                ifelse(
                                    is_pathogenic == 1, 
                                    "This variant IS classified as pathogenic according to GP2 criteria", 
                                    "This variant IS NOT classified as pathogenic according to GP2 criteria"
                                )
                            )
                        ),
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
                            ),
                            tags$p(
                                "**Individual ancestry assignments are not present for clinical exome data, and are thus indicated as 'N/A'.",
                                style = "font-size: 11px; font-style: italic; color: #555; margin-top: 5px; text-align: center;"
                            ),
                            tags$p(
                                "***Some columns may be dropped from the data table if this variant is not present in the data for a given data modality.",
                                style = "font-size: 11px; font-style: italic; color: #555; margin-top: 5px; text-align: center;"
                            )
                        )
                    ),

                    # Conditionally show warning OR charts
                    if (show_warning) {
                        tags$div(
                            id = ns("charts_warning_section"),
                            # Warning banner
                            tags$div(
                                style = paste0(
                                    "background-color: #fff3cd; ",
                                    "border: 3px solid #ffc107; ",
                                    "border-radius: 10px; ",
                                    "padding: 30px; ",
                                    "margin: 20px 0; ",
                                    "text-align: center;"
                                ),
                                tags$p(
                                    style = "font-size: 15px; color: #856404; margin-bottom: 15px;",
                                    tags$strong("⚠️ This variant IS NOT classified as pathogenic according to GP2 criteria.")
                                ),
                                tags$p(
                                    style = "font-size: 14px; color: #666; margin-bottom: 20px;",
                                    "Family history and age-at-onset data for non-pathogenic variants should be ",
                                    "interpreted with caution, as these may represent benign variation or variants ",
                                    "of uncertain significance."
                                ),
                                actionButton(
                                    ns("acknowledge_warning"),
                                    "I Understand - Show Charts",
                                    class = "btn-warning",
                                    style = "font-size: 16px; padding: 10px 30px;"
                                )
                            ),
                            tags$div(
                                style = "text-align: center; padding: 40px; color: #999;",
                                tags$p("Click 'I Understand' above to view these charts")
                            )
                        )
                    } else {
                        # For pathogenic variants OR already acknowledged, show charts directly
                        tagList(
                            tags$h3(
                                "Family History of Parkinson's Disease Among Carriers",
                                style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                            ),
                            fluidRow(
                                column(6, plotlyOutput(ns("pd_pie"), height = "300px")),
                                column(6, plotlyOutput(ns("control_pie"), height = "300px"))
                            ),
                            tags$h3(
                                "Age at Onset (AAO) Distribution Among PD-Affected Carriers",
                                style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                            ),
                            plotlyOutput(ns("aao_hist"), height = "300px"),
                            tags$h3(
                                "Age Distribution Among Healthy Control Carriers",
                                style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center;"
                            ),
                            plotlyOutput(ns("age_hist"), height = "300px")
                        )
                    },

                    # Hidden div for charts (shown after acknowledgment)
                    if (show_warning) {
                        tags$div(
                            id = ns("charts_section"),
                            style = "display: none;",
                            tags$h3(
                                "Family History of Parkinson's Disease Among Carriers",
                                style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center; margin-top: 30px;"
                            ),
                            fluidRow(
                                column(6, plotlyOutput(ns("pd_pie"), height = "300px")),
                                column(6, plotlyOutput(ns("control_pie"), height = "300px"))
                            ),
                            tags$h3(
                                "Age at Onset (AAO) Distribution Among PD-Affected Carriers",
                                style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center; margin-top: 30px;"
                            ),
                            plotlyOutput(ns("aao_hist"), height = "300px"),
                            tags$h3(
                                "Age Distribution Among Healthy Control Carriers",
                                style = "color: #0C8DC3; border-bottom: 2px solid #0C8DC3; padding-bottom: 5px; text-align: center; margin-top: 30px;"
                            ),
                            plotlyOutput(ns("age_hist"), height = "300px")
                        )
                    }
                )
            } else {
                tags$p("No details available for this variant.")
            }

            # Prepare family history data for pie charts
            pd_fh <- variant_details[Ancestry == "Combined", c("No family history (PD)", "Family history (PD)", "Unknown family history (PD)"), with = FALSE]
            control_fh <- variant_details[Ancestry == "Combined", c("No family history (Control)", "Family history (Control)", "Unknown family history (Control)"), with = FALSE]
            fh_labels <- c("No", "Yes", "Unknown")
            fh_colors <- c("#8C4E9F", "#34A270", "#D3D3D3")

            # Keep zeros as NA to preserve order but hide from display
            pd_fh_display <- ifelse(pd_fh == 0, NA, pd_fh)
            control_fh_display <- ifelse(control_fh == 0, NA, control_fh)

            # Render PD pie chart
            output$pd_pie <- renderPlotly({
                # If there are no carriers, show a message instead of a chart
                if (all(is.na(pd_fh_display))) {
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
                    hoverinfo = "label+value"
                ) %>% layout(
                    title = "PD Cases",
                    showlegend = TRUE
                )
            })

            # Render Control pie chart
            output$control_pie <- renderPlotly({
                # If there are no carriers, show a message instead of a chart
                if (all(is.na(control_fh_display))) {
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
                    hoverinfo = "label+value"
                ) %>% layout(
                    title = "Healthy Controls",
                    showlegend = TRUE
                )
            })

            # Render AAO histogram for PD-affected carriers
            output$aao_hist <- renderPlotly({
                # Extract the counts for each AAO bin
                aao_cols <- c(
                    "AAO (11-20)", "AAO (21-30)", "AAO (31-40)", 
                    "AAO (41-50)", "AAO (51-60)", "AAO (61-70)", 
                    "AAO (71-80)", "AAO (81-90)", "AAO (91-100)"
                )
                aao_counts <- as.numeric(variant_details[Ancestry == "Combined", ..aao_cols])

                # If there are no carriers, show a message instead of a chart
                if (sum(aao_counts) == 0) {
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
                min_aao <- floor(suppressWarnings(as.numeric(dat_combined$`Minimum AAO`[1])))
                med_aao <- floor(suppressWarnings(as.numeric(dat_combined$`Median AAO`[1])))
                max_aao <- floor(suppressWarnings(as.numeric(dat_combined$`Maximum AAO`[1])))

                # Build the histogram
                plot_ly(
                    data = df_aao,
                    x = ~Range,
                    y = ~Count,
                    type = "bar",
                    marker = list(color = "#0C8DC3"),
                    hovertemplate = "Count: %{y}<extra></extra>"
                ) %>% layout(
                    xaxis = list(
                        title = paste0(
                            "Age at Onset Range (years)<br>",
                            "Min: ", ifelse(is.na(min_aao), "N/A", min_aao),
                            ", Median: ", ifelse(is.na(med_aao), "N/A", round(med_aao, 2)),
                            ", Max: ", ifelse(is.na(max_aao), "N/A", max_aao), 
                            "<br>",
                            "*AAO data available for ", sum(aao_counts), " out of ", sum(unlist(pd_fh_display), na.rm = TRUE), " PD-affected carriers"
                        )
                    ),
                    yaxis = list(title = "Count"),
                    bargap = 0.2
                )
            })

            # Render age histogram for healthy controls
            output$age_hist <- renderPlotly({
                # Extract the counts for each age bin
                age_cols <- c(
                    "Age (0-10)", "Age (11-20)", "Age (21-30)", "Age (31-40)", 
                    "Age (41-50)", "Age (51-60)", "Age (61-70)", "Age (71-80)", 
                    "Age (81-90)", "Age (91-100)", "Age (101-110)", "Age (111-120)"
                )
                age_counts <- as.numeric(variant_details[Ancestry == "Combined", ..age_cols])

                # If there are no carriers, show a message instead of a chart
                if (sum(age_counts) == 0) {
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

                # Prepare labels and data
                age_labels <- gsub("Age \\(|\\)", "", age_cols)
                df_age <- data.frame(
                    Range = factor(age_labels, levels = age_labels),
                    Count = as.numeric(age_counts)
                )

                # Get min/median/max across all ancestries
                dat_combined <- variant_details[Ancestry == "Combined"]
                min_age <- floor(suppressWarnings(as.numeric(dat_combined$`Minimum Age`[1])))
                med_age <- floor(suppressWarnings(as.numeric(dat_combined$`Median Age`[1])))
                max_age <- floor(suppressWarnings(as.numeric(dat_combined$`Maximum Age`[1])))

                # Build the histogram
                plot_ly(
                    data = df_age,
                    x = ~Range,
                    y = ~Count,
                    type = "bar",
                    marker = list(color = "#0C8DC3"),
                    hovertemplate = "Count: %{y}<extra></extra>"
                ) %>% layout(
                    xaxis = list(
                        title = paste0(
                            "Age Range (years)<br>",
                            "Min: ", ifelse(is.na(min_age), "N/A", min_age),
                            ", Median: ", ifelse(is.na(med_age), "N/A", round(med_age, 2)),
                            ", Max: ", ifelse(is.na(max_age), "N/A", max_age), 
                            "<br>",
                            "*Age data available for ", sum(age_counts), " out of ", sum(unlist(control_fh_display), na.rm = TRUE), " healthy control carriers"
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
