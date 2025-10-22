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
            print(variant_details)

            # Indicate columns to display in the table and convert frequencies to scientific notation
            variant_display <- variant_details[, c("Ancestry", "PD frequency", "Control frequency", "Gnomad allele frequency"), drop = FALSE]
            for (col in c("PD frequency", "Control frequency", "Gnomad allele frequency")) {
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
                            tags$p(tags$strong("Functional Consequence: "), variant_details$`Functional consequence`[1])
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
                            )
                        )
                    )
                )
            } else {
                tags$p("No details available for this variant.")
            }
            
            # Show popup with variant details
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