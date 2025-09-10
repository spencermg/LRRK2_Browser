#!/usr/bin/env Rscript

# =========================================================================
# INITIALIZE VARIABLES
# =========================================================================

# Length of line segment between protein subdomains and exons
subdomain_gap <- 20

# Start positions for each protein subdomain and exon
protein_domain_positions <- c(1, 705, 800, 1335, 1511, 1674, 1879, 2142, 2498)
exon_positions <- c(
    1, 287, 373, 483, 572, 707, 842, 974, 1094, 1237, 1317,
    1424, 1554, 1679, 1792, 1937, 2077, 2206, 2377, 2636,
    2825, 2944, 3014, 3232, 3483, 3632, 3726, 3913, 4095, 
    4325, 4453, 4672, 4874, 4963, 5151, 5306, 5453, 5645,
    5792, 5893, 6084, 6245, 6416, 6517, 6712, 6906, 6979,
    7164, 7317, 7526, 7598, 9239
)

# Number of protein subdomains and exons
num_protein_domains <- length(protein_domain_positions) - 1
num_exons <- length(exon_positions) - 1

# Adjusted start positions of each protein subdomain and exon, accounting for line segments 
protein_domain_positions_adj <- protein_domain_positions + (0:(num_protein_domains)) * subdomain_gap
exon_positions_adj <- exon_positions + (0:(num_exons)) * subdomain_gap

# Metadata used for the boxes for each protein subdomain
protein_domain_names <- c("ARM", "ANK", "LRR", "ROC", "COR-A", "COR-B", "KIN", "WD-40")
exon_names <- paste("", 1:num_exons)

# Define colors to use for each subdomain and exon
protein_domain_colors <- c(
    "#AEC6CF", 
    "#CBAACB", 
    "#FFB7B2", 
    "#FF6961", 
    "#FFB347", 
    "#FDFD96", 
    "#77DD77", 
    "#B2FFFF"
)
exon_colors <- c(
    rep(protein_domain_colors[1], 15), # 15 exons corresponding to the first subdomain
    rep(protein_domain_colors[2], 3),  # 3 exons corresponding to the second subdomain
    rep(protein_domain_colors[3], 10), # 10 exons corresponding to the third subdomain
    rep(protein_domain_colors[4], 4),  # 4 exons corresponding to the fourth subdomain
    rep(protein_domain_colors[5], 2),  # 2 exons corresponding to the fifth subdomain
    rep(protein_domain_colors[6], 4),  # 4 exons corresponding to the sixth subdomain
    rep(protein_domain_colors[7], 6),  # 6 exons corresponding to the seventh subdomain
    rep(protein_domain_colors[8], 7)   # 7 exons corresponding to the eighth subdomain
)

# Add tooltips for each box
### TODO: Replace info tooltips with something meaningful
protein_domain_info <- rep("Protein subdomain", num_protein_domains)
exon_info <- paste("Exon", 1:num_exons)

# Create data tables with metadata for each subdomain and exon
protein_domains <- data.frame(
    start = protein_domain_positions_adj[-(num_protein_domains+1)],
    end   = protein_domain_positions_adj[-1],
    label = protein_domain_names,
    color = protein_domain_colors,
    info  = protein_domain_info
)
exons <- data.frame(
    start = exon_positions_adj[-(num_exons+1)],
    end   = exon_positions_adj[-1],
    label = exon_names,
    color = exon_colors,
    info  = exon_info
)


# =========================================================================
# LOAD GENOMIC VARIANT DATA TABLES
# =========================================================================

# Load data for each ancestry separately and also combined
df <- fread("lrrk2_data.tsv")
ancestry_tables <- split(df, df$Ancestry)
combined_table <- fread("lrrk2_grouped.tsv")
all_tables <- c(list(Combined = combined_table), ancestry_tables)

# Define function to process tables
clean_variant_table <- function(tbl) {
    # Find frequency of cases and controls
    tbl <- as.data.table(tbl)
    if (all(c("het_PD", "hom_PD", "total_PD") %in% names(tbl))) {
        tbl[, `PD Frequency` := (het_PD + hom_PD) / total_PD]
    }
    if (all(c("het_HC", "hom_HC", "total_HC") %in% names(tbl))) {
        tbl[, `Control Frequency` := (het_HC + hom_HC) / total_HC]
    }

    # Keep only selected columns
    keep_cols <- c(
        "variant",
        "PD Frequency",
        "Control Frequency",
        "gnomad41_genome_AF",
        "Func.refGene",
        "ExonicFunc.refGene",
        "CADD_phred",
        "eQTLGen_snp_id",
        "CLNSIG",
        "CLNDN",
        "exon",
        "cDNA",
        "prot_change",
        "domain",
        "Consurf_score",
        "Mean_pRAB10/RAB10",
        "SD",
        "Interpretation"
    )
    tbl <- tbl[, intersect(keep_cols, names(tbl)), with = FALSE]

    # Replace "." with NA so they appear blank in the table
    na_cols <- c("ExonicFunc.refGene", "CADD_phred", "eQTLGen_snp_id", "gnomad41_genome_AF", "CLNSIG", "CLNDN")
    for (col in intersect(na_cols, names(tbl))) {
        tbl[[col]][tbl[[col]] == "."] <- NA
    }

    # Convert numeric columns
    for (col in c("CADD_phred", "gnomad41_genome_AF")) {
        if (col %in% names(tbl)) {
            tbl[[col]] <- suppressWarnings(as.numeric(tbl[[col]]))
        }
    }

    # Clean up exon column
    if ("exon" %in% names(tbl)) {
        tbl$exon <- gsub("exon", "", tbl$exon)
        tbl$exon <- suppressWarnings(as.numeric(tbl$exon))
    }

    # Rename columns
    colnames(tbl) <- c(
        "Variant",
        "PD frequency",
        "Control frequency",
        "Gnomad allele frequency",
        "Region",
        "Functional consequence",
        "CADD",
        "rsID",
        "Clinical significance",
        "Clinical disease name",
        "Exon #",
        "cDNA change",
        "AA change",
        "Protein domain",
        "Conservation score",
        "Kinase activity (mean pRAB10/RAB10)",
        "Standard deviation",
        "Interpretation"
    )

    return(tbl)
}

# Process tables, keeping combined tables by default
all_tables_cleaned <- lapply(all_tables, clean_variant_table)
variantTable <- all_tables_cleaned$Combined

# Define exonic variant annotation summary table to display in the UI
variantTable.exonic <- variantTable[variantTable$Region == "exonic"]
variantTable.exonic.counts <- variantTable.exonic[, .N, by = `Functional consequence`]
colnames(variantTable.exonic.counts) <- c(
    "Exonic variant functional consequence",
    paste0("Count (Total: ", dim(variantTable.exonic)[1], ")")
)

# Define noncoding variant annotation summary table to display in the UI
variantTable.noncoding <- variantTable[variantTable$Region != "exonic"]
variantTable.noncoding.counts <- variantTable.noncoding[, .N, by = Region]
colnames(variantTable.noncoding.counts) <- c(
    "Noncoding variant category",
    paste0("Count (Total: ", dim(variantTable.noncoding)[1], ")")
)


# =========================================================================
# SERVER FUNCTION
# =========================================================================

server <- function(input, output) {
    # General metadata about LRRK2
    output$geneOverview <- renderUI({
        fluidRow(
            column(
                width = 3,
                align = "center",
                descriptionBlock(
                    header = tags$b("CHR"),
                    text = 12,
                    rightBorder = TRUE,
                    marginBottom = FALSE
                )
            ),
            column(
                width = 3,
                align = "center",
                descriptionBlock(
                    header = tags$b("BP (GrCh38)"),
                    text = "40224997 - 40369285",
                    rightBorder = TRUE,
                    marginBottom = FALSE
                )
            ),
            column(
                width = 3,
                align = "center",
                descriptionBlock(
                    header = tags$b("BP (T2T)"),
                    text = "40177355 - 40321422",
                    rightBorder = TRUE,
                    marginBottom = FALSE
                )
            ),
            column(
                width = 3,
                align = "center",
                descriptionBlock(
                    header = tags$b("BP (GrCh37)"),
                    text = "40618799 - 40763087",
                    rightBorder = TRUE,
                    marginBottom = FALSE
                )
            )
        )
    })


    # Counts of variants across functional annotation categories 
    output$annotationSummaryTable <- renderUI({
        total_count <- dim(variantTable.exonic)[1]
        tagList(
            div(
                # Left table for exonic variants
                column(
                    width = 6,
                    renderDT({
                        datatable(
                            variantTable.exonic.counts,
                            rownames = F,
                            selection = 'none',
                            options = list(
                                paging = F,
                                dom = 't'
                            )
                        ) %>% formatStyle(
                            columns = "Count",
                            valueColumns = "Exonic variant functional consequence",
                            target = 'cell',
                            color = "black"
                        ) %>% formatStyle(columns="Exonic variant functional consequence", backgroundColor = "white", color = "black")
                    })
                ),
                # Right table for noncoding variants
                column(
                    width = 6,
                    renderDT({
                        datatable(
                            variantTable.noncoding.counts,
                            rownames = F,
                            selection = 'none',
                            options = list(
                                paging = F,
                                dom = 't'
                            )
                        ) %>% formatStyle(
                            columns = "Count",
                            valueColumns = "Noncoding variant category",
                            target = 'cell',
                            color = "black"
                        ) %>% formatStyle(columns="Noncoding variant category", backgroundColor = "white", color = "black")
                    })
                )
            )
        )
    })


    # Protein domain diagram
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


    # cDNA diagram
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

    
    output$geneVartable <- renderUI(tagList(
        fluidRow(
            div(
                # Define dropdown to select ancestry shown in the table
                selectInput(
                    inputId = "vartable.dataset",
                    label   = "Choose ancestry:",
                    choices = names(all_tables_cleaned),
                    selected = "Combined",
                    width = "250px"
                ),

                # Define filtering options at the top of the table
                radioGroupButtons(
                    inputId = "vartable.filter",
                    label = div(
                        "Filter variants:",
                        title = "For more options, please use the \"Search\" box on the right!"
                    ),
                    choices = c(
                        "No filter",
                        "Nonsynonymous",
                        "Frameshift",
                        "Stop gain/loss"
                    ),
                    selected = "No filter",
                    status = "primary",
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon"),
                        no = icon("remove", lib = "glyphicon")
                    )
                ),

                renderDT({
                    # Define table and apply filters
                    dat <- all_tables_cleaned[[ req(input$vartable.dataset) ]]
                    if (!is.na(input$vartable.filter)) {
                        dat <- switch(
                            input$vartable.filter,
                            "No filter" = dat,
                            "Nonsynonymous" = dat[dat$`Functional consequence` %ni% c(".", "synonymous SNV")],
                            "Frameshift" = dat[dat$`Functional consequence` %in% c("frameshift deletion", "frameshift insertion", "frameshift block substitution")],
                            "Stop gain/loss" = dat[dat$`Functional consequence` %in% c("stopgain", "stoploss")]
                        )
                    }

                    # Populate table
                    datatable(
                        dat,
                        extensions = 'Buttons',
                        options = list(
                            dom = 'Blfrtip',
                            buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                            paging = TRUE,
                            pageLength = 25,
                            lengthChange = TRUE,
                            lengthMenu = list(c(10,25,50,100,500,-1), c("10","25","50","100","500","All")),
                            scrollX = T,
                            deferRender = TRUE
                        ),
                        rownames= FALSE,
                        escape = FALSE
                    ) %>% formatStyle(
                        columns=colnames(variantTable),
                        backgroundColor = "#FFFFFF",
                        color = "black"
                    )
                }),
                style = "margin: 12px 50px 50px 12px;"
            )
        )
    ))
}
