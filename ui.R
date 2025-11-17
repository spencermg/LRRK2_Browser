#!/usr/bin/env Rscript

# =========================================================================
# INITIALIZE VARIABLES
# =========================================================================

metadata_box_color         <- "#0C8DC3" # GP2 dark blue
variant_function_box_color <- "#D8EAF8" # GP2 light blue
protein_diagram_box_color  <- "#EAF4FF" # GP2 light grey
cdna_diagram_box_color     <- "#8C4E9F" # GP2 purple
variant_table_color        <- "#CF2FB3" # GP2 pink
box_outline_color          <- "#D3D3D3" # Light grey
box_title_font_style       <- "font-weight: bold; font-size: 28px;"


# =========================================================================
# UI FUNCTION
# =========================================================================

ui <- dashboardPage(
    title = "LRRK2 Browser",
    header = dashboardHeader(),
    sidebar = dashboardSidebar(collapsed = TRUE),
    body = dashboardBody(
        # Title banner
        h2("LRRK2 Browser", style = "text-align:center; font-weight:bold;"),

        # Overview numbers for LRRK2
        box(
            title = tags$div(
                "LRRK2 Metadata",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", metadata_box_color, "; padding: 10px;"),
            tags$div(
                "Leucine Rich Repeat Kinase 2",
                style = "text-align: center; font-weight: bold; font-size: 20px; margin-bottom: 10px;"
            ),
            geneOverviewUI("gene_overview")
        ),

        # Links to external resources
        box(
            title = tags$div(
                "Other resources",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", metadata_box_color, "; padding: 10px;"),
            otherResourcesUI("other_resources")
        ),

        # Summary table for different variant functional annotation categories
        box(
            title = tags$div(
                "Variant function summary",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", variant_function_box_color, "; padding: 10px;"),
            annotationSummaryTableUI("annotation_summary_table")
        ),

        # Protein subdomain diagram
        box(
            title = tags$div(
                "Protein Diagram",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", protein_diagram_box_color, "; padding: 10px;"),
            diagramUI("protein_diagram")
        ),

        # cDNA diagram
        box(
            title = tags$div(
                "cDNA Diagram",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", cdna_diagram_box_color, "; padding: 10px;"),
            diagramUI("cdna_diagram")
        ),

        # cDNA diagram
        box(
            title = tags$div(
                "Kinase Activity by Variant",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", cdna_diagram_box_color, "; padding: 10px;"),
            barChartUI("bar_chart")
        ),

        # Main variant table
        box(
            title = tags$div(
                "Variant Table",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", variant_table_color, "; padding: 10px;"),
            geneVarTableUI("gene_var_table")
        )
    ),
    controlbar = dashboardControlbar()
)
