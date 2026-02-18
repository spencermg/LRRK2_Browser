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
    sidebar = dashboardSidebar(
        width = 350,
        sidebarMenu(
            id = "sidebar"
        ),
        tags$div(
            style = "padding: 20px; overflow-y: auto; height: calc(100vh - 100px);",
            tags$h3("About the LRRK2 Browser", style = "color: #0C8DC3; margin-top: 0;"),
            tags$hr(),
            
            tags$h4("Data Sources"),
            tags$p("This browser contains genotyping data from participants within the Global Parkinson's Genetics Program (GP2) across several modalities:"),
            tags$ul(
                tags$li(tags$strong("Imputed genotyping")),
                tags$li(tags$strong("Whole genome sequencing (WGS)")),
                tags$li(tags$strong("Raw genotyping")),
                tags$li(tags$strong("Clinical exome"))
            ),
            
            tags$h4("Sample Information"),
            tags$p("Total samples: #"),
            tags$ul(
                tags$li("Imputed genotyping: #"),
                tags$li("WGS: #"),
                tags$li("Raw genotyping: #"),
                tags$li("Clinical exome: #")
            ),
            tags$p("Parkinson's Disease cases: #"),
            tags$ul(
                tags$li("Imputed genotyping: #"),
                tags$li("WGS: #"),
                tags$li("Raw genotyping: #"),
                tags$li("Clinical exome: #")
            ),
            tags$p("Healthy controls: #"),
            tags$ul(
                tags$li("Imputed genotyping: #"),
                tags$li("WGS: #"),
                tags$li("Raw genotyping: #"),
                tags$li("Clinical exome: #")
            ),
            
            tags$h4("Annotations"),
            tags$p("Variants are annotated with ANNOVAR using the following libraries:"),
            tags$ul(
                tags$li("CADD deleteriousness scores"),
                tags$li("ClinVar pathogenicity"),
                tags$li("Population frequencies (gnomAD v4.1)")
            ),
            tags$p("As well as conservation scores and kinase activity measurements from: _______"),
            
            tags$h4("Citation"),
            tags$p("For more information, or if you use data from this browser, please refer to:"),
            tags$p(
                tags$em("Citation"),
                style = "background-color: #f5f5f5; padding: 10px; border-left: 3px solid #0C8DC3;"
            ),
            
            tags$h4("Contact"),
            tags$p("For questions: ", tags$a(href = "mailto:example@example.com", "example@example.com")),
            
            tags$p(
                style = "margin-top: 30px; font-size: 12px; color: #666;",
                "Version 1.0 | Last updated: ", format(Sys.Date(), "%B %Y")
            )
        )
    ),
    body = dashboardBody(
        # CSS to fix sidebar colors and hide text when collapsed
        tags$head(
            tags$style(HTML("
                /* Disable bounce/overscroll effect */
                html, body {
                    overscroll-behavior-y: none !important;
                }
                
                /* When sidebar is open, shift navbar by 350px */
                .main-header .navbar {
                    margin-left: 350px !important;
                    transition: margin-left 0.3s ease !important;
                }
                
                /* When sidebar is collapsed, no margin */
                .sidebar-collapse .main-header .navbar {
                    margin-left: 0px !important;
                }
            
                /* Style the toggle button */
                .sidebar-toggle {
                    color: white !important;
                    padding-right: 15px !important;
                }
                
                /* Add 'About the Data' label after the three lines */
                .sidebar-toggle::after {
                    content: 'About the Data';
                    margin-left: 10px;
                    font-size: 14px;              /* Adjust size */
                    font-weight: bold;             /* or normal, 600, 700, etc. */
                    font-family: 'Arial', sans-serif;  /* Change font family */
                    font-style: normal;            /* or italic */
                    letter-spacing: 0.5px;         /* Add spacing between letters */
                    text-transform: none;          /* or uppercase, capitalize, lowercase */
                    vertical-align: middle;
                }
                
                /* Make entire sidebar background match header blue */
                .main-sidebar {
                    background-color: #367fa9 !important;
                    position: fixed !important;
                    top: 0 !important;
                    left: 0 !important;
                    height: 100vh !important;
                    overflow-y: auto !important;
                }
                
                /* The content area inside sidebar should be white/light */
                .main-sidebar .sidebar {
                    background-color: #ffffff !important;
                }
                
                /* Make all text dark/visible */
                .sidebar h3, .sidebar h4, .sidebar p, .sidebar li, .sidebar a {
                    color: #333333 !important;
                }
                
                /* Hide text content when sidebar is collapsed */
                .sidebar-collapse .main-sidebar .sidebar {
                    display: none !important;
                }
                
                /* Fix the header portion of sidebar to be blue */
                .main-sidebar .sidebar-menu {
                    background-color: #367fa9 !important;
                }
            "))
        ),
        
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
    controlbar = NULL
)
