### EDGER EXPLORER

libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/4.1"
if (dir.exists(libDir)) .libPaths(libDir)

library(scales)
library(shiny)
library(shinycssloaders)
library(tidyverse)
library(DT)
source("../shiny_func.R")

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

data <- read_rds("../data_edger.rds")
contrasts <- unique(data$de$contrast)

gene2name <- set_names(data$genes$gene_name, data$genes$gene_id)

max_points <- 500

all_genes <- data$genes$gene_id %>% unique()
all_terms <- sh_prepare_for_enrichment(data, c("go", "re", "kg"))

#######################################################################

ui <- shinyUI(fluidPage(
  
  tags$style(css),
  
  titlePanel("Tfe1/Tfe2 toxin expression in yeast: differential expression"),

  fluidRow(
    column(12,
      fluidRow(
        column(4,
          conditionalPanel(
            condition = "input.plotType != 'Tfe correlation'",
            radioButtons("contrastSel", "Contrast:", choices = contrasts, inline = TRUE) 
          ),
          radioButtons("plotType", "Plot type:", choices = c("Volcano", "MA", "Tfe correlation"), inline = TRUE),
          plotOutput("mainPlot", height = "480px", width = "100%", brush = "plot_brush", hover = "plot_hover")
        ),
        column(3,
          radioButtons("intensityScale", "Intesity scale:", choices = c("lin" = "", "log" = "log"), inline = TRUE),
          plotOutput("genePlot", height = "400px",width = "100%")
        ),
        column(5,
          p("Gene list"),
          div(style = 'height: 200px; overflow-y: scroll', tableOutput("geneInfo")),
          br(),
          radioButtons("enrichment", "Enrichment:", choices = c("GO" = "go", "Reactome" = "re", "KEGG" = "kg"), inline = TRUE),
          div(style = 'height: 400px; overflow-y: scroll', tableOutput("Enrichment") %>% withSpinner(color = "#0dc5c1", type = 5, size = 0.5)),
        )
      ),
      fluidRow(
        DT::dataTableOutput("allGeneTable")
      )
    )
  )
)
)


########################################################################################


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    stopApp()
  })
  
  getXYData <- function() {
    if (input$plotType == "Tfe correlation") {
      xy_data <- data$tfe %>% 
        mutate(x = atanh(corTfe1), y = atanh(corTfe2), FDR = 0)
    } else if (input$plotType == "Volcano") {
      xy_data <- data$de %>% 
        filter(contrast == input$contrastSel) %>% 
        mutate(x = logFC, y = -log10(PValue))
    } else if (input$plotType == "MA") {
      xy_data <- data$de %>% 
        filter(contrast == input$contrastSel) %>% 
        mutate(x = logCPM, y = logFC)
    }
    xy_data
  }
  
  selectGene <- function(max_hover = 1) {
    xy_data <- getXYData()
    sel <- NULL
    tab_idx <- as.numeric(input$allGeneTable_rows_selected)
    if (!is.null(input$plot_brush)) {
      brushed <- brushedPoints(xy_data, input$plot_brush)
      sel <- brushed$gene_id
    } else if (!is.null(input$plot_hover)) {
      near <- nearPoints(xy_data, input$plot_hover, threshold = 20, maxpoints = max_hover)
      sel <- near$gene_id
    } else if (length(tab_idx) > 0) {
      sel <- xy_data[tab_idx, ] %>% pull(gene_id)
    }
    return(sel)
  }
  
  output$geneInfo <- renderTable({
    xy_data <- getXYData()
    sel <- selectGene()
    df <- NULL
    if (!is.null(sel) && length(sel) >= 1 && length(sel) <= max_points) {
      df <- xy_data %>%
        filter(gene_id %in% sel) %>% 
        arrange(gene_name)
      if (input$plotType == "Gradient") {
        df <- df %>% select(gene_name, gene_biotype, description)
      } else {
        df <- df %>% select(gene_name, gene_biotype, description, FDR)
      }
    } else if (length(sel) > max_points) {
      df <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
    }
    df
  })

  enrichmentTable <- function(terms) {
    xy_data <- getXYData()
    sel <- NULL
    fe <- NULL
    if (!is.null(input$plot_brush)) {
      brushed <- brushedPoints(xy_data, input$plot_brush)
      sel <- brushed$gene_id
      n <- length(sel)
      if (n > 0 && n <= max_points) {
        fe <- sh_functional_enrichment(all_genes, sel, terms, gene2name)
      } else if (n > 0) {
        fe <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  }
  
  output$Enrichment <- renderTable({
    d <- all_terms[[input$enrichment]]
    enrichmentTable(d)
  })
  
  
  output$genePlot <- renderPlot({
    sel <- selectGene()
    if (!is.null(sel) && length(sel) > 0 && length(sel) <= max_points) {
      sh_plot_genes(data$dat, data$metadata, sel, input$intensityScale, input$contrastSel)
    }
  })
  
  output$mainPlot <- renderPlot({
    xy_data <- getXYData()
    tab_idx <- as.numeric(input$allGeneTable_rows_selected)
    
    if (input$plotType == "Volcano") {
      g <- sh_plot_volcano(xy_data)
    } else if (input$plotType == "MA") {
      g <- sh_plot_ma(xy_data)
    } else if (input$plotType == "Tfe correlation") {
      g <- sh_plot_tfe(xy_data)
    }
    if (length(tab_idx) >= 1) {
      g <- g + geom_point(data = xy_data[tab_idx, ], colour = "red", size = 2)
    }
    g
  })

  output$allGeneTable <- DT::renderDataTable({
    if (input$plotType == "Tfe correlation") {
      d <- getXYData() %>%
        select(gene_name, description) %>% 
        mutate_if(is.numeric, ~signif(.x, 3))
    } else {
      d <- getXYData() %>%
        select(gene_name, gene_biotype, logFC, FDR, description) %>% 
        mutate_if(is.numeric, ~signif(.x, 3))
    }
    DT::datatable(d, class = 'cell-border strip hover', selection = "single", rownames = FALSE)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
