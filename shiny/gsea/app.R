### GENE SET ENRICHMENT EXPLORER

libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/4.1"
if(dir.exists(libDir)) .libPaths(libDir)

library(shiny)
library(tidyverse)
library(DT)
library(kableExtra)
source("../shiny_func.R")

select <- dplyr::select

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

data <- read_rds("../data_gsea.rds")
contrasts <- levels(data$res$contrast)


#######################################################################

ui <- shinyUI(fluidPage(
  
  #tags$style(css),
  
  titlePanel("Tfe1/Tfe2 toxin expression in yeast: gene set enrichment"),

  radioButtons("terms", "Functional terms", choices = c("GO" = "go", "Reactome" = "re", "KEGG" = "kg"), inline = TRUE),
  radioButtons("contrast", "Contrast:", choices = contrasts, inline = TRUE),

  fluidRow(
    column(6, dataTableOutput("termTable")),
    column(6, dataTableOutput("geneTable"))
  ),
  
  fluidRow(
    column(3, 
      br(),
      radioButtons("intensityScale", "Intesity scale:", choices = c("lin" = "", "log"="log"), inline = TRUE),
      plotOutput("genePlot", height = "400px",width = "50%")
    ),
    column(4, br(), tableOutput("geneInfo"))
  )
  
))


########################################################################################


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  getData <- reactive({
    fg <- data$fgs[[input$terms]] %>% 
      filter(contrast == input$contrast)
    fg
  })
  
  getTermData <- function() {
    fg <- getData()
    fg %>% 
      filter(padj < 0.05)
  }
  
  getGeneData <- function(trm) {
    fg <- getData()
    r <- fg %>% 
      filter(term == trm)
    df <- NULL
    if(nrow(r) > 0) {
      genes <- r$leading_edge[[1]]
      df <- data$res %>% 
        filter(gene_id %in% genes & contrast == input$contrast)
    }
  }
  
  selectTerm <- function() {
    sel <- NULL
    term_data <- getTermData()
    tab_idx <- as.numeric(input$termTable_rows_selected)
    if(length(tab_idx) == 1) {
      sel <- term_data[tab_idx, ] %>% pull(term)
    }
    return(sel)
  }
  
  selectGene <- function() {
    sel <- NULL
    term <- selectTerm()
    if(!is.null(term)) {
      gene_data <- getGeneData(term)
      tab_idx <- as.numeric(input$geneTable_rows_selected)
      if(length(tab_idx) > 0) {
        gd <- gene_data[tab_idx, ]
        if(!is.null(gd)) sel <- gd %>% pull(gene_id)
      }
    }
    return(sel)
  }

  output$termTable <- DT::renderDataTable({
    terms <- getTermData() %>% 
      select(term, term_name, NES) %>% 
      mutate_if(is.numeric, ~signif(.x, 2)) %>% 
      mutate(term = as.character(term))
    terms
  }, selection = "single", escape = FALSE)
  
  
  output$geneTable <- DT::renderDataTable({
    term <- selectTerm()
    df <- NULL
    if(!is.null(term)) {
      gene_data <- getGeneData(term)
      if(!is.null(gene_data)) {
        df <- gene_data %>% 
          select(gene_name, description, gene_biotype, logFC) %>% 
          mutate_if(is.numeric, ~signif(.x, 2))
      }
    }
    df
  }, selection = "single")


  output$genePlot <- renderPlot({
    dat <- getData()
    sel <- selectGene()
    if(!is.null(sel)) {
      sh_plot_genes(data$dat, data$metadata, sel, input$intensityScale)
    }
  })
  
  output$geneInfo <- function() {
    dat <- getData()
    sel <- selectGene()
    if(!is.null(sel)) {
      d <- data$res %>% 
        filter(gene_id %in% sel) %>% 
        mutate_if(is.numeric, ~signif(.x, 2)) %>% 
        select(logFC, logCPM, PValue, FDR, contrast)
      d %>% 
        kable("html") %>% 
        kable_styling("striped", full_width = F) %>%
        row_spec(which(d$FDR < 0.05 & abs(d$logFC) > 1), bold=TRUE, background="lightsalmon")
        
    }
  }
  
}

# Run the application
shinyApp(ui = ui, server = server)
