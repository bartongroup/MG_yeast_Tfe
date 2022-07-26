### EDGER EXPLORER

libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/4.1"
if (dir.exists(libDir)) .libPaths(libDir)

library(scales)
library(shiny)
#library(shinycssloaders)
library(tidyverse)
library(DT)
source("../shiny_func.R")

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

data <- sh_read_edger_data("../data")
initial_contrasts <- unique(data$de$sel$contrast) %>% as.character()

max_points <- 1000


#######################################################################

ui <- shinyUI(fluidPage(
  
  tags$style(css),
  
  titlePanel("Tfe1/Tfe2 toxin expression in yeast: differential expression"),

  fluidRow(
    column(12,
      fluidRow(
        column(4,
          radioButtons("method", "Method:", choices = c("Pairwise" = "sel", "Full model" = "fi", "Tfe correlation" = "tfe"), inline = TRUE),
          conditionalPanel(
            condition = "input.method != 'tfe'",
            selectInput("contrast", "Contrast:", choices = initial_contrasts),
            radioButtons("plot_type", "Plot type:", choices = c("Volcano" = "vol", "MA" = "ma"), inline = TRUE)
          ),
          plotOutput("main_plot", height = "480px", width = "100%", brush = "plot_brush", hover = "plot_hover")
        ),
        column(3,
          radioButtons("intensity_scale", "Intesity scale:", choices = c("lin" = "", "log" = "log"), inline = TRUE),
          plotOutput("gene_plot", height = "400px",width = "100%")
        ),
        column(5,
          p("Gene list"),
          div(style = 'height: 200px; overflow-y: scroll', tableOutput("gene_info")),
          br(),
          radioButtons("enrichment", "Enrichment:", choices = c("GO" = "go", "Reactome" = "re", "KEGG" = "kg"), inline = TRUE),
          div(style = 'height: 400px; overflow-y: scroll', tableOutput("enrichment")),
        )
      ),
      fluidRow(
        DT::dataTableOutput("all_gene_table")
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
  
  # Available contrasts
  get_contrasts <- reactive({
    ctrs <- NULL
    if (input$method != 'tfe') {
      ctrs <- data$de[[input$method]]$contrast %>% 
        as.character() %>% 
        unique()
    }
    ctrs
  })

  # Update contrasts based on method
  observeEvent(input$method, {
    freezeReactiveValue(input, "contrast")
    updateSelectInput(session = session, inputId = "contrast", choices = get_contrasts())
  })
  
  
  get_de <- function() {
    ctr <- input$contrast
    data$de[[input$method]] %>% 
      filter(contrast == ctr)
  }
  
  get_xy_data <- function() {
    if (input$method == "tfe") {
      xy_data <- data$tfe_cor %>% 
        mutate(x = atanh(corTfe1), y = atanh(corTfe2))
    } else {
      de <- get_de()
      if (input$plot_type == "vol") {
        xy_data <- de %>% 
          mutate(x = logFC, y = -log10(PValue))
      } else if (input$plot_type == "ma") {
        xy_data <- de  %>% 
          mutate(x = logCPM, y = logFC)
      }
    }
    xy_data
  }
  
  select_gene <- function(max_hover = 1) {
    xy_data <- get_xy_data()
    sel <- NULL
    tab_idx <- as.numeric(input$all_gene_table_rows_selected)
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
  
  output$gene_info <- renderTable({
    xy_data <- get_xy_data()
    sel <- select_gene()
    df <- NULL
    if (!is.null(sel) && length(sel) >= 1 && length(sel) <= max_points) {
      df <- xy_data %>%
        filter(gene_id %in% sel) %>% 
        arrange(gene_name)
      if (input$method == "tfe") {
        df <- df %>% select(gene_name, gene_biotype, description)
      } else {
        df <- df %>% select(gene_name, gene_biotype, description, FDR)
      }
    } else if (length(sel) > max_points) {
      df <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
    }
    df
  })

  enrichment_table <- function(terms) {
    xy_data <- get_xy_data()
    sel <- NULL
    fe <- NULL
    if (!is.null(input$plot_brush)) {
      brushed <- brushedPoints(xy_data, input$plot_brush)
      sel <- brushed$gene_id
      n <- length(sel)
      if (n > 0 && n <= max_points) {
        fe <- sh_functional_enrichment(data$all_genes, sel, terms, gene2name = data$gene2name)
      } else if (n > 0) {
        fe <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  }
  
  output$enrichment <- renderTable({
    d <- data$terms[[input$enrichment]]
    enrichment_table(d)
  })
  
  
  output$gene_plot <- renderPlot({
    sel <- select_gene()
    if (!is.null(sel) && length(sel) > 0 && length(sel) <= max_points) {
      sh_plot_genes(data$star, sel, input$intensity_scale)
    }
  })
  
  output$main_plot <- renderPlot({
    xy_data <- get_xy_data()
    tab_idx <- as.numeric(input$all_gene_table_rows_selected)
    
    if (input$method == "tfe") {
      g <- sh_plot_tfe(xy_data)
    } else {
      if (input$plot_type == "vol") {
        g <- sh_plot_volcano(xy_data)
      } else if (input$plot_type == "ma") {
        g <- sh_plot_ma(xy_data)
      }
    }
    if (length(tab_idx) >= 1) {
      g <- g + geom_point(data = xy_data[tab_idx, ], colour = "red", size = 2)
    }
    g
  })

  output$all_gene_table <- DT::renderDataTable({
    if (input$method == "tfe") {
      d <- get_xy_data() %>%
        select(gene_name, corTfe1,corTfe2, description) %>% 
        mutate_if(is.numeric, ~signif(.x, 3))
    } else {
      d <- get_xy_data() %>%
        select(gene_name, gene_biotype, logFC, FDR, description) %>% 
        mutate_if(is.numeric, ~signif(.x, 3))
    }
    DT::datatable(d, class = 'cell-border strip hover', selection = "single", rownames = FALSE)
  })
}

# Run the application
shinyApp(ui = ui, server = server)


# For testing
# input <- list(method = "fi", plot_type = "vol", contrast = "strainTfe2", enrichment = "go", intensity_scale = "lin")
