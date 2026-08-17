library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(tidyr)

# ---- Load data ----
env_med <- new.env()
load("all_res.RData", envir = env_med)
res_med <- env_med$shiny_res

env_bb <- new.env()
load("bb_all_res.RData", envir = env_bb)
res_bb <- env_bb$shiny_res

res_med <- res_med %>% mutate(seed = as.character(seed))
res_bb  <- res_bb  %>% mutate(seed = as.character(seed))
res <- bind_rows(res_med, res_bb)

ref_choices <- sort(unique(res$ref[res$method == "GGA-BC"]))
strata_choices <- sort(unique(res$strata))
N_choices      <- sort(unique(res$N))

strata_real_min <- min(res$strata_real, na.rm = TRUE)
strata_real_max <- max(res$strata_real, na.rm = TRUE)

outcome_choices <- c(
  "CV (total)"          = "cv_tot",
  "Sample size (n_tot)" = "n_tot",
  "Time"                = "time",
  "Trace"               = "trace"
)

ui <- fluidPage(
  
  # ---- CSS for sticky footer ----
  tags$head(
    tags$style(HTML("
      html, body {
        height: 100%;
        margin: 0;
      }
      #page-container {
        min-height: 100vh;
        display: flex;
        flex-direction: column;
      }
      #content-wrap {
        flex: 1;
      }
      #footer {
        text-align: center;
        font-size: 12px;
        color: grey;
        padding: 10px;
        border-top: 1px solid #ddd;
        background-color: #fafafa;
      }
    "))
  ),
  
  div(id = "page-container",
      
      # ---- MAIN CONTENT ----
      div(id = "content-wrap",
          
          # ---- Title + Authors ----
          div(
            style = "text-align: center; margin-bottom: 10px;",
            tags$h2("Multivariate stratification results explorer"),
            tags$p(
              "Georgi Borros, Şebnem Er, Sulaiman Salau",
              style = "font-size: 15px; color: #555; margin-top: -10px;"
            )
          ),
          
          sidebarLayout(
            div(
              style = "padding-top: 20px;",
              sidebarPanel(
                selectInput("method", "Method", choices = sort(unique(res$method)), multiple = TRUE,
                            selected = unique(res$method)),
                selectInput("dataset", "Dataset", choices = sort(unique(res$dataset)), multiple = TRUE,
                            selected = unique(res$dataset)),
                selectInput("ref", "Ref (GGA-BC only)", choices = ref_choices, multiple = TRUE,
                            selected = ref_choices),
                hr(),
                
                checkboxInput("filter_strata_real", "Filter by realised number of strata? (varies for GGA-BC)", value = FALSE),
                conditionalPanel(
                  condition = "input.filter_strata_real == true",
                  sliderInput("strata_real", "Number of strata (realised)",
                              min = strata_real_min, max = strata_real_max,
                              value = c(strata_real_min, strata_real_max), step = 1)
                ),
                hr(),
                selectInput("outcome", "Outcome to plot", choices = outcome_choices, selected = "cv_tot"),
                hr(),
                downloadButton("download_filtered", "Download filtered data (.csv)")
              )
            ),
            
            mainPanel(
              width = 7,
              
              tabsetPanel(
                id = "main_tabs",
                
                # ---- PLOTS ----
                tabPanel(
                  "Outcome plots",
                  icon = icon("chart-line"),
                  
                  tags$div(
                    style = "margin-top: 50px; margin-bottom: 20px;",
                    tags$h4("Distribution of outcomes"),
                    tags$p("Top: varies precribed strata (n = 500). Bottom: varies precribed sample size (strata = 5).")
                  ),
                  
                  plotOutput("outcome_plot_strata", height = "500px"),
                  br(),
                  plotOutput("outcome_plot_N", height = "500px")
                ),
                
                # ---- SUMMARY TABLE ----
                tabPanel(
                  "Summary table",
                  icon = icon("table"),
                  
                  tags$div(style = "margin-top: 15px;",
                           selectInput("table_strata", "Filter strata", choices = strata_choices,
                                       selected = strata_choices, multiple = TRUE),
                           selectInput("table_N", "Filter sample size (n)", choices = N_choices,
                                       selected = N_choices, multiple = TRUE)
                  ),
                  
                  tags$p("Aggregated across seeds for the current filter selection."),
                  DTOutput("summary_table"),
                  br(),
                  downloadButton("download_summary", "Download summary table (.csv)")
                ),
                
                # ---- RAW TABLE ----
                tabPanel(
                  "Raw filtered data",
                  icon = icon("list"),
                  
                  tags$div(style = "margin-top: 15px;",
                           selectInput("raw_strata", "Filter strata", choices = strata_choices,
                                       selected = strata_choices, multiple = TRUE),
                           selectInput("raw_N", "Filter sample size (n)", choices = N_choices,
                                       selected = N_choices, multiple = TRUE)
                  ),
                  
                  DTOutput("raw_table")
                )
              )
            )
          )
      ),
      
      # ---- FOOTER ----
      div(
        id = "footer",
        paste("Last updated:", format(Sys.time(), "%d %B %Y, %H:%M"))
      )
  )
)

server <- function(input, output, session) {
  
  base_filtered <- reactive({
    df <- res %>%
      filter(
        method %in% input$method,
        dataset %in% input$dataset,
        method != "GGA-BC" | ref %in% input$ref
      )
    
    if (isTRUE(input$filter_strata_real)) {
      df <- df %>%
        filter(
          strata_real >= input$strata_real[1],
          strata_real <= input$strata_real[2]
        )
    }
    
    df
  })
  
  table_filtered <- reactive({
    base_filtered() %>%
      filter(strata %in% input$table_strata,
             N %in% input$table_N)
  })
  
  raw_filtered <- reactive({
    base_filtered() %>%
      filter(strata %in% input$raw_strata,
             N %in% input$raw_N)
  })
  
  output$outcome_plot_strata <- renderPlot({
    df <- base_filtered() %>% filter(N == 500)
    validate(need(nrow(df) > 0, "No data for N = 500."))
    
    outcome_col <- input$outcome
    outcome_label <- names(outcome_choices)[outcome_choices == outcome_col]
    
    ggplot(df, aes(x = factor(strata), y = .data[[outcome_col]], colour = method)) +
      geom_boxplot(position = position_dodge(width = 0.7)) +
      facet_wrap(~ dataset) +
      scale_colour_brewer(palette = "Set1") +
      labs(x = "Prescribed number of strata (L)", y = outcome_label,
           title = paste(outcome_label, "(N = 500)")) +
      theme_bw(base_size = 15)
  })
  
  output$outcome_plot_N <- renderPlot({
    df <- base_filtered() %>% filter(strata == 5)
    validate(need(nrow(df) > 0, "No data for strata = 5."))
    
    outcome_col <- input$outcome
    outcome_label <- names(outcome_choices)[outcome_choices == outcome_col]
    
    ggplot(df, aes(x = factor(N), y = .data[[outcome_col]], colour = method)) +
      geom_boxplot(position = position_dodge(width = 0.7)) +
      facet_wrap(~ dataset) +
      scale_colour_brewer(palette = "Set1") +
      labs(x = "Presribed sample size (n)", y = outcome_label,
           title = paste(outcome_label, "(strata = 5)")) +
      theme_bw(base_size = 15)
  })
  
  output$summary_table <- renderDT({
    df <- table_filtered()
    validate(need(nrow(df) > 0, "No data for this combination of filters."))
    
    summary_df <- df %>%
      group_by(method, dataset, ref, strata, strata_real, N) %>%
      summarise(
        n_runs      = n(),
        mean_cv_tot = mean(cv_tot, na.rm = TRUE),
        mean_n_tot  = mean(n_tot, na.rm = TRUE),
        mean_time   = mean(time, na.rm = TRUE),
        mean_trace  = mean(trace, na.rm = TRUE),
        .groups = "drop"
      )
    
    numeric_cols <- names(summary_df)[sapply(summary_df, is.numeric)]
    
    datatable(summary_df, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = numeric_cols, digits = 3)
  })
  
  output$raw_table <- renderDT({
    datatable(raw_filtered(), options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })
  
  output$download_filtered <- downloadHandler(
    filename = function() paste0("filtered_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(base_filtered(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)