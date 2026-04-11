library(shiny)
library(bslib)
library(raretrans)
library(popbio)
library(ggplot2)

# Load data into the environment
data("L_elto", package = "raretrans")

# Pre-defined informative prior (from transition_priors vignette)
RLT_Tprior <- matrix(
  c(0.25, 0.025, 0.0,
    0.05, 0.9,   0.025,
    0.01, 0.025, 0.95,
    0.69, 0.05,  0.025), 
  byrow = TRUE, nrow = 4, ncol = 3
)

# Example selections from the vignette
examples <- list(
  "Pop 905, Year 9 (\u03BB = 0)"      = c(pop = 905, year = 9),
  "Pop 914, Year 6 (\u03BB = 1)"      = c(pop = 914, year = 6),
  "Pop 250, Year 5 (\u03BB = 0.93)"   = c(pop = 250, year = 5)
)

ui <- page_sidebar(
  title = "Bayesian Priors for Matrix Population Models",
  theme = bs_theme(version = 5, preset = "united"),
  
  sidebar = sidebar(
    title = "Configuration",
    
    selectInput(
      "example", 
      "Select Example (L. elto):", 
      choices = names(examples)
    ),
    
    radioButtons(
      "prior_type", 
      "Prior Type:", 
      choices = c("Uniform (Uninformative)" = "uniform", 
                  "RLT (Informative)" = "rlt")
    ),
    
    conditionalPanel(
      condition = "input.prior_type == 'rlt'",
      radioButtons(
        "prior_weight", 
        "Prior Weight:",
        choices = c("Weight = 1" = "1",
                    "Weight = 0.5 * N" = "0.5N",
                    "Weight = 1.0 * N" = "N")
      )
    ),
    
    hr(),
    h6("Population Summary"),
    uiOutput("pop_summary")
  ),
  
  card(
    card_header("Posterior Transition Distributions (Beta Density)"),
    plotOutput("density_plot", height = "500px"),
    full_screen = TRUE
  ),
  
  card(
    card_header("Adjusted Transition Matrix (Expected Values)"),
    tableOutput("matrix_table")
  )
)

server <- function(input, output, session) {
  
  # Filter the data for the selected example
  current_data <- reactive({
    req(input$example)
    ex <- examples[[input$example]]
    subset(raretrans::L_elto, POPNUM == ex["pop"] & year == ex["year"])
  })
  
  # Calculate TF and N
  tf_and_n <- reactive({
    dat <- current_data()
    # Ensure correct factors for stages to populate the full matrix
    dat$stage <- factor(dat$stage, levels = c("p", "j", "a", "m"))
    dat$fate  <- factor(dat$next_stage, levels = c("p", "j", "a", "m"))
    
    TF <- popbio::projection.matrix(
      as.data.frame(dat), 
      stage = "stage", fate = "fate", 
      fertility = "fertility", sort = c("p", "j", "a"), TF = TRUE
    )
    N <- raretrans::get_state_vector(
      as.data.frame(dat), 
      stage = "stage", sort = c("p", "j", "a")
    )
    
    list(TF = TF, N = N)
  })
  
  # Extract the correct prior configuration
  prior_settings <- reactive({
    if (input$prior_type == "uniform") {
      list(P = NULL, weight = -1)
    } else {
      w <- switch(input$prior_weight,
                  "1" = -1,    # using raw unscaled P is weight=1
                  "0.5N" = 0.5,
                  "N" = 1)
      list(P = RLT_Tprior, weight = w)
    }
  })
  
  # Observe inputs and render density plot
  output$density_plot <- renderPlot({
    inputs <- tf_and_n()
    priors <- prior_settings()
    
    raretrans::plot_transition_density(
      TF = inputs$TF, 
      N = inputs$N, 
      P = priors$P,
      priorweight = priors$weight,
      stage_names = c("seedling", "juvenile", "adult")
    )
  })
  
  output$matrix_table <- renderTable({
    inputs <- tf_and_n()
    priors <- prior_settings()
    
    # fill_transitions returns a transition matrix T when returnType = "T"
    filled_T <- raretrans::fill_transitions(
      TF = inputs$TF,
      N = inputs$N,
      P = priors$P,
      priorweight = priors$weight,
      returnType = "T"
    )
    
        # Combine with fertility to get the full A matrix
    # for simplicity, we focus mostly on T for differences in this app
    A_mat <- filled_T + inputs$TF$F
    
    rownames(A_mat) <- c("seedling", "juvenile", "adult")
    colnames(A_mat) <- c("seedling", "juvenile", "adult")
    
    A_mat
  }, rownames = TRUE, align = "c", digits = 4)
  
  output$pop_summary <- renderUI({
    inputs <- tf_and_n()
    priors <- prior_settings()
    
    # Calculate lambdas
    lambda_raw <- format(popbio::lambda(inputs$TF$T + inputs$TF$F), digits = 4)
    
    filled_T <- raretrans::fill_transitions(
             TF = inputs$TF, N = inputs$N, 
             P = priors$P, priorweight = priors$weight, returnType = "T"
    )
    unif_gamma <- matrix(rep(c(NA_real_, NA_real_, 0.00001), 3), byrow=TRUE, nrow=3)
    # the uninformative gamma prior is used for F
    filled_F <- raretrans::fill_fertility(
      TF = inputs$TF, N = inputs$N, 
      alpha = unif_gamma, beta = unif_gamma, priorweight = -1, returnType = "F"
    )
    
    lambda_adj <- format(popbio::lambda(filled_T + filled_F), digits = 4)
    
    HTML(paste0(
      "<b>Actual Sample Size (N):</b> ", sum(inputs$N), "<br/>",
      "<b>Raw Asymptotic Growth (\u03BB):</b> ", lambda_raw, "<br/>",
      "<b>Adjusted \u03BB (Expected):</b> ", lambda_adj
    ))
  })
}

shinyApp(ui, server)
