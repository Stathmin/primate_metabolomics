


library('shiny')
library('tidyverse')
library('stringr')
library('ggplot2')
library('multcompView')
library('broom')
library('flextable')
library("here")

set.seed(42)

formulate <- function(RHS, LHF) {
  if (length(LHF) == 0) {
    LHS = '1'
  } else if ((length(LHF) == 1)) {
    LHS = str_interp('1 + ${LHF}')
  } else {
    LHS = str_interp('1 + (${paste(LHF, collapse=" + ")})^2')
  }
  return(paste(RHS, LHS, sep = " ~ "))
}

generate_label_df <- function(TUKEY, variable) {
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][, 4]
  Tukey.labels <-
    data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

merged_df <- read_rds(str_interp('${here()}/.cache/merged_df.rds'))
merged_df %>%
  select(-individual_id) %>%
  select(where((\(x) {
    !is.numeric(x) & (length(x) / n_distinct(x) > 1)
  }))) %>%
  colnames -> grouping_factors

merged_df %>%
  select(where(is.numeric)) %>%
  colnames -> out_variables

# Define UI for application that draws a histogram
ui <- fluidPage(# Application title
  titlePanel("Planet of the Apes"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput(
        'grouping_factors',
        'Grouping factors:',
        grouping_factors,
        multiple = TRUE
      ),
      selectInput(
        'x_axis',
        'X axis:',
        grouping_factors,
        multiple = TRUE,
        selected = grouping_factors[1]
      ),
      selectInput('out_variables',
                  'Observable:',
                  out_variables,
                  multiple = FALSE),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      textOutput("formula"),
      plotOutput("distPlot"),
      tabsetPanel(
        tabPanel(
          "ANOVA",
          uiOutput("ANOVAvarNames"),
          verbatimTextOutput("anovaTable"),
          fluidRow(uiOutput("anovaTable_flex"))
        ),
        tabPanel(
          "Tukey",
          uiOutput("TUKEYvarNames"),
          verbatimTextOutput("tukeyTable"),
          fluidRow(uiOutput("tukeyTable_flex"))
        ),
        tabPanel(
          "Group Letters",
          uiOutput("TUKEYLvarNames"),
          verbatimTextOutput("tukeyLetters"),
          fluidRow(uiOutput("tukeyLetters_flex"))
        ),
      ),
    )
  ))

# Define server logic required to draw a histogram
server <- function(input, output) {
  local_model <- reactive(formulate(input$out_variables,
                                    input$grouping_factors))
  
  local_anova <-
    reactive(aov(as.formula(local_model()), data = merged_df))
  
  tukey_needed <- reactive({
    as.formula(local_model()) %>%
      all.vars() %>%
      length() > 1
  })
  
  letters_needed <- reactive({tukey_needed() & 
      length(vecsets::vsetdiff(input$x_axis,input$grouping_factors)) == 0
  })
  
  
  local_tukey <- reactive(if (tukey_needed()) {
    TukeyHSD(local_anova())
  } else {
    NA
  })
  
  local_names <- reactive(if (letters_needed()) {
    generate_label_df(local_tukey(), paste(input$x_axis,
                                           sep = ':',
                                           collapse = ':'))
  } else {
    NA
  })
  
  output$formula <- renderText(as.character(local_model()))
  
  output$anovaTable_flex <- renderUI({
    flextable(tidy(local_anova())) %>% htmltools_value()
  })
  output$tukeyTable_flex <- renderUI(if (tukey_needed()) {
    flextable(tidy(local_tukey())) %>%
      htmltools_value()
  } else {
    ''
  })
  output$tukeyLetters_flex <- renderUI(if (letters_needed()) {
    flextable(local_names()) %>%
      htmltools_value()
  } else {
    ''
  })
  
  output$distPlot <- renderPlot({
    ggplot(merged_df,
           aes(
             x = merged_df %>%
               select(!!input$x_axis) %>%
               unite(group, !!input$x_axis,
                     sep = ":") %>%
               pull(),
             y = .data[[input$out_variables]]
           )) +
      geom_boxplot(notch = TRUE) +
      xlab("group") +
      theme_minimal()
    
  })
}

# Run the application
shinyApp(ui = ui, server = server)
