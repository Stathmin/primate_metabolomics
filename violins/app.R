


library('shiny')
library('tidyverse')
library('stringr')
library('ggplot2')
library('multcompView')
library('broom')
library('flextable')
library('here')
library('moments')
library('bslib')

set.seed(42)

custom_theme <- bs_theme(
  version = 5,
  bg = "#FFFFFF",
  fg = "#000000",
  primary = "#0199F8",
  secondary = "#FF374B",
  base_font = font_google("Fira Code", local = TRUE))

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
  colnames(Tukey.labels) = c('letter','group')
  return(Tukey.labels)
}

merged_df <- read_rds(str_interp('${here()}/.cache/merged_df.rds'))
special_names_table <- read_rds(
  str_interp('${here()}/.cache/special_names_table.rds')) %>% 
  mutate(special_name = paste(common_name,special_name, sep = '_'))

merged_df <- merged_df %>% 
  rename_with(\(x) special_names_table %>% 
                filter(common_name == !!x) %>% 
                select(special_name) %>% pull,
              !!special_names_table$common_name)

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
  theme = custom_theme,
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
      selectizeInput('out_variables',
                  'Observable:',
                  choices = out_variables[1:5],
                  multiple = FALSE,
                  selected = out_variables[1]),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      uiOutput("formula"),
      plotOutput("distPlot"),
      tabsetPanel(
        tabPanel(
          "Discriptive",
          uiOutput("DISCvarNames"),
          verbatimTextOutput("DiscriptiveTable"),
          fluidRow(uiOutput("discriptiveTable_flex"))
        ),
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
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'out_variables', 
                       choices = out_variables, 
                       selected = out_variables[1],
                       server = TRUE) #debounce_work_around
  
  out_variables_d <- debounce(reactive({input$out_variables}), 500)
  
  local_model <- reactive(formulate(out_variables_d(),
                                    sort(input$grouping_factors)))
  
  local_discriptive <- reactive(
    merged_df %>%
      filter(!is.na(!!out_variables_d())) %>% 
      unite(group, !!sort(c('gender',  'species')),
            sep = ":", remove = FALSE) %>% 
      select(group, !!out_variables_d()) %>% 
      group_by(group) %>% 
      summarise(n = n(),
                median = median(!!sym(out_variables_d())),
                mean = mean(!!sym(out_variables_d())),
                cv_perc = 100*sd(!!sym(out_variables_d()))/median(!!sym(out_variables_d())),
                min = min(!!sym(out_variables_d())),
                skewness = skewness(!!sym(out_variables_d())),
                kurtosis = kurtosis(!!sym(out_variables_d())),
      ) %>% 
      ungroup() %>% 
      arrange(mean)
  )
  
  local_anova <-
    reactive(aov(as.formula(local_model()), data = merged_df))
  
  tukey_needed <- reactive({
    as.formula(local_model()) %>%
      all.vars() %>%
      length() > 1
  })
  
  letters_needed <- reactive({tukey_needed() & 
      length(vecsets::vsetdiff(sort(input$x_axis),sort(input$grouping_factors))) == 0
  })
  
  
  local_tukey <- reactive(if (tukey_needed()) {
    TukeyHSD(local_anova()) %>% 
      purrr::modify_depth(5, \(x) replace_na(x, 
                                              replace = 0),
                          .ragged = TRUE)
  } else {
    NA
  })
  
  local_names <- reactive(if (letters_needed()) {
    generate_label_df(local_tukey(), paste(sort(input$x_axis),
                                           sep = ':',
                                           collapse = ':'))
  } else {
    NA
  })
  
  local_df <- reactive(if (letters_needed()) {
    merged_df %>%
      filter(!is.na(!!out_variables_d())) %>% 
      unite(group, !!sort(input$x_axis),
            sep = ":", remove = FALSE) %>% 
      left_join(local_names(), by=c('group' = 'group'))
  } else {
    merged_df %>%
      unite(group, !!sort(input$x_axis),
            sep = ":", remove = FALSE) %>% 
      mutate(letter = '-')
  }
  )
  
  output$formula <- renderUI(HTML(
    as.character(div(
      style="text-align: center; font-weight: bold;",local_model()))))
  
  output$discriptiveTable_flex <- renderUI({
    local_discriptive() %>%
      flextable() %>% 
      colformat_double(
        big.mark = ",", digits = 2, na_str = "N/A"
      ) %>% 
      autofit() %>% 
      htmltools_value()
  })
  
  output$anovaTable_flex <- renderUI({
    tidy(local_anova()) %>%
      mutate(sig = case_when(p.value <= 0.001 ~ "***",
                             p.value <= 0.01 ~ "**",
                             p.value <= 0.05 ~ "*", 
                             p.value <= 0.1 ~ ".",
                             .default = "")) %>%
      flextable() %>% 
      colformat_double(
        big.mark = ",", digits = 2, na_str = "N/A"
      ) %>% 
      autofit() %>% 
      htmltools_value()
  })
  output$tukeyTable_flex <- renderUI(if (tukey_needed()) {
    tidy(local_tukey()) %>%
      mutate(sig = case_when(adj.p.value <= 0.001 ~ "***",
                             adj.p.value <= 0.01 ~ "**",
                             adj.p.value <= 0.05 ~ "*", 
                             adj.p.value <= 0.1 ~ ".",
                             .default = "")) %>% 
      flextable() %>%
      colformat_double(
        big.mark = ",", digits = 2, na_str = "N/A"
      ) %>% 
      autofit() %>% 
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
    ggplot(local_df(),
           aes(
             x = group,
             y = .data[[out_variables_d()]],
             color = letter
           )) +
      geom_boxplot(notch = FALSE) +
      theme_minimal()
    
  })

}

# Run the application
shinyApp(ui = ui, server = server)

