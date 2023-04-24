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


merged_df <- read_rds(str_interp('.cache/merged_df.rds'))
special_names_table <- read_rds(
  str_interp('.cache/special_names_table.rds')) %>%
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

out_variables %>% write_rds('violins/not_distinct.rds')
special_names_table$special_name %>% write_rds('violins/distinct_special_names.rds')

# # all primates are distinct -----------------------------------------------


check_distinct <- function(x, .any=FALSE) {
  if (length(x) == 1) {
    stop("faulty grouping, group of size 1")
  } else {
  combn(x, 2) %>%
    t() %>%
    data.frame() %>%
    tibble() %>%
    mutate(different =
             purrr::map2_vec(str_split(X1,''),
                             str_split(X2,''),
                             ~length(intersect(.x,.y)) < 1)
    ) %>%
    pull(different) %>%
    {ifelse(.any,
           any(.),
           all(.))}
  }
  }

enlabel_species_find_distinct <- function(x,
                                          GroupInd = '',
                                          AnovaFactors = c('gender',
                                                           'tissue',
                                                           'species'),
                                          .any_trunk = TRUE,
                                          .any_leaf = TRUE) {
  x %>%
    formulate(., AnovaFactors) %>%
    as.formula() %>%
    aov(data = merged_df) %>%
    TukeyHSD() %>%
    generate_label_df(GroupInd) %>%
    as_tibble() %>%
    rename(!!GroupInd := group) %>%
    separate(!!GroupInd,
             sep = ':',
             into = str_split_1(GroupInd, ':')) %>%
    {
      if (length(str_split_1(GroupInd, ':')) > 1) {
        group_by(., across(all_of(!!({
          str_split_1(GroupInd, ':') %>%
            head(., -1)
        }))))
      } else {
        .
      }
    } %>%
    mutate(are_distinct = check_distinct(letter, .any = .any_leaf)) %>%
    ungroup() %>%
    pull(are_distinct) %>%
    {
      ifelse(.any_trunk == TRUE,
             any(.),
             all(.))
    }
}

# #test
# merged_df %>%
#   summarise(lc_pos_peak_041053 = list(
#     enlabel_species_find_distinct('lc_pos_peak_041053',
#                                   GroupInd = 'species',
#                                   AnovaFactors = c('gender', 'tissue', 'species')))
#     ) %>%
#   unnest(lc_pos_peak_041053)


merged_df %>%
  summarise(across(!!out_variables, ~ list(
    enlabel_species_find_distinct(cur_column(),
                                  GroupInd = 'species',
                                  AnovaFactors = c('gender', 'species', 'tissue'),
                                  .any_trunk = FALSE, .any_leaf = FALSE))
    )) %>%
  unnest(!!out_variables) %>%
  select(where(\(x) x[1] == TRUE)) %>%
  colnames() %>% write_rds('violins/distinct_species.rds')

distinct_species_variables <- read_rds('violins/distinct_species.rds')

merged_df %>%
  summarise(across(!!distinct_species_variables, ~ list(
    enlabel_species_find_distinct(cur_column(),
                                  GroupInd = 'tissue:species',
                                  AnovaFactors = c('gender', 'tissue', 'species'),
                                  .any_trunk = TRUE, .any_leaf = FALSE))
  )) %>%
  unnest(!!distinct_species_variables) %>%
  select(where(\(x) x[1] == TRUE)) %>%
  colnames() %>% write_rds('violins/distinct_species_at_all_tissues.rds')

distinct_species_at_all_tissues_variables <- read_rds('violins/distinct_species_at_all_tissues.rds')

file.create('.SUCCESS_TOKEN')
