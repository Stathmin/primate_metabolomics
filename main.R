
# imports -----------------------------------------------------------------
renv::activate()
renv::hydrate(prompt = FALSE)

library('tidyverse')
library('readxl')
library('janitor')
library('purrr')
library('stringi')
library('checkmate')

# custom_fucntions --------------------------------------------------------

pivot_transpose <- function(
    data,
    cols_to_wide = 1,
    cols_to_long = 2:last_col(),
    names_to = "name",
    .drop = FALSE
) {
  
  cols_to_long_expr <- enquos(cols_to_long)
  cols_to_wide_expr <- enquos(cols_to_wide)
  
  # drop unused columns
  data_select <- data
  if (.drop) {
    data_select <- data_select %>%
      select(!!! cols_to_wide_expr, !!! cols_to_long_expr)
  }
  
  data_select %>%
    pivot_longer(
      cols      = (!!! cols_to_long_expr),
      values_to = "...value",
      names_to  = names_to
    ) %>%
    pivot_wider(
      names_from  = (!!! cols_to_wide_expr),
      values_from = ...value
    )
  
}

#functions to find suspicious values
.return_only_nonumeric_strings <- function(x){
  y <- as.character(x)
  if (stri_count_regex(y,"^[0-9\\.]+$")==1){
    return(NA_character_)
  }else{
    return(x)
  }
}

return_only_nonumeric_strings <- function(X){
  map(X, .return_only_nonumeric_strings) %>% unlist
}

# file_reading_expression -------------------------------------------------

expression_file <- './data/pbio.1001871.s022.xlsx'

exp_df <- readxl::read_xlsx(expression_file) # initial row contains bio_sample_data
remove(expression_file)

exp_init_dim <- exp_df %>% dim() # 4700  200 (!)

exp_df <- exp_df %>% pivot_transpose # transpose to make more logical indexing
remove(pivot_transpose)

exp_pivot_dim <- exp_df %>% dim() # 199 4701 (!)

assert_true(sum(exp_init_dim) == sum(exp_pivot_dim))
remove(exp_init_dim)

exp_df %>% slice(1:5) %>% 
  select(1:5) #upper-left chunk

exp_df %>% 
  tail(5) %>% 
  select(tail(colnames(.), 5)) #lower_right chunk

#l_r chunk shows strange last row, needs checking
exp_df %>% 
  tail(1) %>% 
  select(-where(is.na)) %>% 
  unlist(use.names = FALSE) %>% 
  head() #contains alternative names for columns where applicable

#NO IDEA HOW TO WORK WITH ALT COL NAMES!!!
#so lets create separate table for them

special_names_table <- exp_df %>%  #table of names for genes and loci
  tail(1) %>% 
  select(-where(is.na)) %>% 
  pivot_longer(everything(),
               names_to = 'common_name',
               values_to = 'special_name') %>% 
  slice(2:nrow(.)) %>%             #because row 1 is c('name', 'annotation')
  mutate(across(everything(), 
                make_clean_names)) #for compatibility with merged_df

#now drop last row with non-numeric data
square_exp_df <- exp_df %>% 
  head(-1)
remove(exp_df)

#check for strange values in columns after first three
sus_values_mask <- square_exp_df %>% 
  select(-(1:3)) %>%
  mutate(across(where(is.character), \(x) replace_na(x,''))) %>% 
  mutate(across(everything(),return_only_nonumeric_strings))

assert_true(sus_values_mask %>% unlist %>% unname %>% unique %>% is.na)
remove(sus_values_mask)

square_exp_df <- square_exp_df %>% mutate(across(
  -c('name','species','tissue'),
  as.double)) 

square_exp_dim <- dim(square_exp_df)
# assert only annotation removed (stored in special_names_table)
assert_true(sum(exp_pivot_dim - square_exp_dim) == 1)
remove(list = c('square_exp_dim', 'exp_pivot_dim'))
# sample_file_reading -----------------------------------------------------

bio_samples_file <- './data/pbio.1001871.s016.xlsx'

bio_df <- readxl::read_xlsx(bio_samples_file) # initial row contains bio_sample_data
remove(bio_samples_file)

bio_df %>% #`years of age`, `days of age` are numeric
  mutate(
    across(c(`years of age`, `days of age`), as.double),
    `individual id` = as.factor(as.integer(`individual id`))
    ) #creates some NAs, lets check why

bio_df %>% 
  select(c(`years of age`, `days of age`)) %>%
  mutate_if(is.character, ~replace_na(.,"")) %>% 
  mutate(across(everything(),return_only_nonumeric_strings)) %>% 
  unlist %>% unname %>% unique # has "NA", no problems with introduced NAs
remove(return_only_nonumeric_strings)

bio_df <- bio_df %>% 
  mutate_if(is.character, ~na_if(., 'NA')) %>% 
  mutate(
    across(c(`years of age`, `days of age`), as.double)
  )

#some ids are not unique
bio_df %>% group_by(`individual id`) %>% 
  summarise(n = n()) %>% 
  filter(n == max(n)) %>% 
  head(1) %>% .[[1, "individual id"]] -> non_unique_id

bio_df %>% 
  filter(`individual id` == non_unique_id) #OK, some samples come from same sources

bio_df %>% 
  select_if(function(x) any(is.na(x))) %>% 
  summarise(across(everything(), \(x) sum(is.na(x)))) #days of age missing

bio_df <- bio_df %>% 
  group_by(`individual id`) %>% 
  mutate(`days of age` = mean(`days of age`, na.rm = TRUE)) %>% 
  ungroup() %>% 
  replace_na(list(`days of age` = 0))
            
remove(non_unique_id)
# MERGE -------------------------------------------------------------------

bio_df %>% 
  anti_join(square_exp_df, by=c('sample id'='name')) %>% 
  select(`sample id`, outlier)
# all the missing samples in square_exp_df are outliers or mice, so godspeed!

merged_df <- square_exp_df %>% 
  inner_join(bio_df, by=c('name'='sample id'))
assert_true(dim(merged_df)[1] == dim(square_exp_df)[1])
remove(square_exp_df)
remove(bio_df)

assert(merged_df %>% 
  select(outlier) %>% 
  drop_na(.) %>% nrow(.) == 0)
merged_df <- merged_df %>% 
  select(-outlier)

all(merged_df[c("species.x","tissue.x")] == merged_df[c("species.y","tissue.y")]) 
#there are some discrepancies between .x and .y values!

merged_df %>%
  filter((species.x != species.y) | (tissue.x != tissue.y)) %>% 
  select(name, species.x, species.y, tissue.x, tissue.y)

#some names and species.y state rhesus. while species.x state macaque -> drop species.x
#and tissue are identical
merged_df <- merged_df %>% select(-species.x, -tissue.x)


merged_df %>% select(matches('^[^-]+$')) %>% colnames -> special_colnames

merged_df <- merged_df %>% relocate(all_of(special_colnames), 1) %>% rename(species=species.y, tissue=tissue.y)
remove(special_colnames)

merged_df <- clean_names(merged_df)
merged_df <- merged_df %>% 
  mutate(age_group = as.factor((years_of_age + days_of_age/365.25) %/% 10)) %>% 
  relocate(age_group, .after = gender) %>% 
  select(-c(years_of_age, days_of_age))

merged_df %>% write_rds('.cache/merged_df.rds')
special_names_table %>% write_rds('.cache/special_names_table.rds')
remove(list = c('merged_df', 'special_names_table'))
# First operations --------------------------------------------------------

shinyAppDir('violins/')
