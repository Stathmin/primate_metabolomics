
# imports -----------------------------------------------------------------
install.packages('tidyverse')
install.packages('readxl')
install.packages('janitor')
install.packages('purrr')
install.packages('stringi')

library('tidyverse')
library('readxl')
library('janitor')
library('purrr')
library('stringi')

# custom_fucntions --------------------------------------------------------

pivot_transpose <- function(
    data,
    cols_to_wide = 1,
    cols_to_long = 2:last_col(),
    names_to = "name",
    .drop = TRUE
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

# file_reading_expression -------------------------------------------------

expression_file <- './data/pbio.1001871.s022.xlsx'

exp_df <- readxl::read_xlsx(expression_file) # initial row contains bio_sample_data
exp_df %>% dim() # 4700  200 (!)

exp_df <- exp_df %>% pivot_transpose # transpose to make more logical indexing

#overview exp_df
exp_df %>% dim() # 199 4701 (!)

exp_df %>% slice(1:5) %>% 
  select(1:5) #upper-left chunk

exp_df %>% 
  slice(dim(.)[1] - (4:0)) %>% 
  select(dim(.)[2] - (4:0)) #lower_right chunk

#l_r chunk shows strange last row, needs checking
exp_df %>% 
  tail(1) %>% 
  select(\(x) ! is.na(x)) %>% 
  unlist(., use.names = FALSE) %>% 
  head() #contains alternative names for columns where applicable

# NO IDEA HOW TO WORK WITH ALT COL NAMES!!!
# so lets create separate table for them

special_names_table <- exp_df %>% #table of names for genes and loci
  tail(1) %>% 
  select(\(x) ! is.na(x)) %>% 
  pivot_longer(everything(),
               names_to = 'common_name',
               values_to = 'special_name') %>% 
  slice(2:nrow(.)) # because row 1 is c('name',  'annotation')


#now drop last row with non-numeric data
square_exp_df <- exp_df %>% 
  head(-1)

#functions to find suspicious values
.retrun_only_nonumeric_strings <- function(x){
  y <- as.character(x)
  if (stri_count_regex(y,"^[0-9\\.]+$")==1){
    return(NA_character_)
  }else{
    return(x)
  }
}

retrun_only_nonumeric_strings <- function(X){
  map(X, .retrun_only_nonumeric_strings) %>% unlist
}


#check for strange values in columns after first three
sus_values_mask <- square_exp_df %>% 
  select(-(1:3)) %>%
  mutate_if(is.character, ~replace_na(.,"")) %>% 
  mutate(across(everything(),retrun_only_nonumeric_strings))

sus_values_mask %>% unlist %>% unname %>% unique #ONLY NAs, great!

square_exp_df <- square_exp_df %>% mutate(across(colnames(.)[-(1:3)], as.double)) 

dim(square_exp_df) # 198 4701 no data lost, annotation removed (stored in special_names_table)
# sample_file_reading -----------------------------------------------------

bio_samples_file <- './data/pbio.1001871.s016.xlsx'

bio_df <- readxl::read_xlsx(bio_samples_file) # initial row contains bio_sample_data
bio_df %>% dim() # 280  8

#`years of age`, `days of age` are numeric

bio_df %>% 
  mutate(
    across(c(`years of age`, `days of age`), as.double)
    ) #creates some NAs, lets check why

bio_df %>% 
  select(c(`years of age`, `days of age`)) %>%
  mutate_if(is.character, ~replace_na(.,"")) %>% 
  mutate(across(everything(),retrun_only_nonumeric_strings)) %>% 
  unlist %>% unname %>% unique # has "NA", no problems with introduced NAs

bio_df <- bio_df %>% 
  mutate(
    across(c(`years of age`, `days of age`), as.double)
  )


# MERGE -------------------------------------------------------------------

square_exp_df %>% 
  anti_join(bio_df, by=c('name'='sample id')) # 0 × 4,701 YEAH, FULL MERGE!

merged_df <- square_exp_df %>% 
  left_join(bio_df, by=c('name'='sample id')) #198 × 4,708

all(merged_df[c("species.x","tissue.x")] == merged_df[c("species.y","tissue.y")]) #there are some discrepancies between .x and .y values!

merged_df %>%
  filter((species.x != species.y) | (tissue.x != tissue.y)) %>% 
  select(name, species.x, species.y, tissue.x, tissue.y)

#some names and species.y state rhesus. while species.x state macaque -> drop species.x
#and tissue are identical
merged_df <- merged_df %>% select(-species.x, -tissue.x)


merged_df %>% select(matches('^[^-]+$')) %>% colnames -> special_colnames

merged_df <- merged_df %>% relocate(special_colnames, 1) %>% rename(species=species.y, tissue=tissue.y)

merged_df <- clean_names(merged_df)


# First operations --------------------------------------------------------

merged_df %>% 
