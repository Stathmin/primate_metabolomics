
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

# file_reading_expression -------------------------------------------------

expression_file <- './data/pbio.1001871.s022.xlsx'

# sample_file_reading -----------------------------------------------------

bio_samples_file <- './data/pbio.1001871.s016.xlsx'

# MERGE -------------------------------------------------------------------
