library(tidyverse)
library(rmarkdown)
source("./mapping.R")

render_report = function(country, city) {
  rmarkdown::render(
    "input.Rmd", params = list(
      country = country,
      city = city
    ),
    output_file = paste0(country, city, ".pdf")
  )
}

country_list