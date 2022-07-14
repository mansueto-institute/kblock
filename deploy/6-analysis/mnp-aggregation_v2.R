

library(ggplot2)
library(sf)
library(tidyverse)
library(viridis)
library(patchwork)
library(scales)
library(tidygeocoder)
library(readxl)
library(writexl)
library(units)
library(ggsn)
library(arrow)
library(sfarrow)
library(osmdata)
library(viridis)
options(scipen = 999)

# UN population data ------------------------------------------------------

un_zip_url <- 'https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2022_Demographic_Indicators_Medium.zip'
tmp_filepath <- paste0(tempdir(), '/', basename(un_zip_url))
download.file(url = paste0(un_zip_url), destfile = tmp_filepath)
unzip(tmp_filepath, exdir=tempdir())
un_population <- read_csv(gsub(".zip", ".csv", tmp_filepath)) %>%
  select_all(~gsub("\\s+|\\.|\\/|,|\\*", "_", .)) %>%
  rename_all(list(tolower)) %>% 
  filter(time == 2020,
         loctypeid == 4) %>%
  select(iso3_code, location, tpopulation1jan) %>%
  rename(country_name = location, 
         poptotal = tpopulation1jan ) %>%
  mutate(poptotal = poptotal * 1000)
  
# # https://population.un.org/wpp/Download/Standard/CSV/
# location_url = 'https://population.un.org/wpp/Download/Files/4_Metadata/WPP2022_F01_LOCATIONS.XLSX'
# tmp_filepath <- paste0(tempdir(), '/', basename(location_url))
# download.file(url = paste0(location_url), destfile = tmp_filepath)
# location_codes <- read_xlsx(path = tmp_filepath, sheet = 'Location', range = "A17:H306")  %>%
#   select_all(~gsub("\\s+|\\.|\\/|,|\\*|-", "_", .)) %>%
#   rename_all(list(tolower)) %>%
#   filter(code == 4) %>%
#   rename(country_name = region__subregion__country_or_area_) %>%
#   select(location_code, country_name, iso3_alpha_code)

# un_population <- read_csv('https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_TotalPopulationBySex.csv')
# un_population <- un_population %>%
#   select_all(~gsub("\\s+|\\.|\\/|,|\\*", "_", .)) %>%
#   rename_all(list(tolower)) %>%
#   filter(time == 2020,  
#          variant == 'Medium') %>%
#   mutate(poptotal = poptotal * 1000,
#          popdensity = ((popdensity*.3861)/100)*1000) %>%
#   inner_join(., location_codes, by = c('locid' = 'location_code')) %>%
#   select(locid, country_name, iso3_alpha_code, location, time, variant, poptotal, popdensity) 

# Build dataset --------------------------------------------------------------

iso_code_list <- c( 'AGO', 'BDI', 'BEN', 'BFA', 'BWA', 'CAF', 'CIV', 'CMR', 'COD', 'COG', 'COM', 'CPV', 'DJI', 'ERI', 'ESH', 'ETH', 'GAB', 'GHA', 'GIN', 'GMB', 'GNB', 'GNQ', 'KEN', 'LBR', 'LSO', 'MDG', 'MLI', 'MOZ', 'MRT', 'MUS', 'MWI', 'NAM', 'NER', 'NGA', 'RWA', 'SDN', 'SEN', 'SLE', 'SOM', 'SSD', 'STP', 'SWZ', 'SYC', 'TCD', 'TGO', 'TZA', 'UGA', 'ZAF', 'ZMB', 'ZWE')        

df_combined <- purrr::map_dfr(.x = iso_code_list, .f = function(i) {
  print(i)
  df_ghsl <- read_parquet(list.files(paste0('/Users/nm/Downloads/outputs/crosswalk/partitions/country_code=',i), full.names = TRUE))
  df_complexity <- read_parquet(paste0('/Users/nm/Downloads/outputs/complexity/complexity_',i,'.parquet')) %>%
    select(block_id, gadm_code, country_code, block_area, on_network_street_length, off_network_street_length, nearest_external_street, building_area, building_count, building_layers, k_complexity)   
  df_population <- read_parquet(paste0('/Users/nm/Downloads/outputs/population/population_',i,'.parquet')) %>% 
    select(block_id, landscan_population, worldpop_population)
  
  # manual correction
  df_ghsl <- df_ghsl %>% 
    mutate(urban_id = case_when(urban_id == '2544.0' ~ "2565.0", # Combine Abuja Nigeria IDs (split for unknown reasons)
                                TRUE ~ as.character(urban_id)),
           area_type_4 = case_when(urban_id == "2565.0" ~ 'Core urban', # Fix urban hierarchy codings
                                       TRUE ~ as.character(area_type_4)))
  
  df_full <- df_complexity %>% 
    left_join(x =., y = df_population, by = c('block_id' = 'block_id')) %>%
    left_join(x =., y = df_ghsl, by = c('block_id' = 'block_id'))  %>%
    left_join(., un_population %>% select(iso3_code, poptotal), 
              by = c('country_code' = 'iso3_code')) 
  
  df_full <- df_full %>%
    group_by(country_code) %>%
    mutate(landscan_population_total = sum(landscan_population),
           worldpop_population_total = sum(worldpop_population)) %>%
    ungroup() %>%
    mutate(landscan_population_un = round(poptotal * landscan_population/landscan_population_total,0),
           worldpop_population_un = round(poptotal * worldpop_population/worldpop_population_total,0)) %>%
    select(-one_of(c('landscan_population_total','worldpop_population_total','poptotal'))) %>%
    mutate(block_count = 1) %>% 
    mutate(on_network_street_length = replace_na(on_network_street_length, 0), 
           off_network_street_length = replace_na(off_network_street_length, 0)) %>%
    mutate(k_5 = case_when(nearest_external_street > 0 ~ 'Off network', 
                           k_complexity >= 31 ~ "31+",
                           TRUE ~ as.character(k_complexity)),
           k_4 = case_when(nearest_external_street > 0 ~ 'Off network', 
                           k_complexity >= 11 & k_complexity <= 15 ~ "11 to 15", 
                           k_complexity >= 16 & k_complexity <= 20 ~ "16 to 20", 
                           k_complexity >= 21 ~ "21+",
                           TRUE ~ as.character(k_complexity)),
           k_3 = case_when(nearest_external_street > 0 ~ 'Off network', 
                           k_complexity >= 1 & k_complexity <= 3 ~ "1 to 3",
                           k_complexity >= 4 & k_complexity <= 10 ~ "4 to 10",
                           k_complexity >= 11 ~ "11+"),
           k_2 = case_when(nearest_external_street > 0 ~ 'Off network', 
                           k_complexity >= 1 & k_complexity <= 3 ~ "1 to 3",
                           TRUE ~ as.character('4+')),
           k_1 = case_when(nearest_external_street > 0 ~ 'Off network', 
                           TRUE ~ as.character('Total')),
           k_0 = 'Total') %>%
    group_by(country_name, country_code, area_type_2, area_type_3, area_type_4, urban_id, urban_center_name, urban_country_name, conurbation_id, conurbation_area_name, conurbation_country_name, conurbation_area_name_short, k_complexity, k_5, k_4, k_3, k_2, k_1, k_0) %>%
    summarize_at(vars(block_count, block_area, on_network_street_length, off_network_street_length, nearest_external_street, building_area, building_count, landscan_population, worldpop_population, landscan_population_un, worldpop_population_un), list(sum)) %>%
    ungroup() 
})

# Create labels -----------------------------------------------------------

k_5_order = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31+','Off network')        
#k_4_order = c("1","2","3","4","5","6","7","8","9","10","11","21+","Off network")
k_4_order = c("1","2","3","4","5","6","7","8","9","10","11 to 15","16 to 20","21+","Off network")
k_3_order = c("1 to 3",'4 to 10','11+','Off network')
k_2_order = c("1 to 3","4+",'Off network')
k_1_order = c("Total",'Off network')

df_combined_prep <- df_combined %>%
  mutate(k_5 = factor(k_5, levels = k_5_order),
         k_4 = factor(k_4, levels = k_4_order),
         k_3 = factor(k_3, levels = k_3_order),
         k_2 = factor(k_2, levels = k_2_order),
         k_1 = factor(k_1, levels = k_1_order)) %>%
  mutate(continent = 'sub-Saharan Africa',
         country_name_area_type_4 = paste0(country_name,' (', tolower(area_type_4),')'),
         country_name_area_type_3 = paste0(country_name,' (', tolower(area_type_3),')'),
         country_name_area_type_2 = paste0(country_name,' (', tolower(area_type_2),')'),
         urban_center_name_area_type_4 = paste0(urban_center_name,' (', tolower(area_type_4),')'),
         urban_center_name_area_type_3 = paste0(urban_center_name,' (', tolower(area_type_3),')'),
         urban_center_name_area_type_2 = paste0(urban_center_name,' (', tolower(area_type_2),')'),
         conurbation_area_name_short_area_type_4 = paste0(conurbation_area_name_short,' (', tolower(area_type_4),')'),
         conurbation_area_name_short_area_type_3 = paste0(conurbation_area_name_short,' (', tolower(area_type_3),')'),
         conurbation_area_name_short_area_type_2 = paste0(conurbation_area_name_short,' (', tolower(area_type_2),')') ) 

# Crosstab geographies ----------------------------------------------------

geo_labels <- df_combined_prep %>%
  select(continent, area_type_2, area_type_3, area_type_4, 
         country_name, country_name_area_type_2, country_name_area_type_3, country_name_area_type_4, 
         urban_center_name, urban_country_name, urban_center_name_area_type_2, urban_center_name_area_type_3, urban_center_name_area_type_4, 
         conurbation_area_name_short, conurbation_country_name, conurbation_area_name_short_area_type_2, conurbation_area_name_short_area_type_3, conurbation_area_name_short_area_type_4, urban_id, conurbation_id) %>% 
  distinct()

geo_xwalk <- bind_rows(geo_labels %>% select(continent) %>% distinct() %>% 
                         mutate(group_var = '1 - sub-Saharan Africa', group_val = continent) %>% select(group_var),
                       
                       geo_labels %>% select(area_type_2) %>% distinct() %>% 
                         mutate(group_var = '1a - Urban / non-urban', group_val = area_type_2) ,
                       geo_labels %>% select(area_type_3, area_type_2) %>% distinct() %>% 
                         mutate(group_var = '1b - Urban / peri-urban / non-urban', group_val = area_type_3) ,
                       geo_labels %>% select(area_type_4, area_type_3, area_type_2) %>% distinct() %>% 
                         mutate(group_var = '1c - Urban hierarchy', group_val = area_type_4),
                       
                       geo_labels %>% select(country_name) %>% distinct() %>% 
                         mutate(group_var = '2 - Country', group_val = country_name),
                       geo_labels %>% select(country_name_area_type_2, country_name, area_type_2) %>% distinct() %>% 
                         mutate(group_var = '2a - Country by urban / non-urban') %>% rename(group_val = country_name_area_type_2) ,
                       geo_labels %>% select(country_name_area_type_3, country_name, area_type_3, area_type_2) %>% distinct() %>% 
                         mutate(group_var = '2b - Country by urban / peri-urban / nonurban') %>% rename(group_val = country_name_area_type_3),
                       geo_labels %>% select(country_name_area_type_4, country_name, area_type_4, area_type_3, area_type_2) %>% distinct() %>% 
                         mutate(group_var = '2c - Country by urban hierarchy') %>% rename(group_val = country_name_area_type_4),
                       
                       geo_labels %>% select(urban_center_name, urban_country_name, conurbation_area_name_short, urban_id, conurbation_id, area_type_4, area_type_3, area_type_2) %>% rename(country_name = urban_country_name) %>% distinct() %>% 
                         mutate(group_var = '3 - Urban center', group_val = urban_center_name),
                       #geo_labels %>% select(urban_center_name_area_type_2, urban_center_name, urban_country_name, conurbation_area_name_short, area_type_2) %>% rename(country_name = urban_country_name) %>% distinct() %>% 
                       #  mutate(group_var = '4a - Urban center by urban / non-urban') %>% rename(group_val = urban_center_name_area_type_2),
                       geo_labels %>% select(urban_center_name_area_type_3, urban_center_name, urban_country_name, conurbation_area_name_short, urban_id, conurbation_id, area_type_3, area_type_2) %>% rename(country_name = urban_country_name) %>% distinct() %>% 
                         mutate(group_var = '3b - Urban center by urban / peri-urban / nonurban') %>% rename(group_val = urban_center_name_area_type_3),
                       geo_labels %>% select(urban_center_name_area_type_4, urban_center_name, urban_country_name, conurbation_area_name_short, urban_id, conurbation_id, area_type_4, area_type_3, area_type_2) %>% rename(country_name = urban_country_name) %>% distinct() %>% 
                         mutate(group_var = '3c - Urban center by urban hierarchy') %>% rename(group_val = urban_center_name_area_type_4),
                       
                       geo_labels %>% select(conurbation_area_name_short, conurbation_country_name, area_type_2) %>% rename(country_name = conurbation_country_name) %>% distinct() %>%
                         mutate(group_var = '4 - Conurbation', group_val = conurbation_area_name_short),
                       geo_labels %>% select(conurbation_area_name_short_area_type_2, conurbation_area_name_short, conurbation_country_name, conurbation_id, area_type_2) %>% rename(country_name = conurbation_country_name) %>% distinct() %>%
                         mutate(group_var = '4a - Conurbation by urban / non-urban') %>% rename(group_val = conurbation_area_name_short_area_type_2),
                       geo_labels %>% select(conurbation_area_name_short_area_type_3, conurbation_area_name_short, conurbation_country_name, conurbation_id, area_type_3, area_type_2) %>% rename(country_name = conurbation_country_name) %>% distinct() %>%
                         mutate(group_var = '4b - Conurbation by urban / peri-urban / nonurban') %>% rename(group_val = conurbation_area_name_short_area_type_3),
                       geo_labels %>% select(conurbation_area_name_short_area_type_4, conurbation_area_name_short, conurbation_country_name, conurbation_id, area_type_4, area_type_3, area_type_2) %>% rename(country_name = conurbation_country_name) %>% distinct() %>%
                         mutate(group_var = '4c - Conurbation by urban hierarchy') %>% rename(group_val = conurbation_area_name_short_area_type_4)
                       ) %>%
  relocate(group_var, group_val, country_name, area_type_2, area_type_3, area_type_4, urban_center_name, conurbation_area_name_short) %>%
  rename(urban_nonurban = area_type_2,
         urban_periurban_nonurban = area_type_3,
         urban_hierarchy = area_type_4)

# Crosstab aggregation function -------------------------------------------

recode_list = list('continent' = '1 - sub-Saharan Africa',
                   'area_type_2' = '1a - Urban / non-urban',
                   'area_type_3' = '1b - Urban / peri-urban / non-urban',
                   'area_type_4' = '1c - Urban hierarchy',
                   'country_name' = '2 - Country',
                   'country_name_area_type_2' = '2a - Country by urban / non-urban',
                   'country_name_area_type_3' = '2b - Country by urban / peri-urban / nonurban',
                   'country_name_area_type_4' = '2c - Country by urban hierarchy',
                   'urban_center_name' = '3 - Urban center',
                   'urban_center_name_area_type_3' = '3b - Urban center by urban / peri-urban / nonurban',
                   'urban_center_name_area_type_4' = '3c - Urban center by urban hierarchy',
                   'conurbation_area_name_short' = '4 - Conurbation',
                   'conurbation_area_name_short_area_type_2' = '4a - Conurbation by urban / non-urban',
                   'conurbation_area_name_short_area_type_3' = '4b - Conurbation by urban / peri-urban / nonurban',
                   'conurbation_area_name_short_area_type_4' = '4c - Conurbation by urban hierarchy')


process_blocks <- function(data, crosswalk_data, pivot_col, crosstab_list) { 
  
  data <- data %>% 
    mutate(across(c(block_count, block_area, on_network_street_length, off_network_street_length, nearest_external_street, building_area, building_count, landscan_population_un, worldpop_population_un), ~replace_na(.x, 0))) %>%
    mutate(nearest_external_street_block_count = case_when(nearest_external_street > 0 ~ 1, TRUE ~ as.double(0))) %>%
    mutate(k_complexity_weighted_landscan = k_complexity * landscan_population_un, 
           k_complexity_weighted_worldpop = k_complexity * worldpop_population_un) %>%
    group_by_at(vars(all_of(c(pivot_col, crosstab_list)))) %>%
    summarize_at(vars(k_complexity_weighted_landscan, k_complexity_weighted_worldpop, block_count, block_area, on_network_street_length, off_network_street_length, nearest_external_street, nearest_external_street_block_count, building_area, building_count, landscan_population_un, worldpop_population_un), list(sum)) %>%
    ungroup() 
    
  data <- purrr::map_dfr(.x = crosstab_list, .f = ~ data %>%
                           group_by_at(vars(.x, all_of(pivot_col) )) %>%
                           summarize_at(vars(k_complexity_weighted_landscan, k_complexity_weighted_worldpop, block_count, block_area, on_network_street_length, off_network_street_length, nearest_external_street, nearest_external_street_block_count, building_area, building_count, landscan_population_un, worldpop_population_un), 
                                        list(sum)) %>% ungroup() %>%
                           mutate(group_var = paste0(.x)) %>%
                           rename(group_val = all_of(.x))
                         ) %>%
    relocate(group_var, group_val) %>%
    mutate(group_var = case_when(group_var == 'continent' ~ '1 - sub-Saharan Africa',
                                 group_var == 'area_type_2' ~ '1a - Urban / non-urban',
                                 group_var == 'area_type_3' ~ '1b - Urban / peri-urban / non-urban',
                                 group_var == 'area_type_4' ~ '1c - Urban hierarchy',
                                 group_var == 'country_name' ~ '2 - Country',
                                 group_var == 'country_name_area_type_2' ~ '2a - Country by urban / non-urban',
                                 group_var == 'country_name_area_type_3' ~ '2b - Country by urban / peri-urban / nonurban',
                                 group_var == 'country_name_area_type_4' ~ '2c - Country by urban hierarchy',
                                 group_var == 'urban_center_name' ~ '3 - Urban center',
                                 group_var == 'urban_center_name_area_type_3' ~ '3b - Urban center by urban / peri-urban / nonurban',
                                 group_var == 'urban_center_name_area_type_4' ~ '3c - Urban center by urban hierarchy',
                                 group_var == 'conurbation_area_name_short' ~ '4 - Conurbation',
                                 group_var == 'conurbation_area_name_short_area_type_2' ~ '4a - Conurbation by urban / non-urban',
                                 group_var == 'conurbation_area_name_short_area_type_3' ~ '4b - Conurbation by urban / peri-urban / nonurban',
                                 group_var == 'conurbation_area_name_short_area_type_4' ~ '4c - Conurbation by urban hierarchy')) %>%
    arrange(group_var, group_val, !!!all_of(pivot_col)) %>%
    mutate(block_area_hectares = na_if((block_area*0.0001), 0)) %>%
    replace_na(list(block_hectares = 0)) %>% 
    mutate(landscan_pop_density_hectare = landscan_population_un/block_area_hectares,
           worldpop_pop_density_hectare = worldpop_population_un/block_area_hectares,
           k_complexity_weighted_landscan = k_complexity_weighted_landscan/landscan_population_un, 
           k_complexity_weighted_worldpop = k_complexity_weighted_worldpop/worldpop_population_un,
           landscan_population_per_building_area = landscan_population_un/building_area,
           worldpop_population_per_building_area = worldpop_population_un/building_area,
           landscan_population_per_building = landscan_population_un/building_count,
           worldpop_population_per_building = worldpop_population_un/building_count,
           buildings_per_block = building_count/block_count,
           average_on_network_street_length_per_block = on_network_street_length/building_count,
           buildings_per_hectare = building_count/block_area_hectares,
           average_building_area = building_area/building_count) %>%
    group_by(group_var, group_val) %>%
    mutate(sum_block_count = sum(block_count),
           sum_block_area = sum(block_area),
           sum_on_network_street_length = sum(on_network_street_length),
           sum_off_network_street_length = sum(off_network_street_length),
           sum_nearest_external_street = sum(nearest_external_street),
           sum_building_area = sum(building_area),
           sum_building_count = sum(building_count),
           sum_landscan_population_un = sum(landscan_population_un),
           sum_worldpop_population_un = sum(worldpop_population_un)) %>%
    ungroup() %>%
    mutate(share_block_count = block_count  / sum_block_count,
           share_block_area = block_area / sum_block_area,
           share_on_network_street_length = on_network_street_length / sum_on_network_street_length,
           share_off_network_street_length = off_network_street_length / sum_off_network_street_length,
           share_nearest_external_street = nearest_external_street / sum_nearest_external_street,
           share_building_area = building_area / sum_building_area,
           share_building_count = building_count / sum_building_count,
           share_landscan_population_un = landscan_population_un / sum_landscan_population_un,
           share_worldpop_population_un = worldpop_population_un / sum_worldpop_population_un) %>%
    left_join(., crosswalk_data, by = c('group_var'='group_var', 'group_val'='group_val')) %>%
    mutate(geographic_level = case_when(group_var %in% c('1 - sub-Saharan Africa', '1a - Urban / non-urban', '1b - Urban / peri-urban / non-urban', '1c - Urban hierarchy') ~ '1 - sub-Saharan Africa',
                                        group_var %in% c('2 - Country', '2a - Country by urban / non-urban', '2b - Country by urban / peri-urban / nonurban', '2c - Country by urban hierarchy') ~ '2 - Country',
                                        group_var %in% c('3 - Urban center', '3b - Urban center by urban / peri-urban / nonurban', '3c - Urban center by urban hierarchy') ~ '3 - Urban center',
                                        group_var %in%  c('4 - Conurbation', '4a - Conurbation by urban / non-urban', '4b - Conurbation by urban / peri-urban / nonurban', '4c - Conurbation by urban hierarchy') ~ '4 - Conurbation')) %>%
    relocate(any_of(c('geographic_level','group_var', 'group_val','country_name', 'urban_nonurban', 'urban_periurban_nonurban', 'urban_hierarchy', 'urban_center_name', 'conurbation_area_name_short', pivot_col, 'k_complexity_weighted_landscan', 'landscan_population_un', 'landscan_pop_density_hectare', 'landscan_population_per_building_area', 'landscan_population_per_building', 'block_count', 'block_area', 'on_network_street_length', 'off_network_street_length', 'nearest_external_street', 'nearest_external_street_block_count', 'building_area', 'building_count', 'block_area_hectares', 'buildings_per_block', 'average_on_network_street_length_per_block', 'buildings_per_hectare', 'average_building_area', 'share_landscan_population_un', 'share_block_count', 'share_block_area', 'share_on_network_street_length', 'share_off_network_street_length', 'share_nearest_external_street', 'share_building_area', 'share_building_count', 'k_complexity_weighted_worldpop', 'worldpop_population_un', 'worldpop_pop_density_hectare', 'worldpop_population_per_building_area', 'worldpop_population_per_building', 'share_worldpop_population_un', 'urban_id', 'conurbation_id', 'above_500k_population'))) %>%
    rename(geography_group = group_var,
           geography_area = group_val) %>%
    group_by(geography_group, country_name, urban_center_name, conurbation_area_name_short) %>%
    mutate(area_landscan_population_un = sum(landscan_population_un),
           area_worldpop_population_un = sum(worldpop_population_un)) %>%
    ungroup() %>%
    arrange(geography_group, desc(area_landscan_population_un), !!!all_of(pivot_col), geography_area) %>%
    select(-one_of(c('sum_block_count', 'sum_block_area', 'sum_on_network_street_length', 'sum_off_network_street_length', 'sum_nearest_external_street', 'sum_building_area', 'sum_building_count')))  

  return(data)
}

# Run crosstabs ---------------------------------------------------------------


crosstab_list_input = c('continent', 'area_type_2', 'area_type_3', 'area_type_4', 'country_name', 'country_name_area_type_2', 'country_name_area_type_3', 'country_name_area_type_4', 'urban_center_name', 'urban_center_name_area_type_3', 'urban_center_name_area_type_4', 'conurbation_area_name_short', 'conurbation_area_name_short_area_type_2', 'conurbation_area_name_short_area_type_3', 'conurbation_area_name_short_area_type_4')

out_k_0 <- process_blocks(data = df_combined_prep, crosswalk_data = geo_xwalk, pivot_col = 'k_0', crosstab_list = crosstab_list_input) 
out_k_1 <- process_blocks(data = df_combined_prep, crosswalk_data = geo_xwalk, pivot_col = 'k_1', crosstab_list = crosstab_list_input) 
out_k_2 <- process_blocks(data = df_combined_prep, crosswalk_data = geo_xwalk, pivot_col = 'k_2', crosstab_list = crosstab_list_input) 
out_k_3 <- process_blocks(data = df_combined_prep, crosswalk_data = geo_xwalk, pivot_col = 'k_3', crosstab_list = crosstab_list_input) 
out_k_4 <- process_blocks(data = df_combined_prep, crosswalk_data = geo_xwalk, pivot_col = 'k_4', crosstab_list = crosstab_list_input) 
out_k_5 <- process_blocks(data = df_combined_prep, crosswalk_data = geo_xwalk, pivot_col = 'k_5', crosstab_list = crosstab_list_input) 


file_labels = data.frame(
  column_names = c('geographic_level', 'geography_group', 'geography_area', 'country_name', 'urban_nonurban', 'urban_periurban_nonurban', 'urban_hierarchy', 'urban_center_name', 'conurbation_area_name', 'k_', 'k_complexity_weighted_landscan', 'landscan_population_un', 'landscan_pop_density_hectare', 'landscan_population_per_building_area', 'landscan_population_per_building', 'block_count', 'block_area', 'on_network_street_length', 'off_network_street_length', 'nearest_external_street', 'nearest_external_street_block_count', 'building_area', 'building_count', 'block_area_hectares', 'buildings_per_block', 'average_on_network_street_length_per_block', 'buildings_per_hectare', 'average_building_area', 'share_landscan_population_un', 'share_block_count', 'share_block_area', 'share_on_network_street_length', 'share_off_network_street_length', 'share_nearest_external_street', 'share_building_area', 'share_building_count', 'k_complexity_weighted_worldpop', 'worldpop_population_un', 'worldpop_pop_density_hectare', 'worldpop_population_per_building_area', 'worldpop_population_per_building', 'share_worldpop_population_un', 'urban_id', 'conurbation_id', 'sum_landscan_population_un', 'sum_worldpop_population_un', 'area_landscan_population_un', 'area_worldpop_population_un'),
  column_labels = c('Geography level', 'Geography grouping', 'Geography area', 'Country', 'Urban or non-urban', 'Urban, peri-urban, or non-urban', 'Urban hierarchy', 'Urban center', 'Conurbation', 'K value', 'K complexity (weighted by population)', 'Population', 'Population per hectare', 'Population per building area', 'Population per building', 'Block count', 'Block area', 'On network street length', 'Off network street length', 'Nearest external street', 'Nearest external street block count', 'Building area', 'Building count', 'Block area (hectares)', 'Buildings per block', 'Average on network street length', 'Buildings per hectare', 'Average building area', 'Share of population', 'Share of block count', 'Share of block area', 'Share of on network street length', 'Share of off network street length', 'Share of nearest external street', 'Share of building area', 'Share of building count', 'K complexity (weighted by population)', 'Population', 'Population per hectare', 'Population per building area', 'Population per building', 'Share of population', 'Urban ID', 'Conurbation ID', 'Sum of population by geography area', 'Sum of population by geography area', 'Sum of population by county or urban / conurbation area', 'Sum of population by county or urban / conurbation area')
)

write_xlsx(list("k3" = out_k_3, "k2" = out_k_2, "k1" = out_k_1, "k0" = out_k_0, 'labels' = file_labels), #"k_complexity" = out_k_complexity, 
           col_names = TRUE, format_headers = TRUE, path = '/Users/nm/Desktop/mnp-analysis.xlsx')


# Scatter / dot plots for sub-Saharan Africa ----------------------------------------

area_4 = c("Core urban","Peripheral urban","Peri-urban","Non-urban")

africa_dviz <- out_k_0 %>% filter(geography_group == "2c - Country by urban hierarchy") %>%
  mutate(urban_hierarchy = factor(urban_hierarchy, levels = area_4)) %>% 
  mutate(k_reweight = landscan_population_un*k_complexity_weighted_landscan) %>%
  group_by(country_name) %>%
  mutate(k_reweight = sum(k_reweight)) %>%
  ungroup() %>%
  mutate(k_reweight = k_reweight/area_landscan_population_un)
  
# X population Size k-complexity for sub-Saharan Africa
(africa_dots_xpop_sizek <- ggplot() +
  geom_point(data = africa_dviz %>% filter(landscan_population_un >= 10000), 
             aes(x = reorder(country_name, (area_landscan_population_un)), y = log10(landscan_population_un),
                 fill = urban_hierarchy, color = urban_hierarchy, size = k_complexity_weighted_landscan), alpha = .75) +
  scale_size(range = c(1,10)) + 
  geom_text(data = africa_dviz %>% filter(landscan_population_un >= 10000), 
            aes(x =  reorder(country_name, (area_landscan_population_un)), y = log10(landscan_population_un), label = round(k_complexity_weighted_landscan , 1)), size = 3,
            vjust = .5, color = '#333333', fontface='bold', check_overlap = TRUE) +
  scale_y_continuous(oob = scales::squish, breaks= c(1,2,3,4,5,6,7,8,9), labels = c('0',"100","1K","10K","100K","1M","10M","100M","1B"))+
  scale_fill_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B'), guide = 'none') +
  scale_color_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B'), guide = 'none') +
  labs(x= '', y = 'Population', size = 'Average k-complexity',
        fill = 'Urban hierarchy', color =  'Urban hierarchy') +
  #guides(color = guide_legend(override.aes = list(size = 6))) +
  coord_flip() +
  theme(legend.position = 'bottom',
        legend.key=element_blank()))

# X k-complexity Size population for sub-Saharan Africa
(africa_dots_xk_sizepop <- ggplot() +
  geom_point(data = africa_dviz, 
             aes(x =  reorder(country_name, (k_reweight)), y = k_complexity_weighted_landscan,
                 fill = urban_hierarchy, color = urban_hierarchy, size = landscan_population_un), alpha = .75) +
  geom_text(data = africa_dviz %>% filter(landscan_population_un >= 10000), 
            aes(x =  reorder(country_name, (k_reweight)), y =  k_complexity_weighted_landscan,
                label = ifelse(landscan_population_un >= 500000, paste0(round(landscan_population_un/1000000,1),'M'),'' ) ), 
            size = 3,
            vjust = .5, color = '#333333', fontface='bold', check_overlap = TRUE) +
  labs(x= '', y = 'Average k-complexity', size = 'Population',
       fill = 'Urban hierarchy', color =  'Urban hierarchy') +
  coord_flip() +
  scale_fill_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
  scale_color_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
  scale_size_continuous(range = c(1,10), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
  #scale_y_continuous( expand = c(.01,.01), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
  guides(color = guide_legend(override.aes = list(size = 8, alpha = 1))) +
  theme(legend.position = 'bottom',
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.key=element_blank()))  
  
# Y k-complexity X population for sub-Saharan Africa
(africa_x_y_scatter <- ggplot() +
  geom_point(data = africa_dviz %>% filter(landscan_population_un >= 10000), 
             aes(x = k_complexity_weighted_landscan, y = log10(landscan_population_un),
                 fill = urban_hierarchy, color = urban_hierarchy,
                 size = k_complexity_weighted_landscan), alpha = .7) +
  coord_flip() +
  labs( x= 'Average k-complexity', y = 'Population', size = 'k-complexity',
        fill = 'Urban hierarchy', color =  'Urban hierarchy') +
  #scale_y_continuous( expand = c(.005,.005), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) + 
  scale_fill_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B'), guide = 'none') + # 
  scale_color_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B'), guide = 'none') + # 
  scale_y_continuous(oob = scales::squish, breaks= c(1,2,3,4,5,6,7,8,9), labels = c('0',"100","1K","10K","100K","1M","10M","100M","1B"))+
  geom_text(data = africa_dviz %>% filter(landscan_population_un >= 10000), aes(x = k_complexity_weighted_landscan, y = log10(landscan_population_un),
                label = country_name), check_overlap = TRUE,
            size = 2.5, vjust =.5, color = '#333333', fontface='bold') +
  #guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.position = 'bottom',
        legend.key=element_blank()))

# Histogram for sub-Saharan Africa ----------------------------------------

k_5_order = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31+','Off\nnetwork')        
africa_hist <- out_k_5 %>% filter(geography_group == '1c - Urban hierarchy') %>%
  mutate(k_5 = str_wrap(k_5,width=3),
         k_5 = factor(k_5, levels = k_5_order)) %>%
  group_by(k_5) %>%
  mutate(k_sum = sum(landscan_population_un)) %>%
  ungroup() %>%
  mutate(share = landscan_population_un / k_sum ) %>%
  mutate(urban_hierarchy = factor(urban_hierarchy, levels = area_4)) %>%
  arrange(factor(urban_hierarchy, levels = rev(area_4)), k_5) %>%
  group_by(k_5) %>%
  mutate(pos_id_val = (cumsum(landscan_population_un) - 0.5*landscan_population_un)) %>%
  ungroup()

(africa_hist_k_5 <- ggplot() +
  geom_bar(data = africa_hist , 
           aes(y = k_5, x = landscan_population_un, fill = urban_hierarchy), color = 'white', size = .3, stat="identity") + 
    coord_flip() + 
  labs(y= 'k-complexity', x = 'Population', fill = 'Urban hierarchy', color =  'Urban hierarchy') + 
  geom_text(data = africa_hist , aes(y = k_5, x = pos_id_val, label = ifelse(landscan_population_un >= 5000000, paste0(round(share*100,0),"%"),'')),
            size = 3, vjust = .5, color = '#333333', fontface='bold') +
  scale_fill_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
  scale_x_continuous(expand = c(.005,.005), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
  theme(legend.position = 'bottom'))
  

# Conurbations ------------------------------------------------------------

unique(out_k_0$geography_group)
area_4 = c("Core urban","Peripheral urban","Peri-urban","Non-urban")

conurbation_dviz <- out_k_0 %>% filter(geography_group == "4c - Conurbation by urban hierarchy") %>%
  filter(urban_periurban_nonurban != 'Non-urban') %>%
  mutate(country_name =  (gsub('Democratic Republic of the Congo', 'DR Congo', as.character(country_name))),
         conurbation_area_name_short = gsub('Ambatolampy Tsimahafotsy', 'Ambatolampy', as.character(conurbation_area_name_short))
         ) %>% 
  mutate(conurbation_name = paste0(conurbation_area_name_short,', ',country_name)) %>%
  mutate(conurbation_first_name = conurbation_area_name_short) %>%
  separate(col =conurbation_first_name , sep = ' — ', into = c('conurbation_first_name'), extra = 'drop') %>%
  filter(area_landscan_population_un >= 1500000) %>%
  mutate(urban_hierarchy = factor(urban_hierarchy, levels = area_4)) %>% 
  mutate(k_reweight = landscan_population_un*k_complexity_weighted_landscan) %>%
  group_by(country_name) %>%
  mutate(k_reweight = sum(k_reweight)) %>%
  ungroup() %>%
  mutate(k_reweight = k_reweight/area_landscan_population_un) 

(conurbation_dots_xpop_sizek <- ggplot() +
    geom_point(data = conurbation_dviz %>% filter(landscan_population_un >= 10000), 
               aes(x =  reorder(conurbation_name, (area_landscan_population_un)), y = log10(landscan_population_un),
                   fill = urban_hierarchy, color = urban_hierarchy, size = k_complexity_weighted_landscan), alpha = .75) +
    scale_size(range = c(1,10)) + 
    geom_text(data = conurbation_dviz %>% filter(landscan_population_un >= 10000) , 
              aes(x =  reorder(conurbation_name, (area_landscan_population_un)), y = log10(landscan_population_un), label = round(k_complexity_weighted_landscan , 1)), size = 3,
              vjust = .5, color = '#333333', fontface='bold', check_overlap = TRUE) +
    scale_y_continuous(oob = scales::squish, breaks= c(1,2,3,4,5,6,7,8,9), labels = c('0',"100","1K","10K","100K","1M","10M","100M","1B"))+
    scale_fill_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
    scale_color_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
    labs(x= '', y = 'Population', size = 'Average K-complexity',
         fill = 'Urban hierarchy', color =  'Urban hierarchy') +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    coord_flip() +
    theme(legend.position = 'bottom',
          legend.key=element_blank()))

(conurbation_dots_xk_sizepop <- ggplot() +
    geom_point(data = conurbation_dviz, 
               aes(x =  reorder(conurbation_name, (k_reweight)), y = k_complexity_weighted_landscan,
                   fill = urban_hierarchy, color = urban_hierarchy,
                   size = landscan_population_un
               ), alpha = .75) +
    geom_text(data = conurbation_dviz, 
              aes(x =  reorder(conurbation_name, (k_reweight)), y =  k_complexity_weighted_landscan,
                  label = ifelse(landscan_population_un >= 100000, paste0(round(landscan_population_un/1000000,1),'M'),'' ) ), 
              size = 3,
              vjust = .5, color = '#333333', fontface='bold', check_overlap = TRUE) +
    labs(x= '', y = 'Average K-complexity', size = 'Population',
         fill = 'Urban hierarchy', color =  'Urban hierarchy') +
    coord_flip() +
    scale_fill_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
    scale_color_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
    scale_size_continuous(range = c(1,10), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
    #scale_y_continuous( expand = c(.01,.01), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    theme(legend.position = 'bottom',
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          legend.key=element_blank()))

(conurbation_x_y_scatter <- ggplot() +
    geom_point(data = conurbation_dviz %>% filter(landscan_population_un >= 10000), 
               aes(x = k_complexity_weighted_landscan, y = log10(landscan_population_un),
                   fill = urban_hierarchy, color = urban_hierarchy,
                   size = k_complexity_weighted_landscan), alpha = .7) +
    coord_flip() +
    labs( x= 'Average K-complexity', y = 'Population', size = 'k-complexity',
          fill = 'Urban hierarchy', color =  'Urban hierarchy') +
    #scale_y_continuous( expand = c(.005,.005), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) + 
    scale_fill_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
    scale_color_manual(values = c('#4D96FF','#6BCB77','#FFD93D','#FF6B6B')) +
    scale_y_continuous(oob = scales::squish, breaks= c(1,2,3,4,5,6,7,8,9), labels = c('0',"100","1K","10K","100K","1M","10M","100M","1B"))+
    geom_text(data = conurbation_dviz %>% filter(landscan_population_un >= 10000), 
              aes(x = k_complexity_weighted_landscan, y = log10(landscan_population_un),
                  label = conurbation_first_name), check_overlap = TRUE,
              size = 3, vjust =.5, color = '#333333', fontface='bold') +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    theme(legend.position = 'bottom',
          legend.key=element_blank()))

# Urban center bars -------------------------------------------------------

urban_bchart <- out_k_4 %>%
  filter(geography_group == '3 - Urban center',
         urban_hierarchy %in% c('Core urban','Peripheral urban'),
         area_landscan_population_un >= 1000000) %>%
  mutate(country_name =  (gsub('Democratic Republic of the Congo ', 'DR Congo', as.character(country_name)))) %>%
  mutate(urban_name = paste0(urban_center_name,', ',country_name)) %>%
  mutate(k_reweight = landscan_population_un*k_complexity_weighted_landscan) %>%
  group_by(geography_area) %>%
  mutate(k_reweight = sum(k_reweight)) %>%
  ungroup() %>%
  mutate(k_reweight = k_reweight/area_landscan_population_un) 

grey2 <- c('#414141','#777777')
kdist = max(as.integer(urban_bchart$k_4))
colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(urban_bchart$k_4))-2)

(a <- ggplot() +
  geom_bar(data = urban_bchart, aes(y = landscan_population_un, x = reorder(urban_name, area_landscan_population_un), 
                              fill = k_4), 
           stat="identity") +
  coord_flip() +
  scale_fill_manual(values = c(grey2, colorhexes)) +
  scale_y_continuous( expand = c(0, 0), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
  theme_bw() + 
  labs(y = 'Population', x = 'k-complexity', subtitle = '') + 
  theme(text = element_text(color = "#333333"),
        plot.margin=unit(c(t=3,r=3,b=5,l=5), "pt"),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(face="bold", size = 11),
        plot.subtitle = element_text(size = 15, face="bold", hjust=.5)))

(b <- ggplot() +
  geom_bar(data = urban_bchart, aes(y = share_landscan_population_un, 
                                    x = reorder(urban_name, k_reweight), fill = k_4), 
           stat="identity") +
  coord_flip() +
  scale_fill_manual(values = c(grey2, colorhexes)) +
  scale_y_continuous( expand = c(0, 0), labels = scales::percent ) +
  theme_bw() + 
  labs(y = 'Population', x = 'k-complexity', subtitle = '') + #'Population distribution across k-complexity levels'
  theme(text = element_text(color = "#333333"),
        plot.margin=unit(c(t=3,r=3,b=5,l=5), "pt"),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(face="bold", size = 11),
        plot.subtitle = element_text(size = 15, face="bold", hjust=.5)))
  
(urban_bars <- a+b)

# Conurbation bars --------------------------------------------------------

conurbation_bchart <- out_k_4 %>%
  filter(geography_group == "4 - Conurbation",
         urban_nonurban == "Core, peripheral, & peri-urban",
         area_landscan_population_un >= 1500000) %>%
  mutate(country_name =  (gsub('Democratic Republic of the Congo ', 'DR Congo', as.character(country_name))),
         conurbation_area_name_short = gsub('Ambatolampy Tsimahafotsy', 'Ambatolampy', as.character(conurbation_area_name_short))) %>%
  mutate(conurbation_name = paste0(conurbation_area_name_short,', ',country_name)) %>%
  mutate(k_reweight = landscan_population_un*k_complexity_weighted_landscan) %>%
  group_by(geography_area) %>%
  mutate(k_reweight = sum(k_reweight)) %>%
  ungroup() %>%
  mutate(k_reweight = k_reweight/area_landscan_population_un) 

grey2 <- c('#414141','#777777')
kdist = max(as.integer(conurbation_bchart$k_4))
colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(conurbation_bchart$k_4))-2)

(a <- ggplot() +
    geom_bar(data = conurbation_bchart, aes(y = landscan_population_un, x = reorder(conurbation_name, area_landscan_population_un), 
                                fill = k_4), 
             stat="identity") +
    coord_flip() +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous( expand = c(0, 0), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k-complexity', subtitle = '') + 
    theme(text = element_text(color = "#333333"),
          plot.margin=unit(c(t=3,r=3,b=5,l=5), "pt"),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          legend.title = element_blank(),
          legend.position = 'none',
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(face="bold", size = 11),
          plot.subtitle = element_text(size = 15, face="bold", hjust=.5)))

(b <- ggplot() +
    geom_bar(data = conurbation_bchart, aes(y = share_landscan_population_un, 
                                x = reorder(conurbation_name, k_reweight), fill = k_4), 
             stat="identity") +
    coord_flip() +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous( expand = c(0, 0), labels = scales::percent ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k-complexity', subtitle = '') + #'Population distribution across k-complexity levels'
    theme(text = element_text(color = "#333333"),
          plot.margin=unit(c(t=3,r=3,b=5,l=5), "pt"),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          # axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(face="bold", size = 11),
          plot.subtitle = element_text(size = 15, face="bold", hjust=.5)))

(conurbation_bars <- a+b)

ggsave(plot = africa_hist_k_5, filename = '/Users/nm/Desktop/africa_hist.png',dpi = 300, width = 11, height = 7)
ggsave(plot = africa_x_y_scatter, filename = '/Users/nm/Desktop/africa_scatter.png',dpi = 300, width = 7, height = 7)
ggsave(plot = africa_dots_xpop_sizek, filename = '/Users/nm/Desktop/africa_xpop_sizek.png',dpi = 300, width = 8, height = 8)
ggsave(plot = africa_dots_xk_sizepop, filename = '/Users/nm/Desktop/africa_xk_sizepop.png',dpi = 300, width = 8, height = 8)
ggsave(plot = conurbation_dots_xpop_sizek, filename = '/Users/nm/Desktop/conurb_xpop_sizek.png',dpi = 300, width = 8, height = 8)
ggsave(plot = conurbation_dots_xk_sizepop, filename = '/Users/nm/Desktop/conurb_xk_sizepop.png',dpi = 300, width = 8, height = 8)
ggsave(plot = urban_bars, filename = '/Users/nm/Desktop/urban_bars.png',dpi = 300, width = 12, height = 9)
ggsave(plot = conurbation_bars, filename = '/Users/nm/Desktop/conurbation_bars.png',dpi = 300, width = 12, height = 9)


# -------------------------------------------------------------------------




 

















#africa_x_y_scatter + 
#  guides(color = guide_legend(FALSE), fill = guide_legend(FALSE) )

# africa_hist_k_5 + africa_x_y_scatter # remove fill color legend from scatter, add title / caption in patchwork
# conurbation_x_y_scatter # add histogram for only conurbations, add title / caption in patchwork
# africa_dots_xpop_sizek + africa_dots_xk_sizepop # gather fill color and put at top, add title / caption in patchwork
# # + plot_layout(guides = "collect")
# conurbation_dots_xpop_sizek + conurbation_dots_xk_sizepop # gather fill color and put at top, add title / caption in patchwork
# urban_bars # add title / caption in patchwork
# conurbation_bars # add title / caption in patchwork

# histograms by k by urban hierarchy for each conurbation over 5 million in a facet plot 12

# ridge plots /

# maps of lagos conurbation segmenting different urban hierarchy types
# map of 1 block zoomed in to show K
# show K, show pop, show pop density



# country_name =  (gsub('Democratic Republic of the Congo', 'DR Congo', as.character(country_name)))
# country_name =  (gsub('Morocco', "Moroccan Sahara", as.character(country_name)))
# conurbation_area_name_short = gsub('Ambatolampy Tsimahafotsy', 'Ambatolampy', as.character(conurbation_area_name_short))
# urban_center_name = gsub('Ambatolampy Tsimahafotsy', 'Ambatolampy', as.character(urban_center_name))
# conurbation_area_name_short = gsub('MBour', "M'Bour" , as.character(conurbation_area_name_short))
# urban_center_name = gsub('MBour', "M'Bour" , as.character(urban_center_name))





# Collect legends
# maps
# Violin plots of distributions by country or urban area

# Bar charts side by side using k4


# Urban Centers -----------------------------------------------------------

# Bar charts side by side using k4


# -------------------------------------------------------------------------

#urban_id, conurbation_id, country_code

# Ranking -----------------------------------------------------------------

#mutate_all(~replace(., is.nan(.), 0))

# group_by_at(vars(group_var, all_of(pivot_col))) %>%
# mutate(rank_sum_landscan_population_un = row_number(desc(sum_landscan_population_un)),
#        rank_sum_worldpop_population_un = row_number(desc(sum_worldpop_population_un))) %>%
# ungroup() %>%
# mutate(largest_100_sum_landscan_population_un = case_when(rank_sum_landscan_population_un <= 100 ~ 1, TRUE ~ as.double(0)),
#        largest_100_sum_worldpop_population_un = case_when(rank_sum_worldpop_population_un <= 100 ~ 1, TRUE ~ as.double(0))) %>%
# mutate(above_500k_population = case_when(sum_landscan_population_un >= 500000 ~ 1, TRUE ~ as.double(0)))# %>%
#   mutate(subgroup_var = case_when( 
#     group_var == '5 - Country by urban / non-urban' & urban_nonurban == 'Urban' ~ '5a - Country by urban',
#     group_var == '5 - Country by urban / non-urban' & urban_nonurban == 'Non-urban' ~ '5b - Country by non-urban',
#     group_var == '6 - Urban area' & urban_nonurban == 'Urban' ~ '6a - Urban area (urban only)', 
#     group_var == '7 - Urban area (core / periphery)' & core_periphery == 'Core' ~ '7a - Urban area by core',
#     group_var == '7 - Urban area (core / periphery)' & core_periphery == 'Periphery' ~ '7b - Urban area by periphery',
#     group_var == '7 - Urban area (core / periphery)' & urban_nonurban == 'Non-urban' ~ '7c - Urban area by non-urban',
#     group_var == '8 - Urban agglomeration' & urban_nonurban == 'Urban' ~ '8 - Urban agglomeration (urban only)', 
#     TRUE ~ as.character(group_var))) %>% 
#   mutate(ranking_group = case_when(group_var == '4 - Country' ~ 'A - Country - Rank',
#                                    subgroup_var == '5a - Country by urban' ~ 'B1 - Country (urban) - Rank',
#                                    subgroup_var == '5b - Country by non-urban' ~ 'B2 - Country (non-urban) - Rank',
#                                    subgroup_var == '6a - Urban area (urban only)' & above_500k_population  == 1 ~ 'C - Above 500k urban areas - Rank',
#                                    subgroup_var == '7a - Urban area by core' & above_500k_population  == 1 ~ 'E1 - Above 500k urban cores - Rank',
#                                    subgroup_var == '7b - Urban area by periphery' & above_500k_population  == 1 ~ 'E2 - Above 500k urban peripheries - Rank',
#                                    subgroup_var == '7c - Urban area by non-urban' & above_500k_population == 1 ~ 'E3 - Above 500k urban non-urban - Rank',
#                                    subgroup_var == '8 - Urban agglomeration (urban only)' & above_500k_population == 1 ~ 'D - Above 500k urban agglomerations - Rank',
#                                    TRUE ~ as.character(''))) 

# data <- data %>%
#   group_by_at(vars(ranking_group, subgroup_var, all_of(pivot_col))) %>%
#   mutate(ranking_population_landscan = row_number(desc(landscan_population_un)),
#          ranking_share_landscan = row_number(desc(share_landscan_population_un))) %>%
#   ungroup() %>%
#   mutate(ranking_population_landscan = case_when(ranking_group == '' ~ NA_integer_, TRUE ~ as.integer(ranking_population_landscan)),
#          ranking_share_landscan = case_when(ranking_group == '' ~ NA_integer_, TRUE ~ as.integer(ranking_share_landscan))) 

# -------------------------------------------------------------------------







#scale_fill_distiller(palette = 'Spectral', name = 'Population', oob = scales::squish, breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) + 
#scale_color_distiller(palette = 'Spectral', oob = scales::squish, breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) 

# limits= c(1, max(zoom_zone$landscan_population_log )), 
# limits= c(1, max(zoom_zone$landscan_population_log )), 

bar_ranking_charts <- function(data, k_col) { 
  
  
  
  road_color = '#ffffff'
  grey2 <- c('#414141','#777777')
  kdist = max(as.integer(data$k_col))
  colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(data$k_col))-2)
  a <- ggplot() +
    geom_bar(data = data  , aes_string(y = landscan_population_un, x = reorder(geography_area, landscan_population_un), fill = k_col), 
             stat="identity") +
    coord_flip() +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous( expand = c(0, 0), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k-complexity', subtitle = '') + 
    theme(text = element_text(color = "#333333"),
          plot.margin=unit(c(t=3,r=3,b=5,l=5), "pt"),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          legend.title = element_blank(),
          legend.position = 'none',
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(face="bold", size = 11),
          plot.subtitle = element_text(size = 15, face="bold", hjust=.5))
  
  road_color = '#ffffff'
  grey2 <- c('#414141','#777777')
  kdist = max(as.integer(data $k_col))
  colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(data$k_col))-2)
  b <- ggplot() +
    geom_bar(data = data , aes_string(y = share_landscan_population_un, 
                                      x = reorder(geography_area, landscan_population_un), fill = k_col), 
             stat="identity") +
    coord_flip() +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous( expand = c(0, 0), labels = scales::percent ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k-complexity', subtitle = '') + #'Population distribution across k-complexity levels'
    theme(text = element_text(color = "#333333"),
          plot.margin=unit(c(t=3,r=3,b=5,l=5), "pt"),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          # axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(face="bold", size = 11),
          plot.subtitle = element_text(size = 15, face="bold", hjust=.5))
  
  c <- a + b
  return(c)
}
















NGA.15.2_1_0 #2565
NGA.15.3_1_1009 # 2544

df_x <- read_parquet('/Users/nm/Downloads/outputs/crosswalk/ghsl_crosswalk.parquet')

df_x <- df_x %>% filter(country_code == 'NGA') %>%
  filter(conurbation_id == '1099') %>%
  filter(urban_center_name == "Abuja")




"2565.0" "2544.0"
df_x %>% group_by(urban_id) %>% tally()

#add urban_id conurb id country code

unique(df_x$urban_center_name)

  # abuja double
