

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
options(scipen = 999)

location_url = 'https://population.un.org/wpp/Download/Files/4_Metadata/WPP2019_F01_LOCATIONS.XLSX'
tmp_filepath <- paste0(tempdir(), '/', basename(location_url))
download.file(url = paste0(location_url), destfile = tmp_filepath)
location_codes <- read_xlsx(path = tmp_filepath, sheet = 'Location', range = "A17:H306") 

location_codes <- location_codes %>%
  select_all(~gsub("\\s+|\\.|\\/|,|\\*|-", "_", .)) %>%
  rename_all(list(tolower)) %>%
  filter(name == 'Country/Area') %>%
  rename(country_name = region__subregion__country_or_area_) %>%
  select(location_code, country_name, iso3_alpha_code)

un_population <- read_csv('https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_TotalPopulationBySex.csv')

un_population <- un_population %>%
  select_all(~gsub("\\s+|\\.|\\/|,|\\*", "_", .)) %>%
  rename_all(list(tolower)) %>%
  filter(time == 2020,  
         variant == 'Medium') %>%
  mutate(poptotal = poptotal * 1000,
         popdensity = ((popdensity*.3861)/100)*1000) %>%
  inner_join(., location_codes, by = c('locid' = 'location_code')) %>%
  select(locid, country_name, iso3_alpha_code, location, time, variant, poptotal, popdensity) 


# Paths --------------------------------------------------------------

#blocks <-  st_read_parquet(paste0('/Users/nm/Downloads/production/outputs/blocks/blocks_',iso_code,'.parquet'))
# streets <-  st_read('/Users/nm/Downloads/production/outputs/streets/streets_',iso_code,'.parquet')
iso_code_list <- c( 'AGO', 'BDI', 'BEN', 'BFA', 'BWA', 'CAF', 'CIV', 'CMR', 'COD', 'COG', 'COM', 'CPV', 'DJI', 'ERI', 'ESH', 'ETH', 'GAB', 'GHA', 'GIN', 'GMB', 'GNB', 'GNQ', 'KEN', 'LBR', 'LSO', 'MDG', 'MLI', 'MOZ', 'MRT', 'MUS', 'MWI', 'NAM', 'NER', 'NGA', 'RWA', 'SDN', 'SEN', 'SLE', 'SOM', 'SSD', 'STP', 'SWZ', 'SYC', 'TCD', 'TGO', 'TZA', 'UGA', 'ZAF', 'ZMB', 'ZWE')        

df_combined <- purrr::map_dfr(iso_code_list,  function(i) {
  print(i)
  df_complexity <- read_parquet(paste0('/Users/nm/Downloads/outputs/complexity/complexity_',i,'.parquet')) %>%
    select(block_id, gadm_code, country_code, block_area, on_network_street_length, off_network_street_length, nearest_external_street, building_area, building_count, building_layers, k_complexity)   
  df_population <- read_parquet(paste0('/Users/nm/Downloads/outputs/population/population_',i,'.parquet')) %>% 
    select(block_id, landscan_population, worldpop_population)
  df_ghsl <- read_parquet(paste0('/Users/nm/Downloads/outputs/ghsl/ghsl_main/',i,'_ghsl.parquet')) %>%
    select(block_id, urban_center, urban_agglomeration) 
  #df_blocks = st_read_parquet(dsn = '/Users/nm/Downloads/outputs/blocks/blocks_',i,'.parquet', col_select = c('block_id', 'block_area'))
  
  df_full <- df_complexity %>% 
    left_join(x =., y = df_population, by = c('block_id' = 'block_id')) %>%
    left_join(x =., y = df_ghsl, by = c('block_id' = 'block_id'))  %>%
    left_join(., un_population %>% select(iso3_alpha_code, country_name, poptotal, popdensity), 
              by = c('country_code' = 'iso3_alpha_code')) 
  
  df_full <- df_full %>%
    group_by(country_code) %>%
    mutate(landscan_population_total = sum(landscan_population),
           worldpop_population_total = sum(worldpop_population)) %>%
    ungroup() %>%
    mutate(landscan_population_un = round(poptotal * landscan_population/landscan_population_total,0),
           worldpop_population_un = round(poptotal * worldpop_population/worldpop_population_total,0)) %>%
    select(-one_of(c('landscan_population_total','worldpop_population_total','poptotal','popdensity'))) %>%
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
    group_by(urban_center, urban_agglomeration, country_code, country_name, k_complexity, k_4, k_3, k_2, k_1, k_0) %>%
    summarize_at(vars(block_count, block_area, on_network_street_length, off_network_street_length, nearest_external_street, building_area, building_count, landscan_population, worldpop_population, landscan_population_un, worldpop_population_un), list(sum)) %>%
    ungroup() 
})

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
         k_1 = factor(k_1, levels = k_1_order))  %>%
  mutate(urban_center = case_when(urban_center == 'N/A' ~ paste0('missing'), TRUE ~ as.character(urban_center)),
         urban_agglomeration = case_when(urban_agglomeration == 'N/A' ~ paste0('missing'), TRUE ~ as.character(urban_agglomeration)),
         urban_nonurban = case_when(urban_center == 'missing' ~ 'Non-urban', TRUE ~ as.character('Urban')),
         urban_parts = case_when(urban_agglomeration == urban_center & urban_center != 'missing' ~ urban_agglomeration,
                                 urban_agglomeration != urban_center & urban_center != 'missing' ~ paste0(urban_agglomeration, ' (peripheral)'), 
                                 urban_agglomeration == 'missing' ~ paste0('Rest of ',country_name), TRUE ~ as.character(urban_agglomeration)),
         urban_center = case_when(urban_center == 'missing' ~ paste0('Rest of ',country_name), TRUE ~ as.character(urban_center)),
         urban_agglomeration = case_when(urban_agglomeration == 'missing' ~ paste0('Rest of ',country_name), TRUE ~ as.character(urban_agglomeration))) 

countrylabels_urban_center <- df_combined_prep %>% arrange(country_name) %>% group_by(urban_center) %>% 
  summarize(urban_center_country_name = paste(unique(country_name), collapse = '-')) %>% ungroup()
countrylabels_urban_parts <- df_combined_prep %>% arrange(country_name) %>% group_by(urban_parts) %>% 
  summarize(urban_parts_country_name = paste(unique(country_name), collapse = '-')) %>% ungroup()
countrylabels_urban_agglomeration <- df_combined_prep %>% arrange(country_name) %>% group_by(urban_agglomeration) %>% 
  summarize(urban_agglomeration_country_name = paste(unique(country_name), collapse = '-')) %>% ungroup()

df_combined_prep2 <- df_combined_prep %>%
  left_join(., countrylabels_urban_center, by = c('urban_center'='urban_center')) %>%
  left_join(., countrylabels_urban_parts, by = c('urban_parts'='urban_parts')) %>%
  left_join(., countrylabels_urban_agglomeration, by = c('urban_agglomeration'='urban_agglomeration')) %>%
  mutate(urban_center_country = case_when(urban_nonurban == 'Urban' ~ paste0(urban_center,', ', urban_center_country_name), urban_nonurban == 'Non-urban' ~ urban_center),
         urban_parts_country = case_when(urban_nonurban == 'Urban' ~ paste0(urban_parts,', ', urban_parts_country_name), urban_nonurban == 'Non-urban' ~ urban_parts),
         urban_core_periphery = case_when(grepl('peripheral', urban_parts_country) ~ 'Periphery', grepl('Rest of ', urban_parts_country) ~ 'Non-urban', TRUE ~ as.character('Core')),
         urban_agglomeration_country = case_when(urban_nonurban == 'Urban' ~ paste0(urban_agglomeration,', ', urban_agglomeration_country_name), urban_nonurban == 'Non-urban' ~ urban_agglomeration),
         urban_nonurban_country = paste0(country_name,' (',urban_nonurban,')'),
         continent = 'sub-Saharan Africa')


rm(data,pivot_col)

process_blocks <- function(data, pivot_col,  crosstab_list) { 
  
  # data = df_combined_prep2
  # pivot_col = 'k_1'
  
  geo_labels <- data %>%
    select(continent, country_name, urban_nonurban, urban_core_periphery, urban_nonurban_country, urban_center_country, urban_parts_country, urban_agglomeration_country,
           urban_center_country_name, urban_parts_country_name, urban_agglomeration_country_name) %>%
    distinct()
  
  geo_xwalk <- bind_rows(geo_labels %>% select(urban_nonurban) %>% distinct() %>% mutate(group_var = '2 - Urban / non-urban', group_val = urban_nonurban) ,
                         geo_labels %>% select(urban_core_periphery) %>% distinct() %>% mutate(group_var = '3 - Core / periphery', group_val = urban_core_periphery),
                         geo_labels %>% select(urban_nonurban_country, country_name, urban_nonurban) %>% distinct() %>% rename(group_val = urban_nonurban_country, country_label = country_name) %>% mutate(group_var = '5 - Country by urban / non-urban') ,
                         ##geo_labels %>% select(urban_nonurban_country, country_name, urban_nonurban) %>% distinct() %>% rename(group_val = urban_nonurban_country, country_label = country_name) %>% filter(urban_nonurban == 'Urban') %>% mutate(group_var = '5a - Country by urban') ,
                         #geo_labels %>% select(urban_nonurban_country, country_name, urban_nonurban) %>% distinct() %>% rename(group_val = urban_nonurban_country, country_label = country_name) %>% filter(urban_nonurban == 'Non-urban') %>% mutate(group_var = '5b - Country by non-urban') ,
                         geo_labels %>% select(urban_center_country, urban_agglomeration_country, urban_center_country_name, urban_nonurban, urban_core_periphery) %>% distinct() %>% rename(group_val = urban_center_country, country_label = urban_center_country_name) %>% mutate(group_var = '6 - Urban area') ,
                         geo_labels %>% select(urban_parts_country, urban_agglomeration_country, urban_parts_country_name, urban_nonurban, urban_core_periphery) %>% distinct() %>% rename(group_val = urban_parts_country, country_label = urban_parts_country_name) %>% mutate(group_var = '7 - Urban area (core / periphery)') ,
                         ##geo_labels %>% select(urban_parts_country, urban_agglomeration_country, urban_parts_country_name, urban_nonurban, urban_core_periphery) %>% distinct() %>% rename(group_val = urban_parts_country, country_label = urban_parts_country_name) %>% filter(urban_core_periphery == 'Core') %>% mutate(group_var = '7a - Urban area by core',) ,
                         ##geo_labels %>% select(urban_parts_country, urban_agglomeration_country, urban_parts_country_name, urban_nonurban, urban_core_periphery) %>% distinct() %>% rename(group_val = urban_parts_country, country_label = urban_parts_country_name) %>% filter(urban_core_periphery == 'Periphery') %>% mutate(group_var = '7b - Urban area by periphery') ,
                         ##geo_labels %>% select(urban_parts_country, urban_agglomeration_country, urban_parts_country_name, urban_nonurban, urban_core_periphery) %>% distinct() %>% rename(group_val = urban_parts_country, country_label = urban_parts_country_name) %>% filter(urban_nonurban == 'Non-urban') %>% mutate(group_var =  '7c - Urban area by non-urban') ,
                         geo_labels %>% select(urban_agglomeration_country, urban_agglomeration_country_name, urban_nonurban) %>% distinct() %>% rename(group_val = urban_agglomeration_country, country_label = urban_agglomeration_country_name) %>% mutate(group_var = '8 - Urban agglomeration')) %>%
    rename(urban_agglomeration = urban_agglomeration_country, 
           core_periphery = urban_core_periphery) %>%
    relocate(group_var, group_val, urban_agglomeration, country_label, urban_nonurban, core_periphery)
  
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
    mutate(group_var = case_when(group_var == 'continent' ~ '1 - Continent',
                                 group_var == 'urban_nonurban' ~ '2 - Urban / non-urban',
                                 group_var == 'urban_core_periphery' ~ '3 - Core / periphery',
                                 group_var == 'country_name' ~ '4 - Country',
                                 group_var == 'urban_nonurban_country' ~ '5 - Country by urban / non-urban',
                                 group_var == 'urban_center_country' ~ '6 - Urban area',
                                 group_var == 'urban_parts_country'  ~ '7 - Urban area (core / periphery)',
                                 group_var == 'urban_agglomeration_country' ~ '8 - Urban agglomeration')) %>%
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
    select(-one_of(c('sum_block_count', 'sum_block_area', 'sum_on_network_street_length', 'sum_off_network_street_length', 'sum_nearest_external_street', 'sum_building_area', 'sum_building_count')))  %>%
    left_join(., geo_xwalk, by = c('group_var'='group_var', 'group_val'='group_val')) %>%
    group_by_at(vars(group_var, all_of(pivot_col))) %>%
    # mutate(rank_sum_landscan_population_un = row_number(desc(sum_landscan_population_un)),
    #        rank_sum_worldpop_population_un = row_number(desc(sum_worldpop_population_un))) %>%
    # ungroup() %>%
    # mutate(largest_100_sum_landscan_population_un = case_when(rank_sum_landscan_population_un <= 100 ~ 1, TRUE ~ as.double(0)),
    #        largest_100_sum_worldpop_population_un = case_when(rank_sum_worldpop_population_un <= 100 ~ 1, TRUE ~ as.double(0))) %>%
    mutate(above_500k_population = case_when(sum_landscan_population_un >= 500000 ~ 1, TRUE ~ as.double(0))) %>% 
    mutate(subgroup_var = case_when( 
      group_var == '5 - Country by urban / non-urban' & urban_nonurban == 'Urban' ~ '5a - Country by urban',
      group_var == '5 - Country by urban / non-urban' & urban_nonurban == 'Non-urban' ~ '5b - Country by non-urban',
      group_var == '6 - Urban area' & urban_nonurban == 'Urban' ~ '6a - Urban area (urban only)', 
      group_var == '7 - Urban area (core / periphery)' & core_periphery == 'Core' ~ '7a - Urban area by core',
      group_var == '7 - Urban area (core / periphery)' & core_periphery == 'Periphery' ~ '7b - Urban area by periphery',
      group_var == '7 - Urban area (core / periphery)' & urban_nonurban == 'Non-urban' ~ '7c - Urban area by non-urban',
      group_var == '8 - Urban agglomeration' & urban_nonurban == 'Urban' ~ '8 - Urban agglomeration (urban only)', 
      TRUE ~ as.character(group_var))) %>% 
    mutate(ranking_group = case_when(group_var == '4 - Country' ~ 'A - Country - Rank',
                                     subgroup_var == '5a - Country by urban' ~ 'B1 - Country (urban) - Rank',
                                     subgroup_var == '5b - Country by non-urban' ~ 'B2 - Country (non-urban) - Rank',
                                     subgroup_var == '6a - Urban area (urban only)' & above_500k_population  == 1 ~ 'C - Above 500k urban areas - Rank',
                                     subgroup_var == '7a - Urban area by core' & above_500k_population  == 1 ~ 'E1 - Above 500k urban cores - Rank',
                                     subgroup_var == '7b - Urban area by periphery' & above_500k_population  == 1 ~ 'E2 - Above 500k urban peripheries - Rank',
                                     subgroup_var == '7c - Urban area by non-urban' & above_500k_population == 1 ~ 'E3 - Above 500k urban non-urban - Rank',
                                     subgroup_var == '8 - Urban agglomeration (urban only)' & above_500k_population == 1 ~ 'D - Above 500k urban agglomerations - Rank',
                                     TRUE ~ as.character(''))) 
  
  data <- data %>%
    group_by_at(vars(ranking_group, subgroup_var, all_of(pivot_col))) %>%
    mutate(ranking_population_landscan = row_number(desc(landscan_population_un)),
           ranking_share_landscan = row_number(desc(share_landscan_population_un))) %>%
    ungroup() %>%
    mutate(ranking_population_landscan = case_when(ranking_group == '' ~ NA_integer_, TRUE ~ as.integer(ranking_population_landscan)),
           ranking_share_landscan = case_when(ranking_group == '' ~ NA_integer_, TRUE ~ as.integer(ranking_share_landscan))) %>%
    relocate(any_of(c('group_var', 'subgroup_var', 'group_val', pivot_col, "urban_agglomeration", "country_label", "urban_nonurban", "core_periphery",
                      'ranking_group', 'ranking_population_landscan', 'ranking_share_landscan', 
                      'k_complexity_weighted_landscan', 'share_landscan_population_un', 'landscan_population_un', 'sum_landscan_population_un', 'landscan_pop_density_hectare', 'landscan_population_per_building_area', 'landscan_population_per_building', 
                      'k_complexity_weighted_worldpop', 'share_worldpop_population_un', 'worldpop_population_un', 'sum_worldpop_population_un', 'worldpop_pop_density_hectare', 'worldpop_population_per_building_area', 'worldpop_population_per_building',
                      'block_count', 'share_block_count', 'block_area', 'share_block_area', 'building_count', 'share_building_count', 'building_area', 'share_building_area', 'buildings_per_hectare', 'buildings_per_block', 'average_building_area', 
                      'average_on_network_street_length_per_block', 'on_network_street_length', 'share_on_network_street_length', 'off_network_street_length', 'share_off_network_street_length', 'nearest_external_street', 'share_nearest_external_street', 'nearest_external_street_block_count'))) %>%
    rename(geography_group = group_var, 
           geography_subgroup = subgroup_var,
           geography_area = group_val) %>%
    arrange(geography_group, !!!all_of(pivot_cols), ranking_group, ranking_population_landscan) # %>% 
    #mutate_all(~replace(., is.nan(.), 0))

  return(data)
}

crosstab_list_input = c('continent','country_name','urban_nonurban','urban_core_periphery','urban_nonurban_country','urban_center_country','urban_parts_country','urban_agglomeration_country')

out_k_4 <- process_blocks(data = df_combined_prep2, pivot_col = 'k_4', crosstab_list = crosstab_list_input) 
out_k_3 <- process_blocks(data = df_combined_prep2, pivot_col = 'k_3', crosstab_list = crosstab_list_input) 
out_k_2 <- process_blocks(data = df_combined_prep2, pivot_col = 'k_2', crosstab_list = crosstab_list_input) 
out_k_1 <- process_blocks(data = df_combined_prep2, pivot_col = 'k_1', crosstab_list = crosstab_list_input) 
out_k_0 <- process_blocks(data = df_combined_prep2, pivot_col = 'k_0', crosstab_list = crosstab_list_input) 
#out_k_complexity <- process_blocks(data = df_combined_prep2, pivot_col = 'k_complexity', crosstab_list = crosstab_list_input) 

file_labels = data.frame(
  column_names = c('geography_group', 'geography_subgroup', 'geography_area', 'k_', 'urban_agglomeration', 'country_label', 'urban_nonurban', 'core_periphery', 'ranking_group', 'ranking_population_landscan', 'ranking_share_landscan', 'k_complexity_weighted_landscan', 'share_landscan_population_un', 'landscan_population_un', 'sum_landscan_population_un', 'landscan_pop_density_hectare', 'landscan_population_un_per_building', 'k_complexity_weighted_worldpop', 'share_worldpop_population_un', 'worldpop_population_un', 'sum_worldpop_population_un', 'worldpop_pop_density_hectare', 'worldpop_population_un_per_building', 'block_count', 'share_block_count', 'block_area', 'share_block_area', 'building_count', 'share_building_count', 'building_area', 'share_building_area', 'buildings_per_hectare', 'buildings_per_block', 'average_building_area', 'average_on_network_street_length_per_block', 'on_network_street_length', 'share_on_network_street_length', 'off_network_street_length', 'share_off_network_street_length', 'nearest_external_street', 'share_nearest_external_street', 'nearest_external_street_block_count', 'block_area_hectares', 'above_500k_population'),
  column_labels = c('Geography Group Level', 'Geography Subgroup Level', 'Geographic Area', 'K complexity', 'Urban agglomeration', 'Country', 'Urban / non-urban', 'Core / periphery', 'Ranking Group', 'Ranking of Population (LandScan)', 'Ranking Population Share (LandScan)', 'Average K complexity (LandScan weighted)', 'Population Share (LandScan)', 'Population (LandScan)', 'Total Population in area (LandScan)', 'Population per hectare (LandScan)', 'Population per building (LandScan)', 'Average K complexity (WorldPop weighted)', 'Population Share (WorldPop)', 'Population (WorldPop)', 'Total Population in area (WorldPop)', 'Population per hectare (WorldPop)', 'Population per building (WorldPop)', 'Block count', 'Share of block count', 'Block area (meters squared)', 'Share of block area', 'Building count', 'Share of building count', 'Building area (meters squared)', 'Share of building area', 'Buildings per hectare', 'Buildings per block', 'Average building area (meters squared)', 'Average block on network street length (meters)', 'Total on network street length (meters)', 'Share of on network street length', 'Total off network street length (meters)', 'Share of off network street length', 'Nearest external street (meters)', 'Share of nearest external streets', 'Number of blocks with no on network streets', 'Block area (hectares)', 'Population above 500K (LandScan)')
)


write_xlsx(list("k4" = out_k_4, "k3" = out_k_3, "k2" = out_k_2, "k1" = out_k_1, "k0" = out_k_0, 'labels' = file_labels), #"k_complexity" = out_k_complexity, 
           col_names = TRUE, format_headers = TRUE, path = '/Users/nm/Desktop/mnp-analysis.xlsx')



# Detailed K distributions 
# continent, core / periph, urban / non urban

# Stacked bar distributions
# 50 countries
# 50 largest cities / aggloms




bar_ranking_charts <- function(data, geo_group = "4 - Country", population_filter = 10000000) { 
  dist_viz = data 
  road_color = '#ffffff'
  grey2 <- c('#414141','#777777')
  kdist = max(as.integer(dist_viz$k_3))
  colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(dist_viz$k_3))-2)
  a <- ggplot() +
    geom_bar(data = dist_viz , aes(y = landscan_population_un, x = reorder(geography_area, landscan_population_un), fill = k_3), 
             stat="identity") +
    coord_flip() +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous( expand = c(0, 0), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k complexity', subtitle = '') + 
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
  
  dist_viz = data 
  road_color = '#ffffff'
  grey2 <- c('#414141','#777777')
  kdist = max(as.integer(dist_viz$k_3))
  colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(dist_viz$k_3))-2)
  b <- ggplot() +
    geom_bar(data = dist_viz , aes(y = share_landscan_population_un, x = reorder(geography_area, landscan_population_un), fill = k_3), 
             stat="identity") +
    coord_flip() +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous( expand = c(0, 0), labels = scales::percent ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k complexity', subtitle = '') + #'Population distribution across k-complexity levels'
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

unique(out_k_3$ranking_group)
unique(out_k_3$geography_group)


dir_path = '/Users/nm/Desktop/'
viz <- out_k_3 %>% filter(geography_group == "4 - Country") %>% filter(sum_landscan_population_un >= 10000000)
country_bars <- bar_ranking_charts(data = viz)
country_bars
ggsave(plot = country_bars, filename = paste0(dir_path,'country.png'), dpi = 300, height = 8, width = 14)

viz <- out_k_3 %>% filter(geography_group == "5 - Country by urban / non-urban", ranking_group == "B1 - Country (urban) - Rank"  ) 
countryurb_bars <- bar_ranking_charts(data = viz)
countryurb_bars
ggsave(plot = countryurb_bars, filename = paste0(dir_path,'countryurb.png'), dpi = 300, height = 8, width = 14)

viz <- out_k_3 %>% filter(geography_group == "6 - Urban area", ranking_group == "C - Above 500k urban areas - Rank" ) %>% filter(sum_landscan_population_un >= 1000000)
urban_bars <- bar_ranking_charts(data = viz)
urban_bars
ggsave(plot = urban_bars, filename = paste0(dir_path,'urban.png'), dpi = 300, height = 9, width = 16)


viz <- out_k_3 %>% filter(geography_group == "8 - Urban agglomeration", ranking_group == "D - Above 500k urban agglomerations - Rank" ) %>% filter(sum_landscan_population_un >= 1000000)
agglom_bars <- bar_ranking_charts(data = viz)
agglom_bars
ggsave(plot = agglom_bars, filename = paste0(dir_path,'agglom.png'), dpi = 300, height = 9, width = 16)






viz <- out_k_3 %>% filter(geography_group == "7 - Urban area (core / periphery)" ,
                          ranking_group =="E1 - Above 500k urban cores - Rank"    ) %>%
  filter(sum_landscan_population_un >= 1000000)

urbcore_bars <- bar_ranking_charts(data = viz)
urbcore_bars

viz <- out_k_3 %>% filter(geography_group == "7 - Urban area (core / periphery)" ,
                          ranking_group == "E2 - Above 500k urban peripheries - Rank" )  %>%
  filter(sum_landscan_population_un >= 200000)

urbperiph_bars <- bar_ranking_charts(data = viz)
urbperiph_bars










 






# %>%
#   left_join(., out_k_1 %>% filter(geography_group == "4 - Country", k_1 == 'Total') %>%
#               select(geography_area, k_complexity_weighted_landscan)%>% 
#               rename(order_val = k_complexity_weighted_landscan), by = c('geography_area'='geography_area'))

# geom_text(aes(x = k_complexity_groups, y =landscan_population_un, 
#               label = ifelse(share > .01, paste0(round(share*100,0),"%"),'')),
#           size = 3, vjust =-.5, color = '#333333', fontface='bold') +

  
ggplot(df_urban_thresholds  %>% filter(landscan_population_un_total >= 2000000)) +
  geom_bar(aes(y = landscan_population_un, x = reorder(urban_label, landscan_population_un_total), 
               fill = k_thresholds), 
           stat="identity") +
  coord_flip() +
  #scale_fill_manual() +
  scale_y_continuous( expand = c(0, 0), labels = label_comma(accuracy = 1L, scale =  0.000001, suffix = "M") ) +
  theme_bw() + 
  labs(y = 'Population', x = 'k complexity', subtitle = '') + #'Population distribution across k-complexity levels'
  theme(text = element_text(color = "#333333"),
        plot.margin=unit(c(t=3,r=3,b=5,l=5), "pt"),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        #panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(face="bold", size = 11),
        plot.subtitle = element_text(size = 15, face="bold", hjust=.5))





 ggplot(df_continent) +
    geom_bar(aes(y = landscan_population_un, x = k_complexity_groups, fill = k_complexity_groups), 
             position="dodge",  stat="identity") +
    geom_text(aes(x = k_complexity_groups, y =landscan_population_un, 
                  label = ifelse(share > .01, paste0(round(share*100,0),"%"),'')),
              size = 3, vjust =-.5, color = '#333333', fontface='bold') +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                       expand = expansion(mult = c(0, .1)),
                       limits = c(0, max(df_continent$landscan_population_un)),
                       labels = label_comma(accuracy = 1L, scale = 1e-06, suffix = "M") ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k complexity', subtitle = '') + #'Population distribution across k-complexity levels'
    theme(text = element_text(color = "#333333"),
          legend.position= "none",
          plot.margin=unit(c(t=0,r=5,b=0,l=5), "pt"),
          axis.ticks =element_blank(),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_text(face="bold", size = 11),
          plot.subtitle = element_text(size = 15, face="bold", hjust=.5))


 
# names(df_continent)
# write_csv(df_continent, '/Users/nm/Desktop/cont.csv')




df_urban_wide <- df_urban %>%
  select(urban_label, landscan_population_un_total, k_binary, landscan_population_un_k_binary) %>%
  distinct() %>%
  pivot_wider(id_cols = c(urban_label, landscan_population_un_total),
              names_from = c(k_binary), 
              values_from = c(landscan_population_un_k_binary),
              values_fill = 0) %>%
  mutate(share_1_3 = `1-3`/landscan_population_un_total,
         share_4plus = `4+`/landscan_population_un_total) %>%
  filter(landscan_population_un_total >= 1000000)

df_urban_wide %>% arrange(-share_4plus) %>% select(urban_label, `4+`, share_4plus)
df_urban_wide %>% arrange(-`4+`) %>% select(urban_label,`4+`, share_4plus)

