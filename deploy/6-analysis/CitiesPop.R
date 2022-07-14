
library(sf)
library(dplyr)
library(ggplot2)
library(purrr)
library(tigris)
library(tidycensus)
library(stringr)
library(readxl)
library(viridis)
library(scales)
library(tidyverse)

readRenviron("~/.Renviron")

# States ------------------------------------------------------------------

state_xwalk <- as.data.frame(fips_codes) %>%
  rename(state_fips = state_code,
         state_codes = state,
         county_name = county) %>%
  mutate(county_fips = paste0(state_fips,county_code))
state_fips <- unique(state_xwalk$state_fips)[1:51]
state_codes <- unique(state_xwalk$state_codes)[1:51]

us_states <- get_acs(year = 2020, geography = "state", variables = "B01003_001", geometry = TRUE, keep_geo_vars = TRUE, shift_geo = TRUE)

us_states <- us_states %>%
  rename_all(tolower) %>%
  rename(population = estimate,
         state_fips = geoid)%>%
  mutate_at(vars(state_fips),list(as.character)) %>%
  mutate(state_fips = str_pad(state_fips, width=2, side="left", pad="0")) %>%
  select(state_fips, population) %>%
  rename(state_population = population)

us_states <- left_join(us_states, 
                       state_xwalk %>% 
                         select(state_codes, state_fips, state_name) %>%
                         distinct() %>%
                         filter(!(state_fips %in% c('60','66','69','72','74','78'))),
                       by = c('state_fips'='state_fips')) %>%
  st_transform(crs = st_crs(4326)) %>% 
  st_as_sf()

# Counties ----------------------------------------------------------------

xwalk_url <- 'https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2020/delineation-files/list1_2020.xls'
tmp_filepath <- paste0(tempdir(), '/', basename(xwalk_url))
download.file(url = paste0(xwalk_url), destfile = tmp_filepath)
unzip(tmp_filepath, exdir=tempdir())
cbsa_xwalk <- read_excel(tmp_filepath, sheet = 1, range = cell_rows(3:1919))
cbsa_xwalk <- cbsa_xwalk %>% 
  select_all(~gsub("\\s+|\\.|\\/", "_", .)) %>%
  rename_all(list(tolower)) %>%
  mutate(county_fips = paste0(fips_state_code,fips_county_code)) %>%
  rename(cbsa_fips = cbsa_code,
         area_type = metropolitan_micropolitan_statistical_area) %>%
  select(county_fips,cbsa_fips,cbsa_title,area_type,central_outlying_county) 

us_county <- get_acs(year = 2018, geography = "county", variables = "B01003_001", geometry = TRUE, keep_geo_vars = TRUE, shift_geo = TRUE)
us_county <- us_county %>%
  rename_all(list(tolower)) %>%
  rename(county_fips = geoid,
         county_population = estimate) %>%
  select(county_fips,county_population)

# ggplot(us_county, aes(fill = log(county_population), color = log(county_population))) +
#   geom_sf() + scale_fill_viridis() + scale_color_viridis() 

us_county <- us_county %>% 
  left_join(., cbsa_xwalk, by = c('county_fips'='county_fips') ) %>%
  left_join(., state_xwalk, by = c('county_fips'='county_fips') ) %>%
  mutate(area_type = case_when(is.na(area_type) ~ 'Rural',
                               area_type == 'Metropolitan Statistical Area' ~ 'Metro',
                               area_type == 'Micropolitan Statistical Area' ~ 'Micro'),
         central_outlying_county = ifelse(is.na(central_outlying_county), 'Rural', central_outlying_county)) %>%
  select(county_fips,county_code,county_name,county_population,
         cbsa_fips,cbsa_title,area_type,central_outlying_county,
         state_codes,state_fips,state_name) %>%
  st_transform(crs = st_crs(4326)) %>% 
  st_as_sf()




# Places ------------------------------------------------------------------

  filedir <- paste0(tempdir(), '/places/')
  unlink(filedir, recursive = TRUE)
  dir.create(filedir)
  for (s in state_fips) {
    state_shp <- paste0('https://www2.census.gov/geo/tiger/TIGER2019/PLACE/tl_2019_',s,'_place.zip')
    download.file(url = state_shp, destfile = paste0(filedir, basename(state_shp)))
    unzip(paste0(filedir,basename(state_shp)), exdir= filedir)
  }
  list.files(path = filedir)
  us_places <- st_read(fs::dir_ls(filedir, regexp = "\\.shp$")[1])
  for (f in fs::dir_ls(filedir, regexp = "\\.shp$")[-1] ) {
    state_sf <- st_read(f)
    us_places <- rbind(us_places, state_sf)
  }
 
us_places <- us_places %>% 
  rename_all(tolower) %>%
  mutate_at(vars(geoid, statefp, placefp, placens),list(as.character)) %>%
  mutate(geoid = str_pad(geoid, width=7, side="left", pad="0"),
         statefp = str_pad(statefp, width=2, side="left", pad="0"),
         placefp = str_pad(placefp, width=5, side="left", pad="0"),
         placens = str_pad(placens, width=8, side="left", pad="0")) %>%
  filter(statefp %in% state_fips) %>%
  st_transform(crs = st_crs(4326)) %>% 
  st_as_sf()



places_pop <- readr::read_csv('https://www2.census.gov/programs-surveys/popest/datasets/2010-2020/cities/totals/sub-est2020_all.csv')

places_pop <- places_pop %>%
  rename_all(tolower) %>% 
  filter(sumlev %in% c('061','162','171')) %>%
  select(state, place, name, stname, popestimate2018) %>%
  mutate(state = str_pad(state, width=2, side="left", pad="0"),
         place = str_pad(place, width=5, side="left", pad="0"),
         placeid = paste0(state,place)) %>%
  rename(cityname = name)

us_places <- inner_join(us_places, 
                        places_pop,
                        by = c('geoid'='placeid')) %>%
  st_transform(crs = st_crs(4326)) %>% 
  st_as_sf()

# Write GeoJSONs ----------------------------------------------------------



