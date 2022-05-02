


library(ggplot2)
library(sf)
library(tidyverse)
library(viridis)
library(patchwork)
library(scales)
library(tidygeocoder)
library(readxl)
library(units)
library(ggsn)
library(arrow)
library(sfarrow)

# Paths --------------------------------------------------------------

country_code = c('DJI','COG','ZAF','NGA')

dir_path = '/Users/nm/Downloads/viz/'
iso_code = country_code[1]

blocks <-  st_read_parquet(paste0('/Users/nm/Downloads/production/outputs/blocks/blocks_',iso_code,'.parquet'))
metrics <- read_parquet(paste0('/Users/nm/Downloads/production/outputs/kindex/kindex_',iso_code,'.parquet'))
population <- read_parquet(paste0('/Users/nm/Downloads/production/outputs/population/population_',iso_code,'.parquet'))
#streets <-  st_read('/Users/nm/Downloads/production/outputs/streets/streets_ZAF.parquet')

st_write(blocks, '/Users/nm/Desktop/dji.geojson')
# Checks
sum(population$landscan_population)
sum(population$worldpop_population)


# GHSL City Data ----------------------------------------------------------

ghsl_url <- 'http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_STAT_UCDB2015MT_GLOBE_R2019A/V1-2/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.zip'

filedir <- paste0(tempdir(), '/ghsl/')
unlink(filedir, recursive = TRUE)
dir.create(filedir)
ghsl_dir <- paste0(ghsl_url)
download.file(url = ghsl_url, destfile = paste0(filedir, basename(ghsl_dir)))
unzip(paste0(filedir,basename(ghsl_dir)), exdir= filedir)

city_data <- read_csv(list.files(path = paste0(filedir,list.files(path = filedir)[1]), full.names = TRUE, all.files = TRUE)[6])

country_list <- c('AGO', 'BDI', 'BEN', 'BFA', 'BWA', 'CAF', 'CIV', 'CMR', 'COD', 'COG', 'COM', 'CPV', 'DJI', 'ERI', 'ESH', 'ETH', 'GAB', 'GHA', 'GIN', 'GMB', 'GNB', 'GNQ', 'KEN', 'LBR', 'LSO', 'MDG', 'MLI', 'MOZ', 'MRT', 'MUS', 'MWI', 'NAM', 'NER', 'NGA', 'RWA', 'SDN', 'SEN', 'SLE', 'SOM', 'SSD', 'STP', 'SWZ', 'SYC', 'TCD', 'TGO', 'TZA', 'UGA', 'ZAF', 'ZMB', 'ZWE')

city_data <- city_data %>%
  filter(CTR_MN_ISO %in% country_list) %>%
  group_by(CTR_MN_ISO) %>%
  mutate(rank = row_number(desc(P15))) %>%
  ungroup() %>%
  filter(rank <= 3) %>%
  filter(rank == 1 | P15 >= 500000) %>%
  mutate(UC_NM_MN = case_when(	
    CTR_MN_ISO == 'CPV' & UC_NM_MN == "Jo?o Teves" ~ 'João Teves',
    CTR_MN_ISO == 'CIV' & grepl('Bouak', UC_NM_MN) ~ 'Bouaké',
    CTR_MN_ISO == 'GHA' & UC_NM_MN == "Takoradi [Sekondi-Takoradi]" ~ 'Takoradi',
    CTR_MN_ISO == 'TGO' & grepl('Lom', UC_NM_MN) ~ "Lomé",
    CTR_MN_ISO == 'STP' & grepl('S?o Tom', UC_NM_MN) ~ "São Tomé",
    CTR_MN_ISO == 'SSD' & grepl("N/A", UC_NM_MN) ~ "Juba",
    TRUE ~ as.character(UC_NM_MN))) %>%
  mutate(CTR_MN_NM = case_when(
    CTR_MN_ISO == 'CIV' & grepl("d'Ivoire", CTR_MN_NM) ~ "Côte d'Ivoire",
    CTR_MN_ISO == 'STP' & grepl("S?o Tom", CTR_MN_NM) ~ "São Tomé and Príncipe",
    TRUE ~ as.character(CTR_MN_NM))) %>%
  select(CTR_MN_ISO, CTR_MN_NM, UC_NM_MN, P15, GCPNT_LAT, GCPNT_LON) %>%
  st_as_sf(coords = c("GCPNT_LON", "GCPNT_LAT"), 
           crs = 4326, agr = "constant") %>%
  rename(country_iso = CTR_MN_ISO,
         country_name = CTR_MN_NM, 
         city_name = UC_NM_MN,
         ghsl_pop15 = P15)

# UN Population Data ------------------------------------------------------

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

# Load and prep data ------------------------------------------------------

blocks <- blocks %>%
  left_join(., un_population %>% select(iso3_alpha_code, poptotal, popdensity), 
            by = c('country_code' = 'iso3_alpha_code'))

data <- blocks %>% 
  left_join(x =., y = metrics %>% select(block_id, block_area, building_area, building_count, building_layers, k_complexity), 
            by = c('block_id' = 'block_id')) %>%
  left_join(x =., y = population %>% select(block_id, landscan_population, worldpop_population),
            by = c('block_id' = 'block_id'))

data <- data %>% 
  group_by(country_code) %>%
  mutate(landscan_population_total = sum(landscan_population),
         worldpop_population_total = sum(worldpop_population)) %>%
  ungroup() %>%
  mutate(landscan_population = round(poptotal * landscan_population/landscan_population_total,0),
         worldpop_population = round(poptotal * worldpop_population/worldpop_population_total,0)) %>%
  select(-one_of(c('landscan_population_total','worldpop_population_total','poptotal','popdensity'))) %>%
  mutate(landscan_population_log = replace_na(na_if((log10(landscan_population)), -Inf),1),
         landscan_population_log = ifelse(landscan_population_log < 1, 1, landscan_population_log),
         worldpop_population_log = replace_na(na_if((log10(worldpop_population)), -Inf),1),
         worldpop_population_log = ifelse(worldpop_population_log < 1, 1, worldpop_population_log)) %>%
  mutate(block_hectares = na_if((block_area*0.0001), 0)) %>%
  replace_na(list(block_hectares = 0)) %>% 
  mutate(landscan_pop_density_hectare = landscan_population/block_hectares,
         worldpop_pop_density_hectare = worldpop_population/block_hectares) %>%
  mutate(landscan_pop_density_hectare_log = replace_na(na_if((log10(landscan_pop_density_hectare)), -Inf),1),
         landscan_pop_density_hectare_log = ifelse(landscan_pop_density_hectare_log < 1, 1, landscan_pop_density_hectare_log),
         worldpop_pop_density_hectare_log = replace_na(na_if((log10(worldpop_pop_density_hectare)), -Inf),1),
         worldpop_pop_density_hectare_log = ifelse(worldpop_pop_density_hectare_log  < 1, 1, worldpop_pop_density_hectare_log )) 

data_sum <- data %>% st_drop_geometry() %>%
  mutate(k_complexity = as.factor(as.integer(k_complexity)),
         block_hectares = na_if((block_area*0.0001), 0)) %>%
  replace_na(list(block_hectares = 0)) %>%
  group_by(k_complexity) %>%
  summarize_at(vars(block_hectares, landscan_population, worldpop_population), list(sum), na.rm = TRUE) %>%
  ungroup() %>%
  mutate(landscan_pop_density_hectare = landscan_population/block_hectares,
         worldpop_pop_density_hectare = worldpop_population/block_hectares) %>% 
  replace_na(list(block_pop_density_hectare = 0)) %>%
  mutate(landscan_population_sum = sum(landscan_population),
         worldpop_population_sum = sum(worldpop_population),
         landscan_population_share = landscan_population/landscan_population_sum,
         worldpop_population_share = worldpop_population/worldpop_population_sum)


# City Parameters ---------------------------------------------------------

#city_list = c('Cape Town','Durban','Johannesburg') 
# c('Libreville','Port-Gentil','Franceville')
#country_name = 'South Africa' #'Gabon'
#city_name = city_list[1]
#city_coords = tidygeocoder::geo(city = city_name, country = country_name, method = "osm")

# Generate zoom area ------------------------------------------------------

# target_area <- st_point(c(city_coords$long, city_coords$lat)) %>% 
#   st_geometry() %>% st_set_crs(4326) %>%
#   st_transform(3857) %>%
#   st_buffer(., 25000) %>% st_bbox(.) %>% st_as_sfc()
#   #st_buffer(., 5280 * 1.5) %>% st_bbox(.) %>% st_as_sfc()

target_area <- city_data %>% 
  filter(country_iso %in% iso_code) %>%
  filter(city_name == "Brazzaville") %>%
  st_transform(4326) %>%
  st_transform(3857) %>%
  st_buffer(., 25000) %>% st_bbox(.) %>% st_as_sfc()

blocks_centroids <- blocks %>% st_make_valid() %>% 
  st_transform(3857) %>%
  st_centroid()

#yer_object$geometry <- yer_object$geometry %>%
#  s2::s2_rebuild() %>%
#  sf::st_as_sfc()

zoom_zone <- blocks_centroids %>%
  st_make_valid() %>%
  st_transform(3857) %>%
  mutate(in_zone = ifelse(sf::st_intersects(., target_area, sparse = F), "Yes", "No")) %>% 
  filter(in_zone == 'Yes') %>%
  st_intersection(., target_area)

zoom_zone <- data %>% filter(block_id %in% unique(zoom_zone$block_id))

#947
#'ZAF.9.3.1.98_1_14'
# theories
# buffer out the linestrings, 
# exteriors streets arent registered as streets
# too many buildings
# stranded and defaulting to buffer
z <- data %>% filter(k_complexity >= 7) %>% st_drop_geometry() %>%
  separate(col = building_layers, into = c('first','second','third'), sep = ',', remove = FALSE) %>%
  mutate_at(c('first','second','third'), as.integer) %>%
  filter(first <= 9)
# write_csv(z, '/Users/nm/Desktop/qcblocks_ZAF.csv')
# #12739/2416471
# z <- zoom_zone %>% filter(str_sub(gadm_code, 1, 9) == 'ZAF.9.3.1') %>%
#   filter(k_complexity)
# st_write(zoom_zone %>% filter(k_complexity >= 15),'/Users/nm/Desktop/ct.geojson')


zoom_zone <- zoom_zone %>% 
  mutate(k_complexity_groups = case_when(k_complexity >= 11 & k_complexity <= 15 ~ "11-15",
                                         k_complexity >= 16 & k_complexity <= 20 ~ "16-20",
                                         k_complexity >= 21 ~ "21+",
                                         TRUE ~ as.character(k_complexity))) %>% 
  arrange(k_complexity)

k_order = unique(zoom_zone$k_complexity_groups)
zoom_zone <- zoom_zone %>% 
  mutate(k_complexity_groups = factor(k_complexity_groups,levels = k_order)) 

st_bbox(zoom_zone %>% st_transform(3395))
#zoom_zone_area = zoom_zone %>% select(geometry) %>% st_union() %>% st_transform(3857) %>% st_area()*1e-6 
#units(zoom_zone_area) <- NULL

zoom_zone_sum <- zoom_zone %>% st_drop_geometry() %>%
  mutate(#k_complexity = as.factor(as.integer(k_complexity)),
         k_complexity = k_complexity_groups,
         block_hectares = na_if((block_area*0.0001), 0)) %>%
  replace_na(list(block_hectares = 0)) %>%
  group_by(k_complexity) %>%
  summarize_at(vars(block_hectares, landscan_population, worldpop_population), list(sum), na.rm = TRUE) %>%
  ungroup() %>%
  mutate(landscan_pop_density_hectare = landscan_population/block_hectares,
         worldpop_pop_density_hectare = worldpop_population/block_hectares) %>% 
  replace_na(list(block_pop_density_hectare = 0)) %>%
  mutate(landscan_population_sum = sum(landscan_population),
         worldpop_population_sum = sum(worldpop_population),
         landscan_population_share = landscan_population/landscan_population_sum,
         worldpop_population_share = worldpop_population/worldpop_population_sum) %>%
  mutate()
         

# Visualize zoom area -----------------------------------------------------

road_color = '#ffffff'
grey2 <- c('#414141','#777777')
kdist = max(as.integer(zoom_zone$k_complexity_groups))
colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(zoom_zone$k_complexity_groups))-2)

(plot_k_discrete <- ggplot() +
    geom_sf(data = zoom_zone, aes(fill = as.factor(k_complexity_groups)), color = road_color, size = .0075) +   
    scale_fill_manual(values = c(grey2,colorhexes), name = 'k complexity') + 
    labs(subtitle = '',
         caption = paste0('Population weighted average k complexity in area: ',zoom_zone %>% st_drop_geometry() %>% summarise(wm_var = weighted.mean(as.integer(k_complexity), landscan_population)) %>% pull() %>% round(.,2))) +
    theme_void() + theme(#plot.subtitle = element_text(size = 15, face="bold", vjust = -4, hjust=.5),
                         plot.caption = element_text(size = 9, hjust = .5, vjust = 10,margin=margin(0,0,0,0)),
                         text = element_text(color = "#333333"),
                         legend.position = c(1.1,.5),
                         #legend.box.margin =unit(c(t=5,r=40,b=5,l=10), "pt"),
                         legend.key.height= unit(13, 'pt'),
                         legend.key.width= unit(13, 'pt'),
                         #axis.title.x = element_blank(),
                         legend.text = element_text(size = 8),
                         panel.border = element_blank(),
                         panel.background = element_blank(),
                         plot.margin=unit(c(t=0,r=40,b=0,l=0), "pt"),
                         legend.title = element_text( face="bold", size = 10),
                         axis.text = element_blank()) +
    ggsn::scalebar(y.min = st_bbox(zoom_zone)$ymin - .003,x.min = st_bbox(zoom_zone)$xmin, 
                   y.max = st_bbox(zoom_zone)$ymax, x.max = st_bbox(zoom_zone)$xmax, location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = 1, dist_unit = "km", transform = TRUE, model = "WGS84")
  )

#scale_cut = c(0, k = 1e3, m = 1e6, bn = 1e9, tn = 1e12)


(bar_k_distrib <- ggplot(zoom_zone_sum) +
    geom_bar(aes(y = landscan_population, x = k_complexity, fill = k_complexity), 
             position="dodge",  stat="identity") +
    geom_text(aes(x = k_complexity, y =landscan_population, 
                  label = ifelse(landscan_population_share > .01, paste0(round(landscan_population_share*100,0),"%"),'')),
              size = 3, vjust =-.5, color = '#333333', fontface='bold') +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                       expand = expansion(mult = c(0, .1)),
                       limits = c(0, max(zoom_zone_sum$landscan_population)),
                       labels = label_comma(accuracy = 1L, scale = .001, suffix = "K") ) +
    theme_bw() + 
    labs(y = 'Population', x = 'k complexity', subtitle = '') + #'Population distribution across k-complexity levels'
    theme(text = element_text(color = "#333333"),
          legend.position= "none",
          #plot.margin = unit(c(1,1,1,1),"cm"),
          axis.ticks =element_blank(),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_text(face="bold", size = 11),
          plot.subtitle = element_text(size = 15, face="bold", hjust=.5)))

(plot_populaton <- ggplot() +
    geom_sf(data = zoom_zone,
            aes(fill = landscan_population_log ), 
            color = 'white', size= .0075, alpha = .8) +
    labs(subtitle = '',
         caption = paste0('Total population in area: ',comma(sum(zoom_zone$landscan_population)),'  |  ',
         'Average block population: ',comma(round(mean(zoom_zone$landscan_population),2)))) + 
    scale_fill_viridis(name = 'Population', oob = scales::squish, limits= c(1, max(zoom_zone$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) + 
    scale_color_viridis(name = 'Population', oob = scales::squish, limits= c(1, max(zoom_zone$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) +
    theme_void() + 
    theme(plot.subtitle = element_text(size = 14, face="bold", vjust = -4, hjust=.5),axis.text = element_blank(),
          plot.caption = element_text(size = 9, hjust = .5, vjust = 7),
          legend.position = c(1.07,.8),
          plot.margin=unit(c(t=0,r=50,b=0,l=5), "pt"),
          legend.title = element_text( face="bold"),
          text = element_text(color = "#333333")) +
    ggsn::scalebar(y.min = st_bbox(zoom_zone)$ymin - .003,x.min = st_bbox(zoom_zone)$xmin, 
                   y.max = st_bbox(zoom_zone)$ymax, x.max = st_bbox(zoom_zone)$xmax, location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = 1, dist_unit = "km", transform = TRUE, model = "WGS84"))
  
sd_int = log10(mean(zoom_zone$landscan_pop_density_hectare) + sd(zoom_zone$landscan_pop_density_hectare)*1)
(plot_popdensity_log <- ggplot() +
    geom_sf(data = zoom_zone,
            aes(fill = landscan_pop_density_hectare_log), 
            color = 'white', size= .0075, alpha = .8) +
    labs(subtitle = '',
         caption = paste0('Weighted average population density: ',
                          zoom_zone %>% st_drop_geometry() %>% summarize(pop_dense = weighted.mean(landscan_pop_density_hectare, landscan_population) ) %>% pull() %>% round(.,0),
                          ' people per hectare','\n  1 hectare = 1.4 soccer fields = 2.2 city blocks in Manhattan, New York City')) +
    scale_fill_viridis(name = 'Population\nper hectare', oob = scales::squish, limits= c(1, sd_int), 
                       breaks= c(1,2,3,4,5,6,7), 
                       labels = c('0',"100","1K","10K","100K","1M","10M")) +
                       #breaks= c(1,1.477121,2,2.477121,3,3.69897,4,5,6,7), 
                       #labels = c('0','30',"100",'300',"1K","5K","10K","100K","1M","10M")) +
    scale_color_viridis(name = 'Population\nper hectare', oob = scales::squish, limits= c(1, sd_int), 
                        breaks= c(1,2,3,4,5,6,7), 
                        labels = c('0',"100","1K","10K","100K","1M","10M")) +
                        #breaks= c(1,1.477121,2,2.477121,3,3.69897,4,5,6,7), 
                        #labels = c('0','30',"100",'300',"1K","5K","10K","100K","1M","10M")) +
    theme_void() + 
    theme(plot.subtitle = element_text(size = 14, face="bold", vjust = -4, hjust=.5),
          plot.caption = element_text(size = 9, hjust = .5, vjust = 7),
          legend.position = c(1.07,.8),
          plot.margin=unit(c(t=0,r=70,b=0,l=5), "pt"),
          axis.text = element_blank(),
          legend.title = element_text( face="bold"),
          text = element_text(color = "#333333")) +
    ggsn::scalebar(y.min = st_bbox(zoom_zone)$ymin - .003,x.min = st_bbox(zoom_zone)$xmin, 
                   y.max = st_bbox(zoom_zone)$ymax, x.max = st_bbox(zoom_zone)$xmax, location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = 1, dist_unit = "km", transform = TRUE, model = "WGS84"))

st_bbox(zoom_zone)$ymin
st_bbox(zoom_zone)$xmin

(plot_k <- plot_k_discrete + bar_k_distrib +
    plot_layout(widths = c(1, 1))  +
   plot_annotation(#title = paste0(city_name,', ', country_name),
                   subtitle = paste0(city_name,', ', country_name),
                   theme = theme(#plot.title = element_text(face="bold", size = 18, vjust = -2, hjust = .5),
                                 plot.subtitle = element_text(face="bold", size = 13, vjust = -7, hjust = .5))))
ggsave(plot = plot_k, filename = paste0(dir_path,city_name,'_k.pdf'), dpi = 600, height = 6, width = 12)

(plot_pop <- plot_populaton + plot_popdensity_log +
    plot_layout(widths = c(1, 1)) +
  plot_annotation(#title = paste0(city_name,', ', country_name ),
                  subtitle = paste0(city_name,', ', country_name ),
                  theme = theme(#plot.title = element_text(face="bold", size = 18, vjust = -2, hjust = .5),
                                plot.subtitle = element_text(face="bold", size = 13, vjust = -5, hjust = .5)))
  )
ggsave(plot = plot_pop, filename = paste0(dir_path,city_name,'_pop.pdf'), dpi = 600, height = 7, width = 12)


# Visualize country -------------------------------------------------------

pdf(paste0('/Users/nm/Desktop/viz/plot_',country_name,'.pdf') )
plot(blocks %>% select(block_id), main="", lwd=.01) 
dev.off() 

road_color = '#ffffff'
grey2 <- c('#414141','#777777')
kdist = max(as.integer(data$k_complexity))
colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(data$k_complexity))-2)
country_bbox = st_bbox(data)

plot_k_discrete <- ggplot() +
    geom_sf(data = data, aes(fill = as.factor(as.integer(k_complexity))), color = alpha(c(road_color), .9), size = .005) +   # 
    scale_fill_manual(values = c(grey2,colorhexes), name = 'k complexity') + 
    labs(subtitle = '',
         caption = paste0('Population weighted average k complexity in area: ',data %>% st_drop_geometry() %>% summarise(wm_var = weighted.mean(as.integer(k_complexity), landscan_population)) %>% pull() %>% round(.,2))) +
    theme_void() + theme(#plot.subtitle = element_text(size = 15, face="bold", vjust = -4, hjust=.5),
      plot.caption = element_text(size = 9, hjust = .5, vjust = 10,margin=margin(0,0,0,0)),
      text = element_text(color = "#333333"),
      legend.position = c(1.1,.5),
      #legend.box.margin =unit(c(t=5,r=40,b=5,l=10), "pt"),
      legend.key.height= unit(13, 'pt'),
      legend.key.width= unit(13, 'pt'),
      #axis.title.x = element_blank(),
      legend.text = element_text(size = 8),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.margin=unit(c(t=0,r=40,b=0,l=0), "pt"),
      legend.title = element_text( face="bold", size = 10),
      axis.text = element_blank()) +
    ggsn::scalebar(y.min = st_bbox(data)$ymin - .003, x.min = st_bbox(data)$xmin, 
                   y.max = st_bbox(data)$ymax, x.max = st_bbox(data)$xmax, location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = round(((country_bbox$xmax - country_bbox$xmin)*.1)/1000,-1), dist_unit = "km", transform = TRUE, model = "WGS84")

bar_k_distrib <- ggplot(data_sum) +
    geom_bar(aes(y = landscan_population, x = k_complexity, fill = k_complexity), 
             position="dodge",  stat="identity") +
    geom_text(aes(x = k_complexity, y =landscan_population, 
                  label = ifelse(landscan_population_share > .01, paste0(round(landscan_population_share*100,0),"%"),'')),
              size = 3, vjust =-.5, color = '#333333', fontface='bold') +
    scale_fill_manual(values = c(grey2, colorhexes)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 6), 
                       expand = expansion(mult = c(0, .1)),
                       limits = c(0, max(data_sum$landscan_population)),
                       labels = label_comma(accuracy = 1L, scale = .001, suffix = "K") )+ 
    theme_bw() + 
    labs(y = 'Population', x = 'k complexity', subtitle = '') + #'Population distribution across k-complexity levels'
    theme(text = element_text(color = "#333333"),
          legend.position= "none",
          #plot.margin = unit(c(1,1,1,1),"cm"),
          axis.ticks =element_blank(),
          axis.text = element_text(size = 11),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_text(face="bold", size = 11),
          plot.subtitle = element_text(size = 15, face="bold", hjust=.5))

plot_populaton <- ggplot() +
    geom_sf(data = data,
            aes(fill = landscan_population_log ), 
            color = alpha(c(road_color), .9), size = .005, alpha = .8) +
    labs(subtitle = '',
         caption = paste0('Total population in area: ',comma(sum(data$landscan_population)),'  |  ',
                          'Average block population: ',comma(round(mean(data$landscan_population),2)))) + 
    scale_fill_viridis(name = 'Population', oob = scales::squish, limits= c(1, max(data$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) + 
    scale_color_viridis(name = 'Population', oob = scales::squish, limits= c(1, max(data$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) +
    theme_void() + 
    theme(plot.subtitle = element_text(size = 14, face="bold", vjust = -4, hjust=.5),axis.text = element_blank(),
          plot.caption = element_text(size = 9, hjust = .5, vjust = 7),
          legend.position = c(1.07,.8),
          plot.margin=unit(c(t=0,r=50,b=0,l=5), "pt"),
          legend.title = element_text( face="bold"),
          text = element_text(color = "#333333")) +
    ggsn::scalebar(y.min = st_bbox(data)$ymin - .003,x.min = st_bbox(data)$xmin, 
                   y.max = st_bbox(data)$ymax, x.max = st_bbox(data)$xmax, location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = round(((country_bbox$xmax - country_bbox$xmin)*.1)/1000,-1), dist_unit = "km", transform = TRUE, model = "WGS84")

sd3 = log10(mean(data$landscan_pop_density_hectare) + sd(data$landscan_pop_density_hectare)*5)
plot_popdensity_log <- ggplot() +
    geom_sf(data = data,
            aes(fill = landscan_pop_density_hectare_log), 
            color = alpha(c(road_color), .9), size = .005, alpha = .8) +
    labs(subtitle = '',
         caption = paste0('Weighted average population density: ',
                          data %>% st_drop_geometry() %>% summarize(pop_dense = weighted.mean(landscan_pop_density_hectare, landscan_population) ) %>% pull() %>% round(.,0),
                          ' people per hectare','\n  1 hectare = 1.4 soccer fields = 2.2 city blocks in Manhattan, New York City')) +
    scale_fill_viridis(name = 'Population\nper hectare', oob = scales::squish, limits= c(1, sd3), breaks= c(1,1.477121,2,2.60206,3,3.69897,4,5,6,7), labels = c('0','30',"100",'400',"1K","5K","10K","100K","1M","10M")) +
    scale_color_viridis(name = 'Population\nper hectare', oob = scales::squish, limits= c(1, sd3), breaks= c(1,1.477121,2,2.60206,3,3.69897,4,5,6,7), labels = c('0','30',"100",'400',"1K","5K","10K","100K","1M","10M")) +
    theme_void() + 
    theme(plot.subtitle = element_text(size = 14, face="bold", vjust = -4, hjust=.5),
          plot.caption = element_text(size = 9, hjust = .5, vjust = 7),
          legend.position = c(1.07,.8),
          plot.margin=unit(c(t=0,r=70,b=0,l=5), "pt"),
          axis.text = element_blank(),
          legend.title = element_text( face="bold"),
          text = element_text(color = "#333333")) +
    ggsn::scalebar(y.min = st_bbox(data)$ymin - .003,x.min = st_bbox(data)$xmin, 
                   y.max = st_bbox(data)$ymax, x.max = st_bbox(data)$xmax, location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = round(((country_bbox$xmax - country_bbox$xmin)*.1)/1000,-1), dist_unit = "km", transform = TRUE, model = "WGS84")

plot_k <- plot_k_discrete + bar_k_distrib +
    plot_layout(widths = c(1, 1))  +
    plot_annotation(#title = paste0(country_name),
      subtitle = paste0(country_name),
      theme = theme(#plot.title = element_text(face="bold", size = 18, vjust = -2, hjust = .5),
        plot.subtitle = element_text(face="bold", size = 13, vjust = -8, hjust = .5)))
ggsave(plot = plot_k, filename = paste0('/Users/nm/Desktop/GAB/plotk_',country_name,'.pdf'), dpi = 600, height = 6, width = 12)

plot_pop <- plot_populaton + plot_popdensity_log +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(#title = paste0(country_name ),
      subtitle = paste0(country_name ),
      theme = theme(#plot.title = element_text(face="bold", size = 18, vjust = -2, hjust = .5),
        plot.subtitle = element_text(face="bold", size = 13, vjust = -8, hjust = .5)))
ggsave(plot = plot_pop, filename = paste0('/Users/nm/Desktop/GAB/plotp_',country_name,'.pdf'), dpi = 600, height = 7, width = 12)

bar_k_distrib

data_sum %>% select(k_complexity, landscan_population) %>% print(n = 100)
# after 10 collapse 
# 
# 
# library(osmdata)
# bldgs <- opq(bbox = st_bbox(c(xmin =3.3693, xmax = 3.3923, 
#                               ymax = 6.4385, ymin = 6.4249), crs = st_crs(4326))) %>%
#   add_osm_feature(key = 'buildings') %>%
#   osmdata_sf() 
# bldgs_polygons <- bldgs$osm_polygons %>% dplyr::select(osm_id)









