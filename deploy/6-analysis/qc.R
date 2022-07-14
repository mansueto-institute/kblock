


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


country_list = c('BDI','BWA','CAF','COG')
country_code = country_list[1]

blocks <-  st_read(paste0('/Users/nm/Desktop/qc/blocks_',country_code,'.geojson'))
metrics <- read_csv(paste0('/Users/nm/Desktop/qc/kindex_',country_code,'.csv'))

nrow(blocks)
nrow(metrics)

data <- blocks %>% 
  left_join(x =., y = metrics %>% select(block_id, block_area, building_area, building_count, building_layers, k_complexity), 
            by = c('block_id' = 'block_id')) 

road_color = '#ffffff'
grey2 <- c('#414141','#777777')
kdist = max(as.integer(data$k_complexity))
colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(data$k_complexity))-2)
country_bbox = st_bbox(data %>% st_transform(3857))

plot_k_discrete <- ggplot() +
  geom_sf(data = data, aes(fill = as.factor(as.integer(k_complexity))), color = alpha(c(road_color), .9), size = .005) +   # 
  scale_fill_manual(values = c(grey2,colorhexes), name = 'k complexity') + 
  labs(subtitle = '') +
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
    axis.text = element_blank())

ggsave(plot = plot_k_discrete, filename = paste0('/Users/nm/Desktop/qc/plotk_',country_code,'.pdf'), dpi = 600, height = 6, width = 12)








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



layers <- st_read_parquet('/Users/nm/Downloads/production/outputs/complexity/complexity_DJI_layers.parquet')
blocks <- read_parquet('/Users/nm/Downloads/production/outputs/complexity/complexity_DJI.parquet')

layers2 <- layers %>% st_drop_geometry() %>% 
  filter(block_property == 'building-parcels') %>% 
  group_by(block_id) %>%
  summarize_at(vars(building_count, k_complexity), list(sum, max)) %>%
  ungroup() %>%
  select(block_id, building_count_fn1, k_complexity_fn2)

blocks <- blocks %>%
  left_join(., layers2, by = c('block_id' = 'block_id')) %>%
  mutate(diff_building_count = building_count - building_count_fn1, 
         diff_k_complexity = k_complexity - k_complexity_fn2)

unique(blocks$diff_building_count)
unique(blocks$diff_k_complexity)

lay_viz <- layers %>% filter(block_id == 'DJI.3.1_1_3148')
ggplot() + 
  geom_sf(data = lay_viz  %>% filter(block_property == 'building-parcels'), 
          aes(fill = k_complexity), color = 'black', alpha = .9) + 
  scale_fill_viridis() + 
  geom_sf(data = lay_viz  %>% filter(block_property == 'off-network-streets'),
          color = 'green', fill = 'white', size =1, alpha = .8) + 
  geom_sf(data = lay_viz  %>% filter(block_property == 'on-network-streets'),
          color = 'red', fill = 'white', size =1, alpha = .8) 
lay_viz


lay_viz  %>% filter(block_property == 'on-network-streets') %>%
  st_set_crs(4326) %>%
  st_transform(3395) %>%
  st_length()


# fix connected linestrings







# Paths --------------------------------------------------------------

country_code = c('DJI','COG','ZAF','NGA')

dir_path = '/Users/nm/Downloads/viz/'
iso_code = country_code[1]

metrics_c1 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/run8/complexity_DJI.parquet')
metrics_c2 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/complexity_DJI.parquet')

metrics_c1 <- metrics_c1 %>%
  left_join(.,metrics_c2, by = c('block_id'='block_id') ) %>%
  mutate(in_dif = street_access_length - length_of_internal_street,
         ex_dif = distance_to_nearest_street - distance_to_external_street
  )




names(metrics_c1)

write_csv(metrics_c,'/Users/nm/Desktop/c_DJI.csv')
names(metrics_c)
#blocks <-  st_read_parquet(paste0('/Users/nm/Downloads/production/outputs/blocks/blocks_',iso_code,'.parquet'))
metrics <- read_parquet(paste0('/Users/nm/Downloads/production/outputs/kindex/kindex_',iso_code,'.parquet'))
population <- read_parquet(paste0('/Users/nm/Downloads/production/outputs/population/population_',iso_code,'.parquet'))
#streets <-  st_read('/Users/nm/Downloads/production/outputs/streets/streets_ZAF.parquet')

names(metrics)
names(population)

metrics0 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/run8/complexity_DJI.parquet') %>%
  rename(building_layers0 = building_layers,
         k_complexity0 = k_complexity) %>%
  select(block_id, street_access_length, distance_to_nearest_street, building_layers0,  k_complexity0)

metrics1 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/run9/complexity_DJI.parquet') %>%
  rename(building_layers1 = building_layers,
         k_complexity1 = k_complexity,
         length_of_internal_street1 = length_of_internal_street,
         distance_to_external_street1 = distance_to_external_street) %>%
  select(block_id, length_of_internal_street1, distance_to_external_street1, building_layers1,  k_complexity1)

metrics2 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/complexity_DJI.parquet')
%>%
  rename(building_layers2 = building_layers,
         k_complexity2 = k_complexity,
         length_of_internal_street2 = length_of_internal_street,
         distance_to_external_street2 = distance_to_external_street) %>%
  select(block_id, length_of_internal_street2, distance_to_external_street2, building_layers2,  k_complexity2)

metrics <- metrics0 %>%
  left_join(.,metrics1, by = c('block_id'='block_id')) %>%
  left_join(.,metrics2, by = c('block_id'='block_id')) %>%
  mutate(dif1 = k_complexity1 - k_complexity0 ,
         dif2 = k_complexity2 - k_complexity0 ,
         dif21 = k_complexity2 - k_complexity1 ,
         len21 = length_of_internal_street2 - length_of_internal_street1,
         dist21 = distance_to_external_street2 - distance_to_external_street1)



#   
metrics2 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/run1/complexity_COG.parquet')
names(metrics2)
metrics2 <- metrics2 %>%
  rename(building_layers2 = building_layers,
         k_complexity2 = k_complexity) %>%
  select(block_id, building_layers2,  k_complexity2)

metrics3 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/run2/complexity_COG.parquet')
metrics3 <- metrics3 %>%
  rename(building_layers3 = building_layers,
         k_complexity3 = k_complexity) %>%
  select(block_id, building_layers3,  k_complexity3)

metrics4 <-  read_parquet('/Users/nm/Downloads/production/outputs/complexity/complexity_COG.parquet')
metrics4 <- metrics4 %>%
  rename(building_layers4 = building_layers,
         k_complexity4 = k_complexity) %>%
  select(block_id, building_layers4,  k_complexity4)

metrics <- metrics %>%
  left_join(.,metrics2, by = c('block_id'='block_id')) %>%
  left_join(.,metrics3, by = c('block_id'='block_id')) %>%
  left_join(.,metrics4, by = c('block_id'='block_id')) 

metrics <- metrics %>%
  mutate(dif2 = k_complexity2 - k_complexity ) %>%
  mutate(dif3 = k_complexity3 - k_complexity ) %>%
  mutate(dif4 = k_complexity4 - k_complexity ) %>%
  mutate(dif32 = k_complexity3 - k_complexity2,
         dif42 = k_complexity4 - k_complexity2)

b <- st_read_parquet('/Users/nm/Downloads/production/inputs/buildingpoints/buildings_SYC.parquet')



# +
#   ggsn::scalebar(y.min = st_bbox(data)$ymin - .003, x.min = st_bbox(data)$xmin, 
#                  y.max = st_bbox(data)$ymax, x.max = st_bbox(data)$xmax, location = 'bottomleft',
#                  height = .01, box.fill = c('#333333','#ffffff'),
#                  border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
#                  dist = (country_bbox$xmax - country_bbox$xmin)*0.001*.1, 
#                  transform = FALSE,
#                  dist_unit = "km")
# (country_bbox$xmax - country_bbox$xmin)







library(tidyverse)
library(sf)
library(lwgeom)
library(tidycensus)
library(scales)
library(viridis)
library(DT)
library(shiny)
library(ggplot2)
library(readxl)
library(patchwork)

# Read the .Renviron file (only necessary fi you ran census_api_key()
readRenviron("~/.Renviron")

# Explore ACS application ----------------------------------------------------

# Function to launch a mini Shiny app to look up Census variables
explore_acs_vars <- function () { 
  ui <- basicPage(h2("ACS Variable Search"), 
                  tags$style('#display {height:100px; white-space: pre-wrap;}'),
                  verbatimTextOutput('display', placeholder = TRUE),
                  mainPanel(DT::dataTableOutput(outputId = "acs_table", width = '800px'))
  )
  server <- function(input, output, session) {
    output$acs_table= DT::renderDataTable({ 
      acs5_vars <- acs5_vars 
    }, filter = "top", selection = 'multiple', options = list(columnDefs = list( list(className = "nowrap",width = '100px', targets = c(1,2))), pageLength = 20), server = TRUE) 
    selected_index <- reactive({
      acs5_vars %>% slice(input$acs_table_rows_selected) %>% pull(name)
    })
    output$display = renderPrint({
      s = unique(input$acs_table_rows_selected)
      if (length(s)) {cat(paste0("'",selected_index(),"'",collapse = ","))}
    })
  }
  shinyApp(ui, server)
}

# Census Variables
acs5_vars <- load_variables(year = 2020, dataset = c('acs5'), cache = FALSE) %>% 
  separate(col = 'concept',  into = c('concept_main','concept_part'), sep = c(' BY '), remove = FALSE,extra = "merge") %>%
  mutate(concept_part = case_when(is.na(concept_part) ~ 'TOTAL', TRUE ~ as.character(concept_part)))

explore_acs_vars()






library(osmdata)

aoi_box <- target_area %>% st_transform(4326)

# Download OSM water features 
water <- opq(bbox = st_bbox(aoi_box), memsize = 1e+9) %>%
  add_osm_feature(key = 'water') %>%
  osmdata_sf() 
water_mulitpolygons <- water$osm_multipolygons %>% dplyr::select(osm_id)
water_polygons <- water$osm_polygons %>% dplyr::select(osm_id)
water_lines <- water$osm_lines %>% dplyr::select(osm_id)

# Download OSM waterway features
waterway <- opq(bbox = st_bbox(aoi_box), memsize = 1e+9) %>%
  add_osm_feature(key = 'waterway') %>%
  osmdata_sf() 
waterway_mulitpolygons <- waterway$osm_multipolygons %>% dplyr::select(osm_id)
waterway_polygons <- waterway$osm_polygons %>% dplyr::select(osm_id)
waterway_lines <- waterway$osm_lines %>% dplyr::select(osm_id)
waterway_multilines <- waterway$osm_multilines %>% dplyr::select(osm_id)

# Download OSM coastline features
coastline <- opq(bbox = st_bbox(aoi_box), memsize = 1e+9) %>%
  add_osm_feature(key = 'natural', value = 'coastline') %>%
  osmdata_sf() %>%
  pluck("osm_lines") 

# Parse and combine water linestrings and polygons
water_poly <- rbind(water_mulitpolygons,water_polygons,waterway_mulitpolygons,waterway_polygons) %>%
  st_intersection(.,st_as_sfc(st_bbox(aoi_box))) %>%
  st_transform(4326) %>% dplyr::select(geometry) 
water_line <- rbind(coastline, water_lines,waterway_lines,waterway_multilines) %>%
  st_intersection(.,st_as_sfc(st_bbox(aoi_box))) %>%
  st_transform(4326) %>% dplyr::select(geometry)

# Create natural barrier polygons -----------------------------------------

# Build bounding polygon
print('Build bounding polygon')
crop_box = st_bbox(aoi_box) %>% st_as_sfc()
crop_box = (crop_box - st_centroid(crop_box)) * 0.999 + st_centroid(crop_box)
crop_box <- crop_box %>% st_set_crs(4326)

# Clean up OSM water features
print('Clean up OSM water features')
water_barriers = rbind(water_poly %>% st_as_sf() ) %>% 
  st_union() %>% st_as_sf() %>% 
  sf::st_crop(x = . , y = crop_box) %>%
  st_transform(3857) %>%
  st_collection_extract(. , type = c( "POLYGON")) %>%
  st_cast("POLYGON") %>%
  mutate(area_w = st_area(x)) %>%
  filter(area_w >= units::set_units(1e+07,m^2)) %>%
  smooth(., method = "chaikin") %>%
  st_make_valid()
plot(water_barriers)

# Polygons of water features (no oceans)
print('Polygons of water features (no oceans)')
water_boundaries_poly <- water_barriers %>% st_transform(3857) %>%
  st_buffer(x = ., dist = units::set_units(300,m)) %>%
  smooth(., method = "chaikin") %>%
  st_buffer(x = ., dist = units::set_units(1000,m)) %>%
  st_simplify(x = ., preserveTopology = TRUE, dTolerance = units::set_units(500,m)) %>%
  st_make_valid()
plot(water_boundaries_poly)

# Linestrings of water features
print('Linestrings of water features')
water_boundaries <- water_boundaries_poly %>%
  st_boundary() %>%
  st_collection_extract(. , type = c( "LINESTRING")) %>%
  st_cast(., "LINESTRING") %>%
  rename(geometry = x) %>%
  dplyr::select(geometry)
plot(water_boundaries)

ggplot() +
  geom_sf(data = water_boundaries_poly, color = 'blue') +
  geom_sf(data = water_barriers, color = 'green') 


