


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








# +
#   ggsn::scalebar(y.min = st_bbox(data)$ymin - .003, x.min = st_bbox(data)$xmin, 
#                  y.max = st_bbox(data)$ymax, x.max = st_bbox(data)$xmax, location = 'bottomleft',
#                  height = .01, box.fill = c('#333333','#ffffff'),
#                  border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
#                  dist = (country_bbox$xmax - country_bbox$xmin)*0.001*.1, 
#                  transform = FALSE,
#                  dist_unit = "km")
# (country_bbox$xmax - country_bbox$xmin)


