source("./boundaries.R")

# Visualize zoom area -----------------------------------------------------

road_color = '#ffffff'
grey2 <- c('#414141','#777777')
kdist = max(as.integer(zoom_zone$k_complexity_groups))
colorhexes <- colorRampPalette(c('#93328E','#CF3F80','#F7626B','#FF925A','#FFC556','#F9F871'))(length(unique(zoom_zone$k_complexity_groups))-2)

width = st_distance(st_sf(geom = st_sfc(st_point(c(st_bbox(zoom_zone)$xmin, st_bbox(zoom_zone)$ymin)),
                                        st_point(c(st_bbox(zoom_zone)$xmax, st_bbox(zoom_zone)$ymin)), crs = 4326)))[2] %>%  drop_units()
height = st_distance(st_sf(geom = st_sfc(st_point(c(st_bbox(zoom_zone)$xmin, st_bbox(zoom_zone)$ymin)),
                                         st_point(c(st_bbox(zoom_zone)$xmin, st_bbox(zoom_zone)$ymax)), crs = 4326)))[2] %>%  drop_units()
width_tenth = round((width*.2)/1000,-1)
if (width_tenth < 1) {
  width_tenth = round((width*.2)/1000,0)
}
height_decdegs = abs(unname(st_bbox(zoom_zone)$ymax) - unname(st_bbox(zoom_zone)$ymin))

(plot_k_discrete <- ggplot() +
    geom_sf(data = zoom_zone, aes(fill = as.factor(k_complexity_groups)), color = road_color, size = .0075) +   
    # geom_sf(data = water_poly, fill = 'white', color = 'white', size = .1) +
    geom_sf(data = water_line, fill = 'white', color = 'white', size = .1) +
    scale_fill_manual(values = c(grey2,colorhexes), name = 'k complexity') + 
    labs(caption = paste0('Population-weighted average k complexity: ',zoom_zone %>% st_drop_geometry() %>% summarise(wm_var = weighted.mean(as.integer(k_complexity), landscan_population)) %>% pull() %>% round(.,2))) +
    guides(color = guide_legend(nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1),
           fill =  guide_legend(nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1)) +
    theme_void() + theme(plot.caption = element_text(size = 11, hjust = .5, vjust = 20, margin=margin(0,0,0,0)),
                         text = element_text(color = "#333333"),
                         #legend.position = c(1.05,.7),
                         legend.position = 'bottom',
                         legend.spacing.x = unit(1, 'pt'),
                         #legend.key.height = unit(10, 'pt'), 
                         #legend.key.width = unit(10, 'pt'),
                         legend.text = element_text(size = 10),
                         panel.border = element_blank(),
                         panel.background = element_blank(),
                         plot.margin=unit(c(t=0,r=10,b=0,l=0), "pt"),
                         legend.title = element_blank(),
                         axis.text = element_blank()) +
    ggsn::scalebar(y.min = st_bbox(zoom_zone)$ymin - (height_decdegs*.03), 
                   x.min = st_bbox(zoom_zone)$xmin, 
                   y.max = st_bbox(zoom_zone)$ymax, 
                   x.max = st_bbox(zoom_zone)$xmax, 
                   location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = width_tenth/2, dist_unit = "km", transform = TRUE, model = "WGS84") 
)


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
          plot.margin=unit(c(t=0,r=5,b=0,l=5), "pt"),
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
    # geom_sf(data = water_poly, fill = 'white', color = 'white', size = .1) +
    geom_sf(data = water_line, fill = 'white', color = 'white', size = .1) +
    labs(subtitle = "Population",
         caption = paste0('Total population in area: ',comma(sum(zoom_zone$landscan_population)),'  |  ',
                          'Average block population: ',comma(round(mean(zoom_zone$landscan_population),2)))) + 
    # scale_fill_viridis(name = 'Population', oob = scales::squish, limits= c(1, max(zoom_zone$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) + 
    # scale_color_viridis(name = 'Population', oob = scales::squish, limits= c(1, max(zoom_zone$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) +
    scale_fill_distiller(palette = 'Spectral', name = 'Population', oob = scales::squish, limits= c(1, max(zoom_zone$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) + 
    scale_color_distiller(palette = 'Spectral', oob = scales::squish, limits= c(1, max(zoom_zone$landscan_population_log )), breaks= c(1,2,3,4,5,6,7), labels = c('0',"100","1K","10K","100K","1M","10M")) +
    theme_void() + 
    #guides(color = guide_legend(title.position="top", title.hjust = 0.5),
    #       fill = guide_legend(title.position="top", title.hjust = 0.5)) +
    theme(#plot.subtitle = element_text(size = 14, face="bold", vjust = -4, hjust=.5),axis.text = element_blank(),
      plot.caption = element_text(size = 10, hjust = .5, vjust = 5),
      #legend.position = c(1.05,.7),
      #legend.position = 'bottom',
      #legend.position = 'top',
      #legend.position = c(.5, .01),
      legend.position = c(.5, 1),
      legend.direction = "horizontal",
      legend.key.width=unit(40,"pt"),
      legend.key.height=unit(5,"pt"),
      plot.margin = unit(c(t=15,r=0,b=0,l=0), "pt"),
      plot.subtitle = element_text(size = 11, face="bold", vjust = 5, hjust = .5),
      legend.title = element_blank(),
      #legend.title = element_text(face="bold", hjust = .5),
      text = element_text(color = "#333333")) +
    ggsn::scalebar(y.min = st_bbox(zoom_zone)$ymin - (height_decdegs*.04), 
                   x.min = st_bbox(zoom_zone)$xmin, y.max = st_bbox(zoom_zone)$ymax, x.max = st_bbox(zoom_zone)$xmax, 
                   location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = width_tenth/2, dist_unit = "km", transform = TRUE, model = "WGS84") )


sd_int = log10(mean(zoom_zone$landscan_pop_density_hectare, na.rm = TRUE) + sd(zoom_zone$landscan_pop_density_hectare, na.rm = TRUE)*1)
(plot_popdensity_log <- ggplot() +
    geom_sf(data = zoom_zone,
            aes(fill = landscan_pop_density_hectare_log), 
            color = 'white', size= .0075, alpha = .8) +
    # geom_sf(data = water_poly, fill = 'white', color = 'white', size = .1) +
    geom_sf(data = water_line, fill = 'white', color = 'white', size = .1) +
    labs(subtitle = "Population per hectare",
         caption = paste0('Weighted average population density: ',
                          zoom_zone %>% st_drop_geometry() %>% summarize(pop_dense = weighted.mean(landscan_pop_density_hectare, landscan_population) ) %>% pull() %>% round(.,0),
                          ' people per hectare','\n 1 hectare = 10k m^2 = 1.4 soccer fields = 2.2 Manhattan city blocks')) +
    # scale_fill_viridis(name = 'Population\nper hectare', oob = scales::squish, limits= c(1, 3), 
    #                    breaks= c(1,2,3,4,5,6,7), 
    #                    labels = c('0',"100","1K","10K","100K","1M","10M")) +
    # scale_color_viridis(name = 'Population\nper hectare', oob = scales::squish, limits= c(1, 3), 
    #                     breaks= c(1,2,3,4,5,6,7), 
    #                     labels = c('0',"100","1K","10K","100K","1M","10M")) +
    scale_fill_distiller(direction = -1, palette = 'Spectral', name = 'Population\nper hectare', oob = scales::squish, limits= c(1, 3), 
                         breaks= c(1,2,3,4,5,6,7), 
                         labels = c('0',"100","1K","10K","100K","1M","10M")) +
    scale_color_distiller(direction = -1, palette = 'Spectral', name = 'Population\nper hectare', oob = scales::squish, limits= c(1, 3), 
                          breaks= c(1,2,3,4,5,6,7), 
                          labels = c('0',"100","1K","10K","100K","1M","10M")) +
    theme_void() + 
    theme(#plot.subtitle = element_text(size = 14, face="bold", vjust = -4, hjust=.5),
      #plot.caption = element_text(size = 10, hjust = .5, vjust = 25),
      plot.caption = element_text(size = 10, hjust = .5, vjust = 5),
      #legend.position = c(1.05,.7),
      #legend.position = 'top',
      legend.position = c(.5, 1),
      legend.direction = "horizontal",
      #legend.position = c(.5, .01),
      legend.key.width=unit(40,"pt"),
      legend.key.height=unit(5,"pt"),
      plot.margin=unit(c(t=15,r=0,b=0,l=0), "pt"),
      #axis.text = element_blank(),
      plot.subtitle = element_text(size = 11, face="bold", vjust = 5, hjust = .5),
      legend.title = element_blank(),
      #legend.title = element_text(face="bold", hjust = .5),
      text = element_text(color = "#333333")) +
    ggsn::scalebar(y.min = st_bbox(zoom_zone)$ymin - (height_decdegs*.04), 
                   x.min = st_bbox(zoom_zone)$xmin, y.max = st_bbox(zoom_zone)$ymax, x.max = st_bbox(zoom_zone)$xmax, 
                   location = 'bottomleft',
                   height = .01, box.fill = c('#333333','#ffffff'),
                   border.size = .4, st.color = '#333333', st.size = 2.5, box.color = '#333333',
                   dist = width_tenth/2, dist_unit = "km", transform = TRUE, model = "WGS84") )

if (ghsl_delin == TRUE) {
  full_city = '_city'
} else {
  full_city = ''
}


(plot_k <- plot_k_discrete + bar_k_distrib +#+ plot_spacer() 
    plot_layout(widths = c(1,.8))  + # , height = c(1.5, 1) .01 ,
    plot_annotation(#title = paste0(city_name,', ', country_name),
      subtitle = paste0(city_label,', ', country_label),
      theme = theme(#plot.title = element_text(face="bold", size = 18, vjust = -2, hjust = .5),
        plot.subtitle = element_text(face="bold", size = 13, vjust = -7, hjust = .5))))
ggsave(plot = plot_k, filename = paste0(dir_path,city_label,'_k',full_city,'.pdf'), height = (12*(height/width))/2+1.5, width = 12)
ggsave(plot = plot_k, filename = paste0(dir_path,city_label,'_k',full_city,'.png'), dpi = 300, height = (12*(height/width))/2+1.5, width = 12)

(plot_pop <- plot_populaton + plot_popdensity_log +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(#title = paste0(city_name,', ', country_name ),
      subtitle = paste0(city_label,', ', country_label),
      theme = theme(#plot.title = element_text(face="bold", size = 18, vjust = -2, hjust = .5),
        plot.subtitle = element_text(face="bold", size = 13, vjust = -5, hjust = .5)))
)
ggsave(plot = plot_pop, filename = paste0(dir_path,city_label,'_pop',full_city,'.pdf'), height = (12*(height/width))/2+1.5, width = 12)
ggsave(plot = plot_pop, filename = paste0(dir_path,city_label,'_pop',full_city,'.png'), dpi = 300, height = (12*(height/width))/2+1.5, width = 12)


# Visualize country -------------------------------------------------------

pdf(paste0('./','Angola','.pdf') )
plot(block_AGO %>% select(block_id), main="", lwd=.01)
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
                        ' people per hectare','\n  1 hectare = 10k m^2 = 1.4 soccer fields = 2.2 Manhattan city blocks')) +
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
