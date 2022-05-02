


# Read in cities from GHSL

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



# country_iso = c('SYC'),
# country_name = c('Seychelles'),
# city_name = c("Victoria"),
# ghsl_pop15 = c(9462)



names(city_data)


#Make markdown PDF creation script
# Figure out way to read parquet into R

# Water data to remove water blocks

