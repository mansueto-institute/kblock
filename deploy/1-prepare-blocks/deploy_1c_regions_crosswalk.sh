#!/bin/bash


# iso_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA ZAF ZMB ZWE UGA RWA)

#(GIN GNB LBR SLE)
iso_list=(LBR SLE)

working_directory=/Users/nm/Downloads/update



log_file_arg=$working_directory/jobs/deploy_1c_regions_crosswalk.log
country_chunk_arg=${iso_list[@]}
africapolis_file_arg=$working_directory/inputs/africapolis/AFRICAPOLIS2020.gpkg
ghsl_file_arg=$working_directory/inputs/ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg
blocks_dir_arg=$working_directory/outputs/blocks
output_dir_arg=$working_directory/outputs

python batch_1c_regions_crosswalk.py --log_file $log_file_arg --country_chunk $country_chunk_arg --africapolis_file $africapolis_file_arg --ghsl_file $ghsl_file_arg --blocks_dir $blocks_dir_arg  --output_dir $output_dir_arg
