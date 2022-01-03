#!/bin/bash

#SBATCH --job-name=build_block
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=54000
#SBATCH --output=/project2/bettencourt/mnp/analytics/deployments/build_block_job.out
#SBATCH --error=/project2/bettencourt/mnp/analytics/deployments/build_block_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=16:00:00
#SBATCH --account=pi-bettencourt

country_list=(STP SYC GNQ COM MUS DJI GNB CPV ERI GAB MRT SWZ COG LBR RWA SLE BDI BEN TGO NAM CIV NER AGO BFA SDN GHA GMB TCD CAF ETH BWA SOM SSD GIN MWI LSO MLI ZWE KEN MOZ ESH CMR ZMB MDG ZAF UGA COD NGA TZA SEN)

log_file_arg=/project2/bettencourt/mnp/analytics/deployments/log_build_blocks.log
country_chunk_arg=${country_list[@]}
osm_dir_arg=/project2/bettencourt/mnp/analytics/data/geofabrik
gadm_dir_arg=/project2/bettencourt/mnp/analytics/data/gadm
output_dir_arg=/project2/bettencourt/mnp/analytics/outputs

python batch_etl.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg  --output_dir $output_dir_arg
