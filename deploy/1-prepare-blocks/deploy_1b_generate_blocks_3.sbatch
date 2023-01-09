#!/bin/bash

#SBATCH --job-name=deploy_1b_generate_blocks
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/update/jobs/deploy_1b_generate_blocks_3.out
#SBATCH --error=/project2/bettencourt/mnp/update/jobs/deploy_1b_generate_blocks_3.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
source activate geospatial

iso_list=(NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)

working_directory=/project2/bettencourt/mnp/update

log_file_arg=$working_directory/jobs/deploy_1b_generate_blocks_3.log
country_chunk_arg=${iso_list[@]}
osm_dir_arg=$working_directory/inputs/osm/parquet
gadm_dir_arg=$working_directory/inputs/gadm/parquet
output_dir_arg=$working_directory/outputs

python batch_1b_generate_blocks.py --log_file $log_file_arg --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg  --output_dir $output_dir_arg

