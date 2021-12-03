#!/bin/bash

#SBATCH --job-name=kblock-malaria-nga
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/analytics/deployments/malaria_nga_job.out
#SBATCH --error=/project2/bettencourt/mnp/analytics/deployments/malaria_nga_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=12:00:00
#SBATCH --account=pi-bettencourt

cd /project2/bettencourt/mnp/analytics/scripts/kblock/dev/compute/

gadm_list=(NGA.20.7_1 NGA.20.9_1 NGA.20.10_1 NGA.20.12_1 NGA.20.17_1 NGA.20.18_1 NGA.20.21_1 NGA.20.25_1 NGA.20.28_1 NGA.20.30_1 NGA.20.31_1 NGA.20.38_1 NGA.20.39_1 NGA.20.42_1 NGA.20.43_1 NGA.28.15_1 NGA.30.3_1 NGA.30.4_1 NGA.30.21_1 NGA.30.22_1 NGA.30.23_1 NGA.31.1_1 NGA.31.2_1 NGA.31.5_1 NGA.31.6_1 NGA.31.7_1 NGA.31.8_1 NGA.31.9_1 NGA.31.10_1 NGA.31.12_1 NGA.31.14_1 NGA.31.20_1 NGA.31.25_1 NGA.31.26_1)

log_file_arg=/project2/bettencourt/mnp/analytics/deployments/log_output_malaria_nga.log
country_code_arg=NGA
country_code_file_arg=/project2/bettencourt/mnp/analytics/scripts/kblock/dev/compute/country_codes.csv
gadm_parent_dir_arg=/project2/bettencourt/mnp/analytics/data/gadm
gadm_chunk_arg=${gadm_list[@]}
streets_parent_dir_arg=/project2/bettencourt/mnp/analytics/data/geofabrik
building_parent_dir_arg=/project2/bettencourt/mnp/analytics/data/buildings
output_dir_arg=/project2/bettencourt/mnp/analytics/outputs

python batch_kblock.py --log_file $log_file_arg --country_code $country_code_arg --country_code_file $country_code_file_arg --gadm_chunk $gadm_chunk_arg --gadm_parent_dir $gadm_parent_dir_arg --streets_parent_dir $streets_parent_dir_arg --building_parent_dir $building_parent_dir_arg --output_dir $output_dir_arg
