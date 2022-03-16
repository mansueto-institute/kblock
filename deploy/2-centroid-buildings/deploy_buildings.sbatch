#!/bin/bash

#SBATCH --job-name=building-points
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/building_points_job.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/building_points_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
#source activate pygisparallel

country_list=(SLE AGO NGA TZA ZAF ZMB BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER SDN SEN SOM SSD STP SWZ SYC TCD TGO ZWE UGA RWA)

log_file_arg=/project2/bettencourt/mnp/production/jobs/log_building_points.log
country_chunk_arg=${country_list[@]}
buildings_dir_arg=/project2/bettencourt/mnp/production/inputs/buildings
gadm_dir_arg=/project2/bettencourt/mnp/production/inputs/gadm
output_dir_arg=/project2/bettencourt/mnp/production/inputs/buildingpoints

python batch_buildings.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --buildings_dir $buildings_dir_arg --gadm_dir $gadm_dir_arg  --output_dir $output_dir_arg

