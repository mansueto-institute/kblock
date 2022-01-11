#!/bin/bash

#SBATCH --job-name=build_block
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=54000
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/build_block_job.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/build_block_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=24:00:00
#SBATCH --account=pi-bettencourt


# country_list=(STP SYC GNQ COM MUS DJI GNB CPV ERI GAB MRT SWZ COG LBR RWA SLE BDI BEN TGO)
# country_list=(NAM CIV NER AGO BFA SDN GHA GMB TCD CAF ETH BWA SOM SSD GIN MWI LSO MLI ZWE KEN MOZ ESH CMR ZMB MDG ZAF UGA COD NGA TZA SEN)
country_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA ZAF ZMB ZWE UGA RWA)

log_file_arg=/project2/bettencourt/mnp/production/jobs/log_build_blocks.log
country_chunk_arg=${country_list[@]}
osm_dir_arg=/project2/bettencourt/mnp/production/inputs/osm
gadm_dir_arg=/project2/bettencourt/mnp/production/inputs/gadm
output_dir_arg=/project2/bettencourt/mnp/production/outputs

python batch_blocks.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg  --output_dir $output_dir_arg
