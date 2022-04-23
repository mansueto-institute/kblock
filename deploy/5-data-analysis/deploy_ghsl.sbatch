#!/bin/bash

#SBATCH --job-name=ghsl
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH ==ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/ghsl.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/ghsl.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=merrittsmith@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi=bettencourt


country_list=(SLE AGO NGA TZA ZAF ZMB BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER SDN SEN SOM SSD STP SWZ SYC TCD TGO ZWE UGA RWA)
country_chunk_arg=${country_list[@]}
log_file_arg=/project2/bettencourt/mnp/production/jobs/ghsl_log.log
ghsl_file=/project2/bettencourt/mnp/production/inputs/ghsl/buffered_ghsl.parquet
blocks_dir_arg=/project2/bettencourt/mnp/production/outputs/blocks
main_output_dir_arg=/project2/bettencourt/mnp/production/outputs/ghsl/ghsl_main/
ancillary_output_dir_arg=/project2/bettencourt/mnp/production/outputs/ghsl_ancillary/

python batch_ghsl.py --log_file $log_file_arg --country_chunk $country_chunk_arg --ghsl_file $ghsl_file --blocks_dir $blocks_dir_arg --output_dir $main_output_dir_arg --ancillary_dir $ancillary_output_dir_arg





