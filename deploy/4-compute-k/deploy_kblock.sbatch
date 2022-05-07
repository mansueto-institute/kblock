#!/bin/bash

#SBATCH --job-name=::COUNTRYCODE::
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/job_::COUNTRYCODE::.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/job_::COUNTRYCODE::.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
country_chunk=(::COUNTRYCODE::)

log_file_arg=/project2/bettencourt/mnp/production/jobs/log_::COUNTRYCODE::.log
country_chunk_arg=${country_chunk[@]}
chunk_size_arg=100000
core_count_arg=14
blocks_dir_arg=/project2/bettencourt/mnp/production/outputs/blocks
streets_dir_arg=/project2/bettencourt/mnp/production/outputs/streets
buildings_dir_arg=/project2/bettencourt/mnp/production/inputs/buildingpoints
population_dir_arg=/project2/bettencourt/mnp/production/outputs/population
output_dir_arg=/project2/bettencourt/mnp/production/outputs


python batch_kblock.py --log_file $log_file_arg --country_chunk $country_chunk_arg --chunk_size $chunk_size_arg --core_count $core_count_arg --blocks_dir $blocks_dir_arg --streets_dir $streets_dir_arg --buildings_dir $buildings_dir_arg --population_dir $population_dir_arg --output_dir $output_dir_arg


