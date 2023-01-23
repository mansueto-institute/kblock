#!/bin/bash

#SBATCH --job-name=::COUNTRYCODE::
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/update/jobs/complexity_jobs/::COUNTRYCODE::.out
#SBATCH --error=/project2/bettencourt/mnp/update/jobs/complexity_jobs/::COUNTRYCODE::.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt


country_chunk=(::COUNTRYCODE::)

working_directory=/project2/bettencourt/mnp/update

mkdir -p /project2/bettencourt/mnp/update/jobs/complexity_jobs

log_file_arg=$working_directory/jobs/complexity_jobs/::COUNTRYCODE::.log
country_chunk_arg=${country_chunk[@]}
chunk_size_arg=40000
core_count_arg=10
blocks_dir_arg=$working_directory/outputs/blocks
streets_dir_arg=$working_directory/outputs/streets
buildings_dir_arg=$working_directory/outputs/buildings/points
output_dir_arg=$working_directory/outputs

python /project2/bettencourt/mnp/update/repo/kblock/kblock/batch_4_compute_k.py --log_file $log_file_arg --country_chunk $country_chunk_arg --chunk_size $chunk_size_arg --core_count $core_count_arg --blocks_dir $blocks_dir_arg --streets_dir $streets_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg
