#!/bin/bash

#SBATCH --job-name=kblock-sle
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/analytics/deployments/sle_job.out
#SBATCH --error=/project2/bettencourt/mnp/analytics/deployments/sle_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=12:00:00
#SBATCH --account=pi-bettencourt

cd /project2/bettencourt/mnp/analytics/scripts/kblock/dev/compute/

module load python/anaconda-2020.02
source activate pygeos_env

log_file_arg=/project2/bettencourt/mnp/analytics/deployments/log_output_sle.log
country_code_arg=SLE
country_code_file_arg=/project2/bettencourt/mnp/analytics/scripts/kblock/dev/compute/country_codes.csv
gadm_parent_dir_arg=/project2/bettencourt/mnp/analytics/data/gadm
streets_parent_dir_arg=/project2/bettencourt/mnp/analytics/data/geofabrik
building_parent_dir_arg=/project2/bettencourt/mnp/analytics/data/buildings
output_dir_arg=/project2/bettencourt/mnp/analytics/outputs

python batch_kblock.py --log_file $log_file_arg --country_code $country_code_arg --country_code_file $country_code_file_arg --gadm_parent_dir $gadm_parent_dir_arg --streets_parent_dir $streets_parent_dir_arg --building_parent_dir $building_parent_dir_arg --output_dir $output_dir_arg
