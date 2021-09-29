#!/bin/bash
#SBATCH --job-name=fiona_error
#SBATCH --output=fiona_error.out
#SBATCH --error=fiona_error.err
#SBATCH --time=0:30:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-user=alexfeistritzer@uchicago.edu
#SBATCH --mail-type=ALL

module purge
module load python

base_path='/project2/bettencourt/mnp/prclz/data'
full_name='TZA.10.1.4_1'
aoi=${full_name%%.*}

python $base_path/parallel/utils/block_summary.py --aoi_path $base_path/blocks/Africa/${aoi}/blocks_${full_name}.csv --landscan_path $base_path/LandScan_Global_2018/raw_tif/ls_2018.tif --buildings_dir $base_path/buildings/Africa/${aoi}/ --blocks_dir $base_path/blocks/Africa/${aoi} --gadm_dir $base_path/GADM/${aoi}/ --summary_out_path ${full_name}_summary.geojson