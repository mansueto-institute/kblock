#!/bin/bash
#SBATCH --job-name=parallel_block_summary_test
#SBATCH --output=parallel_block_summary.out
#SBATCH --error=parallel_block_summary.err
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

python mnp-analysis/analytics2/utils/block_summary.py --aoi_path /project2/bettencourt/mnp/prclz/data/blocks/Africa/SLE/blocks_SLE.4.2.1_1.csv --landscan_path /project2/bettencourt/mnp/prclz/data/LandScan_Global_2018/raw_tif/ls_2018.tif --buildings_dir /project2/bettencourt/mnp/prclz/data/buildings/Africa/SLE/ --blocks_dir /project2/bettencourt/mnp/prclz/data/blocks/Africa/SLE --gadm_dir /project2/bettencourt/mnp/prclz/data/GADM/SLE/ --summary_out_path SLE.4.2.1_1_summary.geojson
