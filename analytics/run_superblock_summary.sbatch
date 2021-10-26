#!/bin/bash
#SBATCH --job-name=kblock_summary
#SBATCH --output=parallel_block_summary.out
#SBATCH --error=parallel_block_summary.err
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

python mnp-analysis/analytics2/utils/block_summary.py --superblock_path ../dev/data/kblock_freetown.geojson --landscan_path ../dev/data/ls_2019.tif --buildings_path /project2/bettencourt/mnp/analytics/data/buildings/SLE/buildings_SLE.4.2.1_1.geojson --summary_out_path SLE.4.2.1_1_summary.geojson
