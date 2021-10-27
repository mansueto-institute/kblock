#!/bin/bash
#SBATCH --job-name=kblock_summary
#SBATCH --output=kblock_summary.out
#SBATCH --error=kblock_summary.err
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

python utils/block_summary.py --superblock ../dev/data/kblock_freetown.geojson --landscan_path /project2/bettencourt/mnp/analytics/data/population/ls_2019.tif --superblock_buildings /project2/bettencourt/mnp/analytics/data/buildings/SLE/buildings_SLE.4.2.1_1.geojson --summary_out_path SLE.4.2.1_1_summary.geojson
