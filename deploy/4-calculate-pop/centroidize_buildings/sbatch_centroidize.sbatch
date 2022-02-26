#!/bin/bash

#SBATCH --job-name=centroidize
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/centroidize.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/centroidize.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=merrittsmith@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

python batch_centroidize.py --top_building_dir /project2/bettencourt/mnp/production/inputs/buildings
