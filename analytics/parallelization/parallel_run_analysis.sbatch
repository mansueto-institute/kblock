#!/bin/bash
#SBATCH --job-name=parallel_block_summary
#SBATCH --output=parallel_block_summary.out
#SBATCH --error=parallel_block_summary.err
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-user=alexfeistritzer@uchicago.edu
#SBATCH --mail-type=ALL

# requires file `run_analysis.sbatch` to be in the same directory

module purge
module load parallel
module load python
ulimit -u 10000

# define srun and parallel commands
srun="srun --exclusive -N1 -n1"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog runtask.log --resume"

# base_path is needed by run_analysis.sbatch
# search_path gives location of the `block` .csv files
base_path="/project2/bettencourt/mnp/prclz/data"
search_path="${base_path}/blocks/Africa/TZA"

### for later
#aoi=$1
#search_path="${base_path}/blocks/Africa/${aoi}"
###

### find the completed gadms:
# first 2 commands: find all .geojson files, remove the file locations
# ``while`` command: remove the _summary.geojson and _summary-buildings.geojson from the filename
# last 3 commands: remove duplicate filenames
completed_gadms=($(find . -type f -name '*.geojson' -maxdepth 1 -print | xargs -n 1 basename | while read fname; do echo ${fname%%_s*}; done | tr ' ' '\n' | sort -u | tr '\n' ' '))


### find all possible gadms:
# first while loop removes 'blocks_' from filename
# second removes '.csv'
all_gadms=($(find $search_path -type f -name '*.csv' -print | xargs -n 1 basename | while read fname; do echo ${fname##"blocks_"}; done | while read fname2; do echo ${fname2%%".csv"}; done))

### find the difference between the two and pipe to parallel
# pipes over e.g. 'TZA.1.2.3_1'
echo ${completed_gadms[@]} ${all_gadms[@]} | tr ' ' '\n' | sort | uniq -u | $parallel "$srun ./run_analysis.sbatch $base_path {}"

### we know this works:
#find $search_path -type f -name '*.csv' -print | xargs -n 1 basename | $parallel "$srun ./run_analysis.sbatch $base_path {}"


### when completed, print out the empty files corresponding to missing building files
# line 1 finds all empty files, removes path names, saves to file
# line 2 removes the filename from this list, couldn't fix it in L1
find . -name '*summary.geojson' -maxdepth 1 -size 0c | xargs -n 1 basename > 'missing_blocks.txt'
grep -v 'missing_blocks.txt' missing_blocks.txt > tmpfile && mv tmpfile missing_blocks.txt 