#!/bin/bash

# Parameters
max_number_of_jobs=10
cnetid=nmarchio

# Directory to place filled SBATCH templates (exclude trailing slash)
template_dir=/project2/bettencourt/mnp/production/jobs/filled_templates
# SBATCH template
sbatch_file=/project2/bettencourt/mnp/production/repos/kblock/deploy/4-compute-k/deploy_kblock.sbatch

# Path of input files and prefix
input_start=/project2/bettencourt/mnp/production/inputs/buildingpoints/buildings_
# Path of input files and suffix
input_end=.parquet

# Path of output files and prefix
output_start=/project2/bettencourt/mnp/production/outputs/kindex/kindex_
# Path of output files and suffix
output_end=.parquet

mkdir -p "$template_dir"

# Add input files to array
input_list=()
for i in "$input_start"*"$input_end"; do
    i="${i/$input_start/}"
    i="${i/"$input_end"/}"
    input_list+=("$i")
done
# echo "${input_list[@]}"

# Add output files to array
output_list=()
for i in "$output_start"*"$output_end"; do
    i="${i/$output_start/}"
    i="${i/"$output_end"/}"
    output_list+=("$i")
done
# echo "${output_list[@]}"

# List residual countries between output and input files
residual_list=()
residual_list=(`echo ${input_list[@]} ${output_list[@]} | tr ' ' '\n' | sort | uniq -u`)

# Print residual countries and limit to max jobs
country_list=()
country_list=(`echo ${residual_list[@]} | tr ' ' '\n' | sort | cut -d'/' -f 4 | uniq`)
printf '%s\n' "${country_list[@]}"
country_list=(${country_list[@]:0:$max_number_of_jobs})

# Remove countries currently deployed
job_list=()
for job in ${country_list[@]}; do
check_job=$(squeue --user="$cnetid" --name="$job" --Format=name --noheader)
if [ "$job" "==" "$check_job" ]; then
  echo "$job job is running"
else
  job_list+=("$job")
fi
done
echo "${job_list[@]}"

# Submit a jobs
for countrycode in "${job_list[@]}"; do
echo "$countrycode"
< "$sbatch_file" sed -e "s/::COUNTRYCODE::/${countrycode}/g" > "$template_dir"/"${countrycode}_job.sbatch"
sbatch "$template_dir"/"${countrycode}_job.sbatch"
done



