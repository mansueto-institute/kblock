#!/bin/bash

max_number_of_jobs=10
cnetid=nmarchio

template_dir=/project2/bettencourt/mnp/production/jobs/filled_templates
sbatch_file=/project2/bettencourt/mnp/production/repos/kblock/deploy/4-compute-k/deploy_kblock.sbatch

input_dir=/project2/bettencourt/mnp/production/inputs/buildingpoints/
input_prefix=buildings_
input_suffix=.parquet

output_dir=/project2/bettencourt/mnp/production/outputs/complexity/
output_prefix=kindex_
output_suffix=.parquet

valid_codes=(NGA ZAF ETH SDN TZA KEN COD UGA MOZ BFA ZWE AGO MLI GHA NER MDG TCD ZMB MWI CMR CIV SEN GIN BEN RWA SSD TGO SOM NAM BDI MRT BWA ERI CAF COG SLE LSO LBR SWZ GAB GMB MUS GNB GNQ CPV COM ESH DJI STP SYC)

mkdir -p "$template_dir"

input_list=()
for i in $input_dir$input_prefix*$input_suffix; do
    i="${i/$input_dir/}"
    i="${i/$input_prefix/}"
    i="${i/$input_suffix/}"
    if [[ "${valid_codes[@]}" =~ "$i" ]]; then
        input_list+=("$i")
    fi
done

output_list=()
for i in $output_dir$output_prefix*$output_suffix; do
    i="${i/$output_dir/}"
    i="${i/$output_prefix/}"
    i="${i/$output_suffix/}"
    if [[ "${valid_codes[@]}" =~ "$i" ]]; then
        output_list+=("$i")
    fi
done

residual_list=()
residual_list=(`echo ${input_list[@]} ${output_list[@]} | tr ' ' '\n' | sort | uniq -u`)

country_list=()
country_list=(`echo ${residual_list[@]} | tr ' ' '\n' | sort | cut -d'/' -f 4 | uniq`)
country_list=(${country_list[@]:0:$max_number_of_jobs})

job_list=()
for job in ${country_list[@]}; do
check_job=$(squeue --user="$cnetid" --name="$job" --Format=name --noheader)
check_job="$(echo -e "${check_job}" | tr -d '[:space:]')"
if [ "$job" "=" "$check_job" ]; then
  echo "$job job is running"
else
  echo "$job job is not running"
  job_list+=("$job")
fi
done
echo "Queue: ${job_list[@]}"

for countrycode in "${job_list[@]}"; do
echo "Deploy $countrycode"
< "$sbatch_file" sed -e "s/::COUNTRYCODE::/${countrycode}/g" > "$template_dir"/"${countrycode}_job.sbatch"
sbatch "$template_dir"/"${countrycode}_job.sbatch"
done

