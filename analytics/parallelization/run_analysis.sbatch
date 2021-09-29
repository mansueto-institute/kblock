#!/bin/bash

base_path=$1 # '/project2/bettencourt/mnp/prclz/data'
full_name=$2 # 'TZA.30.1.40_1'
aoi=${full_name%%.*}

#fn=$2 # 'blocks_TZA.30.1.40_1.csv'
#not_blocks=${fn##"blocks_"} # 'TZA.30.1.40_1.csv'
#aoi=${not_blocks%%.*} # 'TZA'
#full_name=${not_blocks%%".csv"} # 'TZA.30.1.40_1'

python $base_path/parallel/utils/block_summary.py --aoi_path $base_path/blocks/Africa/${aoi}/blocks_${full_name}.csv --landscan_path $base_path/LandScan_Global_2018/raw_tif/ls_2018.tif --buildings_dir $base_path/buildings/Africa/${aoi}/ --blocks_dir $base_path/blocks/Africa/${aoi} --gadm_dir $base_path/GADM/${aoi}/ --summary_out_path ${full_name}_summary.geojson