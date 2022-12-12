#!/bin/bash

while getopts ":o:" opt; do
  case $opt in
    o) outpath="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

cwd=$(pwd)

mkdir -p $outpath
mkdir -p $outpath/osm/pbf
mkdir -p $outpath/osm/geojson
mkdir -p $outpath/osm/parquet
mkdir -p $outpath/gadm/shp

geofabrik_list=(angola burundi benin burkina-faso botswana central-african-republic ivory-coast cameroon congo-democratic-republic congo-brazzaville comores cape-verde djibouti eritrea morocco ethiopia gabon ghana guinea senegal-and-gambia guinea-bissau equatorial-guinea kenya liberia lesotho madagascar mali mozambique mauritania mauritius malawi namibia niger nigeria rwanda sudan sierra-leone somalia south-sudan sao-tome-and-principe swaziland seychelles chad togo tanzania uganda south-africa zambia zimbabwe)
iso_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)

# Create array of countries that haven't been downloaded from https://www.geofabrik.de/
osm_finished=()
shopt -s nullglob
for i in $outpath/osm/pbf/*-latest.osm.pbf; do
    osm_finished+=("$i")
done
osm_finished=("${osm_finished[@]//$outpath\/osm\/pbf\//}")
osm_finished=("${osm_finished[@]//-latest.osm.pbf/}")
echo "${osm_finished[@]}"

residual_osm_list=()
residual_osm_list=(`echo ${geofabrik_list[@]} ${osm_finished[@]} | tr ' ' '\n' | sort | uniq -u`)
echo "${residual_osm_list[@]}"

# Download countries from Geofabrik and extract data from PBF to GeoJSON
for country in ${residual_osm_list[@]}; do
echo "$country"
    wget -O $outpath/osm/pbf/${country}-latest.osm.pbf https://download.geofabrik.de/africa/${country}-latest.osm.pbf
    osmium export $outpath/osm/pbf/${country}-latest.osm.pbf -o $outpath/osm/geojson/${country}-latest-linestring.geojson -c $cwd/osm-linestring-config.json --overwrite --geometry-types=linestring
    osmium export $outpath/osm/pbf/${country}-latest.osm.pbf -o $outpath/osm/geojson/${country}-latest-polygon.geojson -c $cwd/osm-water-config.json --overwrite --geometry-types=polygon
done

# Create array of countries that haven't been downloaded from https://gadm.org/
gadm_finished=()
shopt -s nullglob
for i in $outpath/gadm/shp/*; do
    gadm_finished+=("$i")
done
gadm_finished=("${gadm_finished[@]//$outpath\/gadm\/shp\//}")
echo "${gadm_finished[@]}"

residual_gadm_list=()
residual_gadm_list=(`echo ${iso_list[@]} ${gadm_finished[@]} | tr ' ' '\n' | sort | uniq -u`)
echo "${residual_gadm_list[@]}"

# Create array of OSM country GeoJSONs that haven't been converted to parquet
osm_parquet_finished=()
shopt -s nullglob
for i in $outpath/osm/parquet/*-polygon.parquet; do
    osm_parquet_finished+=("$i")
done
osm_parquet_finished=("${osm_parquet_finished[@]//$outpath\/osm\/parquet\//}")
osm_parquet_finished=("${osm_parquet_finished[@]//-polygon.parquet/}")
echo "${osm_parquet_finished[@]}"

residual_osm_parquet_list=()
residual_osm_parquet_list=(`echo ${iso_list[@]} ${osm_parquet_finished[@]} | tr ' ' '\n' | sort | uniq -u`)
echo "${residual_osm_parquet_list[@]}"

# Download GADM and convert OSM
gadm_list_arg=${residual_gadm_list[@]}
osm_list_arg=${residual_osm_parquet_list[@]}
osm_dir_arg=$outpath/osm/geojson
output_dir_arg=$outpath

python $cwd/subjob_ingestion.py --gadm_chunk $gadm_list_arg --osm_chunk $osm_list_arg --osm_dir $osm_dir_arg --output_dir $output_dir_arg

# # Preprocess / clean up GADM data
# gadm_list_finished=()
# shopt -s nullglob
# for i in $outpath/gadm/geojson/gadm_*.geojson; do
#     gadm_list_finished+=("$i")
# done
# gadm_list_finished=("${gadm_list_finished[@]//$outpath\/gadm\/geojson\//}")
# gadm_list_finished=("${gadm_list_finished[@]//.geojson/}")
# gadm_list_finished=("${gadm_list_finished[@]//gadm_/}")
# 
# residual_gadm_list=()
# residual_gadm_list=(`echo ${iso_list[@]} ${residual_gadm_list[@]} | tr ' ' '\n' | sort | uniq -u`)
# echo "${residual_gadm_list[@]}"
# 
# country_list_arg=${residual_gadm_list[@]}
# gadm_dir_arg=$outpath/gadm/geojson
# daylight_dir_arg=$outpath/daylight/coastlines-v1.19
# osm_parquet_dir_arg=$outpath/osm/parquet
# output_gadm_dir_arg=$outpath/gadm
# 
# python $cwd/subjob_standardize.py --country_chunk $country_list_arg --gadm_dir $gadm_dir_arg --daylight_dir $daylight_dir_arg --osm_dir $osm_parquet_dir_arg --output_dir $output_gadm_dir_arg



