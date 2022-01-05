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
mkdir -p $outpath/gadm/shp

# Create array of countries that haven't been downloaded from https://www.geofabrik.de/
osm_finished=()
shopt -s nullglob
for i in $outpath/osm/pbf/*-latest.osm.pbf; do
    osm_finished+=("$i")
done
osm_finished=("${osm_finished[@]//$outpath\/osm\/pbf\//}")
osm_finished=("${osm_finished[@]//-latest.osm.pbf/}")
echo "${osm_finished[@]}"

geofabrik_list=(angola burundi benin burkina-faso botswana central-african-republic ivory-coast cameroon congo-democratic-republic congo-brazzaville comores cape-verde djibouti eritrea morocco ethiopia gabon ghana guinea senegal-and-gambia guinea-bissau equatorial-guinea kenya liberia lesotho madagascar mali mozambique mauritania mauritius malawi namibia niger nigeria rwanda sudan senegal-and-gambia sierra-leone somalia south-sudan sao-tome-and-principe swaziland seychelles chad togo tanzania uganda south-africa zambia zimbabwe)

residual_osm_list=()
residual_osm_list=(`echo ${geofabrik_list[@]} ${osm_finished[@]} | tr ' ' '\n' | sort | uniq -u`)
echo "${residual_osm_list[@]}"

for country in ${residual_osm_list[@]}; do
echo "$country"
    wget -O $outpath/osm/pbf/${country}-latest.osm.pbf https://download.geofabrik.de/africa/${country}-latest.osm.pbf
    osmium export $outpath/osm/pbf/${country}-latest.osm.pbf -o $outpath/osm/${country}-latest-linestring.geojson -c $cwd/export-config.json --overwrite --geometry-types=linestring
done

# Create array of countries that haven't been downloaded from https://gadm.org/
gadm_finished=()
shopt -s nullglob
for i in $outpath/gadm/shp/*; do
    gadm_finished+=("$i")
done
gadm_finished=("${gadm_finished[@]//$outpath\/gadm\/shp\//}")
echo "${gadm_finished[@]}"

gadm_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)

residual_gadm_list=()
residual_gadm_list=(`echo ${gadm_list[@]} ${gadm_finished[@]} | tr ' ' '\n' | sort | uniq -u`)
echo "${residual_gadm_list[@]}"

gadm_list_arg=${residual_gadm_list[@]}
output_dir_arg=$outpath/gadm

python $cwd/download_gadm.py --country_chunk $gadm_list_arg --output_dir $output_dir_arg
