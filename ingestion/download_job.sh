#!/bin/bash

mkdir -p osm/pbf
mkdir -p gadm/shp

# Create array of countries that haven't been downloaded from https://www.geofabrik.de/
osm_finished=()
for i in osm/pbf/*-latest.osm.pbf; do
    osm_finished+=("$i")
done
osm_finished=("${osm_finished[@]//osm\/pbf\//}")
osm_finished=("${osm_finished[@]//-latest.osm.pbf/}")
echo "${osm_finished[@]}"

geofabrik_list=(angola burundi benin burkina-faso botswana central-african-republic ivory-coast cameroon congo-democratic-republic congo-brazzaville comores cape-verde djibouti eritrea morocco ethiopia gabon ghana guinea senegal-and-gambia guinea-bissau equatorial-guinea kenya liberia lesotho madagascar mali mozambique mauritania mauritius malawi namibia niger nigeria rwanda sudan senegal-and-gambia sierra-leone somalia south-sudan sao-tome-and-principe swaziland seychelles chad togo tanzania uganda south-africa zambia zimbabwe)

residual_osm_list=()
residual_osm_list=(`echo ${osm_finished[@]} ${geofabrik_list[@]} | tr ' ' '\n' | sort | uniq -u`)
echo "${residual_osm_list[@]}"

for country in ${residual_osm_list[@]}; do
echo "$country"
    wget -O osm/pbf/${country}-latest.osm.pbf https://download.geofabrik.de/africa/${country}-latest.osm.pbf
    osmium export osm/pbf/${country}-latest.osm.pbf -o osm/${country}-latest-linestring.geojson -c /Users/nm/Desktop/download/export-config.json --overwrite --geometry-types=linestring
done

# Create array of countries that haven't been downloaded from https://gadm.org/
gadm_finished=()
for i in gadm/shp/*; do
    gadm_finished+=("$i")
done
gadm_finished=("${gadm_finished[@]//gadm\/shp\//}")
echo "${gadm_finished[@]}"

gadm_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)

residual_gadm_list=()
residual_gadm_list=(`echo ${gadm_finished[@]} ${gadm_list[@]} | tr ' ' '\n' | sort | uniq -u`)
echo "${residual_gadm_list[@]}"

gadm_list_arg=${residual_gadm_list[@]}
output_dir_arg=gadm

python /Users/nm/Desktop/download/download_gadm.py --country_chunk $gadm_list_arg --output_dir $output_dir_arg

