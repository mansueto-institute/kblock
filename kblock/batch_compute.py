import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
from typing import List, Union
gpd.options.use_pygeos = True

import prepare 
import compute

from pathlib import Path
from urlpath import URL
import psutil
import contextlib
import warnings
import logging
import sys
import time
import re
import os
import io
import argparse

import pyarrow
import dask 
import dask.bag as db
import dask_geopandas
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

#from dask_mpi import initialize
#initialize() #local_directory=os.path.dirname(Path(log_file))

#from dask.distributed import Client, LocalCluster
#dask.config.set(scheduler='processes')
#partition_count = 10

# Set 30 MB per partition (Dask recommends 100 MB)
#mem_estimate = osm_gpd.memory_usage(index=True, deep=True).sum()/1000000 + gadm_gpd.memory_usage(index=True, deep=True).sum()/1000000
#partition_count = round(mem_estimate/30)

def mem_profile() -> str: 
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use

country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}    

def bulk_compute_k(gadm_code: str, country_code: Union[str, Path], blocks_dir: Union[str, Path], streets_dir: Union[str, Path], buildings_dir: Union[str, Path], building_file_list: list) -> gpd.GeoDataFrame:

    print(gadm_code)
    
    gadm_blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True, filters = [('gadm_code', '=', gadm_code)]).to_crs(3395)
    
    try:
        street_net = gpd.read_parquet(path = Path(streets_dir) / f'streets_{country_code}.parquet', memory_map = True, filters = [('gadm_code', '=', gadm_code)]).to_crs(3395)
        street_net = prepare.from_shapely_srid(geometry = street_net, srid = 3395) 
    except Warning:
        print(f"Warning: No street data available for: {gadm_code}")

    buildings_file = list(filter(re.compile(str("%s" % gadm_code +'.geojson')).findall, sorted(building_file_list)))[0]
    buildings = gpd.read_file(Path(buildings_dir) / country_code / buildings_file).to_crs(3395) 
    buildings = compute.index_buildings(gadm_block_data = gadm_blocks, bldg_data = buildings)
    buildings = buildings.to_crs(3395)

    block_list = list(gadm_blocks['block_id'].unique())
    k_gadm = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'block_area': pd.Series(dtype='float'), 'building_area': pd.Series(dtype='float'), 'building_count': pd.Series(dtype='int'), 'building_layers': pd.Series(dtype='object'),  'k_complexity': pd.Series(dtype='int')})    
        
    #for x in block_list: 
    #    k_block = kblock.compute_k(block_data = gadm_blocks, bldg_data = buildings, block_col = 'block_id', block_id = x, street_linestrings = street_net, buffer_streets = True, include_geometry = False)
    #    k_gadm = pd.concat([k_gadm, k_block], ignore_index=True)
    k_block = list(map(lambda x: compute.compute_k(block_data = gadm_blocks, bldg_data = buildings, block_col = 'block_id', block_id = x, street_linestrings = street_net, buffer_streets = True, include_geometry = False), block_list))
    for i,j in enumerate(k_block): k_gadm = pd.concat([k_gadm, k_block[i]], ignore_index=True)
    
    return k_gadm


def main(log_file: Path, country_chunk: list, blocks_dir: Path, streets_dir: Path, buildings_dir: Path, output_dir: Path):

    #initialize(local_directory=os.path.dirname(Path(log_file))) 
    #client = Client() 

    # List of files in the given directory 
    #list_of_files = filter(os.path.isfile,glob.glob('/Users/nm/Downloads/production/inputs/buildings/SLE/' + '*') )
    # Sort list of files in directory by size 
    #list_of_files = sorted(list_of_files,key =  lambda x: os.stat(x).st_size)

    # Make directory
    metrics_dir =  str(output_dir) + '/metrics'
    Path(metrics_dir).mkdir(parents=True, exist_ok=True)
    logging.info(f"metrics_dir: {metrics_dir}")

    # Check if country is completed
    output_file_list = list(filter(re.compile("metrics_").match, sorted(list(os.listdir(Path(metrics_dir))))))
    output_country_list = [(re.sub('metrics_', '', re.sub('.parquet', '', i))) for i in output_file_list] 
    logging.info(f"Completed countries in output directory: {output_country_list}")

    # Subset to country list (remove invalid country codes and remove completed countries)
    if country_chunk: 
        country_list = list(country_dict.keys())
        country_list = [x for x in country_chunk if x in set(country_list)]
        country_list = [x for x in country_chunk if x not in set(output_country_list)]
        if not country_list: 
            raise ValueError('Empty country list')
        
    logging.info(f"Remaining countries: {country_list}")

    # Process country list
    for country_code in country_list:
        
        logging.info(f"Processing: {country_code}")
        t0_country = time.time()

        # Create list of GADMs based on block file
        country_blocks = dask_geopandas.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)    
        gadm_list = list(country_blocks['gadm_code'].unique())
        del country_blocks
        logging.info(f"GADMs to process: {gadm_list}")

        # Validate GADM codes from country file against building directory codes
        building_file_list = list(filter(re.compile("buildings_").match, sorted(list(os.listdir(Path(buildings_dir) / country_code)))))  
        building_file_gadm_list = [(re.sub('buildings_', '', re.sub('.geojson', '', i))) for i in building_file_list] 
        gadm_list_match = [x for x in gadm_list if x in set(building_file_gadm_list)]
        gadm_in_block_not_in_building = [x for x in gadm_list if x not in set(building_file_gadm_list)]
        gadm_in_building_not_in_block = [x for x in building_file_gadm_list if x not in set(gadm_list)]
        logging.info(f"GADMs in block file and buildings directory: {gadm_list_match}")
        logging.info(f"GADMs in block file not in buildings directory: {gadm_in_block_not_in_building}")
        logging.info(f"GADMs in buildings directory not in block file: {gadm_in_building_not_in_block}")

        # 
        gadm_list = gadm_list_match

        # Initialize
        k_bulk = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'block_area': pd.Series(dtype='float'), 'building_area': pd.Series(dtype='float'), 'building_count': pd.Series(dtype='int'), 'building_layers': pd.Series(dtype='object'),  'k_complexity': pd.Series(dtype='int')})    

        # Process each GADM / parallelize
        logging.info(f"Iterate over GADMs")
        #bag_sequence = db.from_sequence(gadm_list, npartitions=partition_count) 
        ##compute_sequence = bag_sequence.map_partitions(lambda x: bulk_compute_k(gadm_code = x, country_code = country_code, blocks_dir = blocks_dir, streets_dir = streets_dir, buildings_dir = buildings_dir, building_file_list = building_file_list))
        #compute_sequence = bag_sequence.map(lambda x: bulk_compute_k(gadm_code = x, country_code = country_code, blocks_dir = blocks_dir, streets_dir = streets_dir, buildings_dir = buildings_dir, building_file_list = building_file_list))
        #output_sequence = compute_sequence.compute() # persist()
        #for i,j in enumerate(output_sequence): k_bulk = pd.concat([k_bulk, output_sequence[i]], ignore_index=True)

        output_map = list(map(lambda x: bulk_compute_k(gadm_code = x, country_code = country_code, blocks_dir = blocks_dir, streets_dir = streets_dir, buildings_dir = buildings_dir, building_file_list = building_file_list), gadm_list))
        for i,j in enumerate(output_map): k_bulk = pd.concat([k_bulk, output_map[i]], ignore_index=True)
        
        # Write parquet file (build file incrementally and check how much is done)
        k_bulk.to_parquet(path = Path(metrics_dir) / f'metrics_{country_code}.parquet') 
        t1_country = time.time()
        logging.info(f"Finished processing metrics: {country_code}, {str(round((t1_country-t0_country)/60,3))} minutes")


def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute K.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Blocks directory")
    parser.add_argument('--streets_dir', required=True, type=Path, dest="streets_dir", help="Streets directory")
    parser.add_argument('--buildings_dir', required=True, type=Path, dest="buildings_dir", help="Buildings directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Outputs directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))

# capture logs

# use map instead of loops
# option to skip over checks 
# figure out speed ups with list comprehensions *
# replace appends with newlist = map(str.upper, wordlist) *

# incremental file writes with parquet or CSV (write to "optimal GADM level") 
#    after workflow ends append all the CSVs together, write to parquet, then delete the CSVs
# check for completed CSVs (not just completed countries and pick up after they complete)

# dask_geopandas where it calls geopandas
# use dask bag with mpirun
