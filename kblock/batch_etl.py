import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
from typing import List, Union
gpd.options.use_pygeos = True

import prepare 
import download

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
np.set_printoptions(suppress=True)

def mem_profile() -> str: 
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use

def main(log_file: Path, country_chunk: list, osm_dir: Path, gadm_dir: Path, output_dir: Path):

    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    #print = logging.getLogger().info
    logging.info(f"Countries to process: {country_chunk}")
    
    # Make directories
    block_dir =  str(output_dir) + '/blocks'
    Path(block_dir).mkdir(parents=True, exist_ok=True)
    street_dir =  str(output_dir) + '/streets'
    Path(street_dir).mkdir(parents=True, exist_ok=True)
    logging.info(f"block_dir: {block_dir}")
    logging.info(f"street_dir: {street_dir}")

    # Check if country is completed
    output_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(block_dir))))))
    output_country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in output_file_list] 
    logging.info(f"Finished countries: {output_country_list}")
    
    country_list = ['DZA','AGO','BEN','BWA','BFA','BDI','CPV','CMR','CAF','TCD','COM','COG','CIV','COD','DJI','EGY','GNQ','ERI','SWZ','ETH','GAB','GMB','GHA','GIN','GNB','KEN','LSO','LBR','LBY','MDG','MWI','MLI','MRT','MUS','MAR','ESH','MOZ','NAM','NER','NGA','RWA','SHN','STP','SEN','SYC','SLE','SOM','ZAF','SSD','SDN','TZA','TGO','TUN','UGA','ZMB','ZWE','AFG','ARM','AZE','BHR','BGD','BTN','BRN','MYS','SGP','KHM','CHN','IND','IDN','IRQ','IRN','PSE','ISR','JPN','JOR','KAZ','KGZ','LAO','LBN','MDV','MNG','MMR','NPL','PRK','PAK','PHL','KOR','LKA','SYR','TWN','TJK','THA','TKM','UZB','VNM','YEM','ALB','AND','AUT','BLR','BEL','BIH','BGR','HRV','CYP','CZE','DNK','EST','FRO','FIN','FRA','GEO','DEU','GRC','HUN','ISL','IRL','IMN','ITA','KOS','LVA','LIE','LTU','LUX','MLT','MDA','MCO','MNE','NLD','MKD','NOR','POL','PRT','ROU','RUS','SRB','SVK','SVN','ESP','SWE','CHE','TUR','UKR','GBR','CAN','GRL','MEX','USA','AUS','COK','FJI','KIR','MHL','FSM','NRU','NCL','NZL','NIU','PLW','PNG','WSM','SLB','TON','TUV','VUT','ARG','BOL','BRA','CHL','COL','ECU','PRY','PER','SUR','URY','VEN']

    # Subset to country list (otherwise will process all countries)
    if country_chunk: 
        country_list = [x for x in country_list if x in set(country_chunk)]
        if not country_list: 
            raise ValueError('Empty country list')
        
    # Subset to countries with no existing output file
    if output_country_list:
        country_list = [x for x in country_chunk if x in set(country_list)]
    logging.info(f"Remaining countries: {country_list}")
    logging.info(f"Generate block geometries")

    # Process country list
    for country_code in country_list: 
        logging.info(f"Processing: {country_code}")
        
        # Download OSM files
        t0 = time.time()
        osm_gpd = download.get_osm_lines(country_code = country_code, download_dir = osm_dir) 
        osm_gpd = osm_gpd.explode(ignore_index = True)
        osm_pygeos = prepare.from_shapely_srid(geometry = osm_gpd, srid = 4326) 
        t1 = time.time()
        logging.info(f"Download OSM: osm_gpd: {osm_gpd.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
    
        # Download GADM files
        t0 = time.time()
        gadm_gpd = download.get_gadm_data(country_code = country_code, download_dir = gadm_dir) 
        gadm_col = max(list(filter(re.compile("GID_*").match, list(gadm_gpd.columns))))
        t1 = time.time()
        logging.info(f"Download GADM: gadm_gpd: {gadm_gpd.shape}, gadm_col: {gadm_col}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        
        # Trim coastline if applicable
        if osm_gpd[osm_gpd['natural'].isin(['coastline'])].shape[0] > 0:
            t0 = time.time()
            gadm_gpd_trim = prepare.trim_coastline(gadm_data = gadm_gpd, osm_data = osm_gpd)
            trimmed_area_percent = (sum(gadm_gpd.to_crs(3395).area/10**6)-sum(gadm_gpd_trim.to_crs(3395).area/10**6))/sum(gadm_gpd_trim.to_crs(3395).area/10**6)
            logging.info(f"Removed area with trim_coastline(): {round(trimmed_area_percent,3)} %")
            gadm_gpd_trim = gpd.clip(gdf = gadm_gpd['geometry'], mask = gadm_gpd_trim.unary_union, keep_geom_type=True)
            gadm_gpd_trim = gpd.GeoDataFrame(geometry=gpd.GeoSeries(gadm_gpd_trim)).reset_index(drop=True)
            gadm_gpd_trim = gpd.overlay(df1 = gadm_gpd_trim, df2 = gadm_gpd, keep_geom_type=True)
            gadm_dropped = [x for x in gadm_gpd[gadm_col].unique() if x not in set(gadm_gpd_trim[gadm_col].unique())]
            gadm_gpd_trim = gadm_gpd_trim.dissolve(by=[gadm_col], as_index = False)
            t1 = time.time()
            logging.info(f"Dropped GADMs from trim_coastline(): {gadm_dropped}")
            logging.info(f"Trimmed gadm_gpd: {gadm_gpd.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
            gadm_gpd_union = gadm_gpd_trim.unary_union
            gadm_gpd_explode = gadm_gpd.explode(ignore_index = True)
            gadm_gpd_residual = gadm_gpd_explode.disjoint(gadm_gpd_union, align=True)
            gadm_gpd = gadm_gpd_trim.append(gadm_gpd_explode[gadm_gpd_residual], ignore_index=True)
            gadm_gpd = gadm_gpd.dissolve(by=[gadm_col], as_index = False)
            del gadm_gpd_trim, gadm_gpd_union, gadm_gpd_explode, gadm_gpd_residual
            logging.info(f"Trimmed with residuals gadm_gpd: {gadm_gpd.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
            
        # Write street geometries
        t0 = time.time()
        osm_streets = prepare.index_streets(osm_data = osm_gpd, gadm_data = gadm_gpd, gadm_column = gadm_col)
        osm_streets.to_parquet(Path(street_dir) / f'streets_{country_code}.parquet', compression='snappy')
        t1 = time.time()
        logging.info(f"Writing osm_streets ({osm_streets.shape}) from osm_gpd ({osm_gpd.shape}) ./streets_{country_code}.parquet, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        del osm_gpd
        
        gadm_list = list(gadm_gpd[gadm_col].unique())
        
        # Initialize
        block_init = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 
        block_boundaries = []
        
        # Build blocks for each GADM
        logging.info(f"Iterate over GADMs")
        t0 = time.time()
        for i in gadm_list: 
            logging.info(f"GADMs: {i}")
            gadm_blocks = prepare.build_blocks(gadm_data = gadm_gpd, osm_data = osm_pygeos, gadm_column = gadm_col, gadm_code = i)
            block_boundaries.append(gadm_blocks)
            check_area = np.sum(gadm_blocks['geometry'].to_crs(3395).area)/np.sum(gadm_gpd[gadm_gpd[gadm_col] == i]['geometry'].to_crs(3395).area)
            if (check_area/1) < .99: logging.info(f"Area proportion: {round(check_area,4)}")
        country_blocks = block_init.append(block_boundaries, ignore_index=True)
        t1 = time.time()
        logging.info(f"build_blocks() function: country_blocks: {country_blocks.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        
        # Write block geometries
        t0 = time.time()
        country_blocks.to_parquet(Path(block_dir) / f'blocks_{country_code}.parquet', compression='snappy')
        t1 = time.time()
        logging.info(f"Writing ./blocks_{country_code}.parquet, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        logging.info(f"Finished {country_code}.")
        
    logging.info(f"Finished.")

def setup(args=None):
    parser = argparse.ArgumentParser(description='Download and build blocks.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format")
    parser.add_argument('--osm_dir', required=True, type=Path, dest="osm_dir", help="OSM directory")
    parser.add_argument('--gadm_dir', required=True, type=Path, dest="gadm_dir", help="GADM directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Output directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))
    
