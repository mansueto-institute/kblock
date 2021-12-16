
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
gpd.options.use_pygeos = True

import kblock

import pyarrow
import time
import re
import sys
import os
import io
import argparse
from typing import List, Union
from pathlib import Path
from contextlib import redirect_stderr, redirect_stdout
import logging
import warnings
import block_summary
import multiprocessing

def gadm_dir_to_path(gadm_dir: Union[str, Path]) -> str:
    """
    For a given country, the GADM dir contains multiple file levels
    so this convenience function just returns the path to the highest
    resolution gadm file within the directory, which changes from 
    country to country
    Args:
        gadm_dir: directory containing all gadm files for a country
    Returns:
        Path to specific gadm file
    """
    sort_fn = lambda p: int(p.stem.split("_")[-1])
    gadm_dir = Path(gadm_dir)
    files = [p for p in gadm_dir.iterdir() if ".shp" in p.name]
    files.sort(key=sort_fn)
    return files[-1]

def main(log_file: Path, country_code: str, country_code_file: Path, gadm_parent_dir: Path, 
         gadm_chunk: list, streets_parent_dir: Path, building_parent_dir: Path, population_raster_path: Path,
         output_dir: Path):

    logging.getLogger().setLevel("DEBUG")
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f"Starting: {country_code}")

    # Output Directory
    output_dir_country = Path(output_dir) / country_code
    output_dir_country.mkdir(parents=True, exist_ok=True)
    logging.info(f"output_dir_country: {output_dir_country}")
    output_file_list = list(filter(re.compile("kblock_").match, sorted(list(os.listdir(output_dir_country)))))
    output_gadm_list = [(re.sub('kblock_', '', re.sub('.geojson', '', i))) for i in output_file_list] 

    # GADMs
    gadm_dir = Path(gadm_parent_dir) / country_code
    gadm_file = gadm_dir_to_path(gadm_dir = gadm_dir)

    t0 = time.time()
    logging.info(f"gadm_file: {gadm_file}")
    gadm_gpd = gpd.read_file(gadm_file).to_crs(4326)
    logging.info(f"gadm_gpd.shape: {gadm_gpd.shape}")
    gadm_col = max(list(filter(re.compile("GID_*").match, list(gadm_gpd.columns))))
    logging.info(f"gadm_col: {gadm_col}")
    gadm_list = list(gadm_gpd[gadm_col].unique())
    t1 = time.time()
    logging.info(f"Read GADM time: {round(t1-t0,5)}")

    # Subset to GADM list
    if gadm_chunk: 
        gadm_list = [x for x in gadm_list if x in set(gadm_chunk)]
        logging.info(f"Subset GADM list: {gadm_list}")
        if not gadm_list: 
            raise ValueError('Empty GADM list')
    logging.info(f"Processed GADM chunk")

    if output_gadm_list:
        gadm_list = [x for x in gadm_list if x not in set(output_gadm_list)]
        logging.info(f"Remove completed GADMs: {gadm_list}")
    logging.info(f"Processed GADM list")

    # OSM Directory
    country_metadata = pd.read_csv(country_code_file)
    logging.info(f"Processed country_code_file")
    geofabrik_name = list(country_metadata[country_metadata['country_code'] == country_code]['geofabrik_name'])
    osm_file = Path(streets_parent_dir) / str(geofabrik_name[0]+'_lines.geojson')

    # OSM linestrings
    logging.info(f"osm_file: {osm_file}")
    t0 = time.time()
    osm_gpd = gpd.read_file(osm_file).to_crs(4326)
    osm_pygeos = kblock.from_shapely_srid(geometry = osm_gpd, srid = 4326) 
    t1 = time.time()
    logging.info(f"Read OSM time: {round(t1-t0,5)}")
    logging.info(f"osm_gpd.shape: {osm_gpd.shape}")

    # Trim coastline
    # if osm_gpd[osm_gpd['natural'].isin(['coastline','water'])].shape[0] > 0:
    #     t0 = time.time()
    #     gadm_gpd_trim = kblock.trim_coastline(gadm_data = gadm_gpd, osm_data = osm_gpd)
    #     trimmed_area_percent = (sum(gadm_gpd.to_crs(3395).area/10**6)-sum(gadm_gpd_trim.to_crs(3395).area/10**6))/sum(gadm_gpd_trim.to_crs(3395).area/10**6)
    #     gadm_gpd = gpd.clip(gdf = gadm_gpd, mask = gadm_gpd_trim, keep_geom_type=True)
    #     t1 = time.time()
    #     logging.info(f"Trim coastline time: {round(t1-t0,5)}")
    #     logging.info(f"Trimmed percentage: {trimmed_area_percent}")
    #     logging.info(f"gadm_gpd_trim.shape: {gadm_gpd.shape}")
    #     with open(Path(os.path.dirname(Path(log_file))) / '_log_trim_coast.txt', 'a') as f: 
    #         with redirect_stdout(f):
    #             print(f'{country_code}, {round(trimmed_area_percent,5)}')

    # Building directory
    logging.info(f"Check building directory")
    building_file_list = list(filter(re.compile("buildings_").match, sorted(list(os.listdir(Path(building_parent_dir) / country_code)))))  
    building_file_gadm_list = [(re.sub('buildings_', '', re.sub('.geojson', '', i))) for i in building_file_list] 
    gadm_list_match = [x for x in gadm_list if x in set(building_file_gadm_list)]
    gadm_list_no_match = [x for x in gadm_list if x not in set(building_file_gadm_list)]
    building_file_list_no_match = [x for x in building_file_gadm_list if x not in set(gadm_list)]
    logging.info(f"Units in GADM file not found in building directory: {gadm_list_no_match}")
    logging.info(f"Units in building directory not found in GADM file: {building_file_list_no_match}")

    # Iterate through GADMs:
    gadm_data_list = [{'gadm': i, 'gadm_gpd': gadm_gpd, 'osm_pygeos': osm_pygeos, 'gadm_column': gadm_col,
                       'log_file': log_file, 'building_parent_dir': building_parent_dir,
                       'country_code': country_code, 'building_file_list': building_file_list,
                       'output_dir_country': output_dir_country, 'osm_gpd': osm_gpd,
                       'population_raster_path': population_raster_path} for i in gadm_list]
    logging.info(f'Running over GADMs: {gadm_list}')

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        pool.map(main_helper, gadm_data_list)
    logging.info('Finished')


def main_helper(gdd: dict) -> None:
    # Initialize GeoDataFrame
    k_init = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 
        'block_area': pd.Series(dtype='float'), 'building_area': pd.Series(dtype='float'), 
        'building_count': pd.Series(dtype='int'), 'building_layers': pd.Series(dtype='object'),  'k_complexity': pd.Series(dtype='int'), 
        'block_pop_ls': pd.Series(dtype='float'), 'block_pop_wp': pd.Series(dtype='float'),
        'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 
    logging.info(f'k_init: {k_init}')

    logging.getLogger().setLevel("DEBUG")
    logging.basicConfig(filename=Path(gdd['log_file']), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f"Starting: {gdd['country_code']}")
    logging.info(f"GADM: {gdd['gadm']}")

    t0 = time.time()
    gadm_blocks = kblock.build_blocks(gadm_data = gdd['gadm_gpd'], osm_data = gdd['osm_pygeos'], 
                                      gadm_column = gdd['gadm_column'], gadm_code = gdd['gadm'])
    gadm_blocks_list = list(gadm_blocks['block_id'].unique())
    logging.info(f'Build blocks success')
    t1 = time.time()
    check_area = np.sum(gadm_blocks['geometry'].to_crs(3395).area)/np.sum(gdd['gadm_gpd'][gdd['gadm_gpd'][gdd['gadm_column']] == gdd['gadm']]['geometry'].to_crs(3395).area)
    logging.info(f"Subdivide GADM into blocks: {round(t1-t0,5)}")
    logging.info(f"gadm_blocks area / gadm_gpd area: {check_area}")
    logging.info(f"gadm_blocks.shape: {gadm_blocks.shape}")
    with open(Path(os.path.dirname(Path(gdd['log_file']))) / '_log_block_integrity.txt', 'a') as f: 
        with redirect_stdout(f):
            print(f'{gdd["gadm"]}, {round(check_area,5)}')

    t0 = time.time()
    building_file = list(filter(re.compile(str("%s" % gdd['gadm'] +'.geojson')).findall, sorted(gdd['building_file_list'])))[0]
    building_gpd = gpd.read_file(Path(gdd['building_parent_dir']) / gdd['country_code'] / building_file).to_crs(4326) 
    t1 = time.time()
    logging.info(f"Building file read: {round(t1-t0,5)}")
    logging.info(f"building_gpd.shape: {building_gpd.shape}")
    
    t0 = time.time()
    block_coded_buildings = kblock.index_buildings(gadm_block_data = gadm_blocks, bldg_data = building_gpd)
    t1 = time.time()
    logging.info(f"Index building time: {round(t1-t0,5)}")
    logging.info(f"block_coded_buildings.shape: {block_coded_buildings.shape}")
    
    block_coded_buildings = block_coded_buildings.to_crs(3395)
    gadm_blocks = gadm_blocks.to_crs(3395)
    #osm_highways = pygeos.from_shapely(osm_gpd[osm_gpd['highway'].notnull()]['geometry'].to_crs(epsg=3395))
    osm_highways = kblock.from_shapely_srid(geometry = gdd['osm_gpd'][gdd['osm_gpd']['highway'].notnull()].to_crs(3395), srid = 3395) 
    logging.info(f"osm_highways.shape: {gdd['osm_gpd'][gdd['osm_gpd']['highway'].notnull()].shape}")

    # Iterate through blocks:
    logging.info("Compute k")
    block_metrics = []
    for x in gadm_blocks_list: 
        t0 = time.time()
        with open(Path(os.path.dirname(Path(gdd['log_file']))) / '_log_street_buffer.txt', 'a') as f: 
            with redirect_stdout(f):
                df_k = kblock.compute_k(block_data = gadm_blocks,
                                        bldg_data = block_coded_buildings, 
                                        block_col = 'block_id', 
                                        block_id = x,
                                        street_linestrings = osm_highways,
                                        buffer_streets=True)
        block_metrics.append(df_k)
        t1 = time.time()
        logging.info(f"block_id: {x} - {round(t1-t0,5)}")

    t0 = time.time()
    kblock_w_pop = block_summary.make_summary(block_metrics, gdd['population_raster_path'], building_gpd, gdd['log_file'])
    t1 = time.time()
    kblock_output_w_pop = k_init.append(kblock_w_pop, ignore_index=True)
    logging.info(f"Block statistics time: {round(t1-t0,5)}")
    kblock_w_pop.to_file(Path(output_dir_country) / str('kblock_'+gdd['gadm']+'.geojson'), driver='GeoJSON')



def setup(args=None):
    parser = argparse.ArgumentParser(description='kblock computation')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="path to write log file") 
    parser.add_argument('--country_code', required=True, type=str, dest="country_code", help="3 digit country code")
    parser.add_argument('--country_code_file', required=True, type=Path, dest="country_code_file", help="csv with 3 digit country codes")
    parser.add_argument('--gadm_chunk', required=False, type=str, dest="gadm_chunk", nargs='+', help="subset list of GADM codes")
    parser.add_argument('--gadm_parent_dir', required=True, type=Path, dest="gadm_parent_dir", help="GADM parent directory")
    parser.add_argument('--streets_parent_dir', required=True, type=Path, dest="streets_parent_dir", help="list of GADM codes")
    parser.add_argument('--building_parent_dir', required=True, type=Path, dest="building_parent_dir", help="building file parent directory")
    parser.add_argument('--population_raster_path', required=True, type=Path, dest="population_raster_path", help="Filepath for population raster data")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="output directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))
