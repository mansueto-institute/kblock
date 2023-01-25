
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True
import pygeos

import re
import os
from pathlib import Path
import argparse
import warnings; warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

def main(blocks_dir: Path, output_dir: Path):

    # Create list of countries based on block directory
    input_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_file_list]

    #country_list = ['SLE','BEN']

    # Blocks
    all_blocks = gpd.GeoDataFrame({
        'block_id': pd.Series(dtype='str'), 
        'block_geohash': pd.Series(dtype='str'), 
        'gadm_code': pd.Series(dtype='str'), 
        'country_code': pd.Series(dtype='str'), 
        'geometry': gpd.GeoSeries(dtype='geometry')}).set_crs(epsg=4326) 

    all_blocks = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': gpd.GeoSeries(dtype='geometry')}).set_crs(epsg=4326)
    for country_code in country_list:
        print(country_code)
        blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
        all_blocks = pd.concat([all_blocks, blocks], ignore_index=True)

    # Population
    all_population = pd.DataFrame({
        'block_id': pd.Series(dtype='str'),
        'landscan_population': pd.Series(dtype= 'float64'), 
        'landscan_population_un': pd.Series(dtype= 'float64'), 
        'worldpop_population': pd.Series(dtype= 'float64'), 
        'worldpop_population_un': pd.Series(dtype= 'float64')
        })

    for country_code in country_list:
        print(country_code) 
        population = pd.read_parquet(path = Path(output_dir) / f'population/population_{country_code}.parquet')
        all_population = pd.concat([all_population, population], ignore_index=True)

    all_population.to_parquet(path = Path(output_dir) / f'combined/all_population.parquet')
    all_blocks = pd.merge(left = all_blocks, right = all_population, how = 'left', on = 'block_id')    
    all_blocks['landscan_population'] = all_blocks['landscan_population'].fillna(0)
    all_blocks['landscan_population_un'] = all_blocks['landscan_population_un'].fillna(0)
    all_blocks['worldpop_population'] = all_blocks['worldpop_population'].fillna(0)
    all_blocks['worldpop_population_un'] = all_blocks['worldpop_population_un'].fillna(0)
    del population, all_population 

    # Buildings
    all_buildings = pd.DataFrame({
        'block_id': pd.Series(dtype='str'),
        'building_area': pd.Series(dtype= 'float64')
        })

    for country_code in country_list:
        print(country_code) 
        buildings = gpd.read_parquet(path = Path(output_dir) / f'buildings/points/buildings_points_{country_code}.parquet')
        buildings = buildings[['block_id', 'building_area']]
        buildings = buildings.groupby('block_id').sum('building_area').reset_index().rename(columns={'building_area': 'building_area'})
        all_buildings = pd.concat([all_buildings, buildings], ignore_index=True)

    all_buildings.to_parquet(path = Path(output_dir) / f'combined/all_buildings.parquet')
    all_blocks = pd.merge(left = all_blocks, right = all_buildings, how = 'left', on = 'block_id')
    del buildings, all_buildings

    # Complexity
    all_complexity = pd.DataFrame({
        'block_id': pd.Series(dtype='str'), 
        'block_area': pd.Series(dtype='float64'), 
        'on_network_street_length': pd.Series(dtype='float64'),  
        'off_network_street_length': pd.Series(dtype='float64'), 
        'nearest_external_street': pd.Series(dtype='float64'), 
        'building_area': pd.Series(dtype='float64'), 
        'building_count': pd.Series(dtype='int'), 
        'building_layers': pd.Series(dtype='object'), 
        'k_complexity': pd.Series(dtype='int')
        })   

    for country_code in country_list:
        print(country_code) 
        complexity = pd.read_parquet(path = Path(output_dir) / f'complexity/complexity_{country_code}.parquet')
        all_complexity = pd.concat([all_complexity, complexity], ignore_index=True)

    all_complexity.to_parquet(path = Path(output_dir) / f'combined/all_complexity.parquet')
    all_blocks = pd.merge(left = all_blocks, right = all_complexity, how = 'left', on = 'block_id')
    all_blocks['block_area'] = all_blocks['block_area'].fillna(0)
    all_blocks['on_network_street_length'] = all_blocks['on_network_street_length'].fillna(0)
    all_blocks['off_network_street_length'] = all_blocks['off_network_street_length'].fillna(0)
    # all_blocks['nearest_external_street'] = all_blocks['nearest_external_street'].fillna(0)
    all_blocks['building_area'] = all_blocks['building_area'].fillna(0)
    all_blocks['building_count'] = all_blocks['building_count'].fillna(0)
    all_blocks['building_layers'] = all_blocks['building_layers'].fillna(['0'])
    all_blocks['k_complexity'] = all_blocks['k_complexity'].fillna(0)
    del complexity, all_complexity

    # Write files
    all_buildings.to_parquet(path = Path(output_dir) / f'combined/all_data_geometry.parquet')
    all_buildings.drop(columns=['geometry']).to_parquet(path = Path(output_dir) / f'combined/all_data.parquet')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute crosswalk.')    
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to outputs directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))

