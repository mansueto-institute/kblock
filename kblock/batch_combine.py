
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

import re
import os
import shutil
from typing import Union
from pathlib import Path
import argparse
import warnings; warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')
import pyogrio

def main(blocks_dir: Path, output_dir: Path):

    # Create list of countries based on block directory
    input_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_file_list]

    #country_list = ['SLE','BEN']

    all_blocks = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326)
 
    for country_code in country_list:
        print(country_code)
        blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
        all_blocks = pd.concat([all_blocks, blocks], ignore_index=True)

    print('Joining...')
    data = pd.read_parquet(path = Path(output_dir) / 'all_blocks_data.parquet')
    data = data.drop(columns=['country_code'])
    all_blocks = all_blocks.merge(data, how = 'left', left_on=['block_id'], right_on = ['block_id'])
    
    print('Writing...')
    #pyogrio.write_dataframe(df = all_blocks, driver = 'GeoJSON', path = Path(output_dir) / f'subsaharan_blocks.geojson')
    all_blocks.to_parquet(path = Path(output_dir) / f'subsaharan_blocks.parquet')
    # all_blocks.to_file(filename = Path(output_dir) / f'subsaharan_blocks.gpkg', driver='GPKG')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute crosswalk.')    
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to outputs directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))