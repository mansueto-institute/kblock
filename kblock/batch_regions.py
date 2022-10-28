
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
import dask_geopandas

def main(ghsl_dir: Path, blocks_dir: Path, output_dir: Path):

    # Setup data
    ghsl_data = gpd.read_file(Path(ghsl_dir))
    ghsl_data = ghsl_data[ghsl_data['GRGN_L2'].isin(['Western Africa', 'Northern Africa', 'Middle Africa', 'Southern Africa', 'Eastern Africa'])][['ID_HDC_G0','GCPNT_LAT', 'GCPNT_LON', 'XBRDR', 'XCTR_NBR', 'CTR_MN_ISO', 'CTR_MN_NM', 'UC_NM_MN','UC_NM_LST','XC_NM_LST','XC_ISO_LST','GRGN_L1','GRGN_L2','geometry']].reset_index()
    ghsl_data = ghsl_data.to_crs(3395)
    ghsl_data['urban_area'] = ghsl_data['geometry'].area
    ghsl_data['urban_radius'] = pygeos.minimum_bounding_radius(pygeos.from_shapely(ghsl_data['geometry']))
    col_recode = {'ID_HDC_G0':'urban_id', 'GCPNT_LAT': 'latitude' , 'GCPNT_LON': 'longitude', 'XBRDR':'crossborder', 'XCTR_NBR':'num_country_crossborder', 'CTR_MN_ISO':'primary_country_code', 'CTR_MN_NM':'primary_country_name', 'UC_NM_MN':'urban_center_name', 'UC_NM_LST':'list_of_urban_center_names', 'XC_NM_LST':'intersected_countries', 'XC_ISO_LST':'intersected_codes', 'GRGN_L1':'geographical_region', 'GRGN_L2':'geographical_subregion'}
    ghsl_data = ghsl_data.rename(columns=col_recode)
    ghsl_data = ghsl_data[['urban_id','urban_center_name','list_of_urban_center_names','primary_country_code','primary_country_name','geographical_subregion','geographical_region','urban_area','urban_radius','latitude', 'longitude', 'geometry']]
    ghsl_data = ghsl_data[ghsl_data['urban_center_name'] != 'N/A']
    ghsl_data = ghsl_data.explode(index_parts=False)
    assert ghsl_data[ghsl_data['urban_id'].duplicated()].shape[0] == 0
    
    # Urban centers
    urban_centers = ghsl_data[['urban_id','geometry']].to_crs(4326)
    
    # Periurban zones based on minimum bounding circle
    periurban_inner = pygeos.from_shapely(ghsl_data['geometry'].buffer(ghsl_data['urban_radius']*1))
    periurban_inner = pygeos.get_parts(pygeos.union_all(periurban_inner))
    periurban_inner = gpd.GeoDataFrame.from_dict({'geometry': pygeos.to_shapely(periurban_inner)}).set_crs(epsg=3395)
    periurban_inner = periurban_inner.assign(conurbation_id = [str(x) for x in list(periurban_inner.index)])
    periurban_inner = periurban_inner[['conurbation_id','geometry']].to_crs(4326)
    
    periurban_residual = periurban_inner.overlay(urban_centers, how='difference')
    urban_centers = pd.concat([urban_centers,periurban_residual])
    
    # Create list of countries based on block directory
    input_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_file_list]

    # country_list = ['ZMB']
 
    for country_code in country_list:
        print(country_code)
    
        blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
    
        # Urban centers 
        urban_overlay = gpd.overlay(df1 = blocks, df2 = urban_centers, how='intersection')
        urban_overlay['area'] = urban_overlay['geometry'].to_crs(3395).area
        urban_overlay['rank'] = urban_overlay.groupby('block_id')['area'].rank(method='first', ascending=False)
        urban_overlay = urban_overlay[(urban_overlay['urban_id'].notnull()) & (urban_overlay['rank'] == 1)]
        urban_overlay = urban_overlay[['block_id','urban_id','rank']]
    
        blocks_xwalk = blocks.merge(urban_overlay, how = 'left', left_on=['block_id'], right_on = ['block_id'])
        blocks_xwalk = blocks_xwalk[['block_id','gadm_code','country_code','geometry','urban_id']]
    
        # Conurbations
        conurbation_overlay = gpd.overlay(df1 = blocks, df2 = periurban_inner, how='intersection')
        conurbation_overlay['area'] = conurbation_overlay['geometry'].to_crs(3395).area
        conurbation_overlay['rank'] = conurbation_overlay.groupby('block_id')['area'].rank(method='first', ascending=False)
        conurbation_overlay = conurbation_overlay[(conurbation_overlay['conurbation_id'].notnull()) & (conurbation_overlay['rank'] == 1)]
        conurbation_overlay = conurbation_overlay[['block_id','conurbation_id','rank']]
    
        blocks_xwalk = blocks_xwalk.merge(conurbation_overlay, how = 'left', left_on=['block_id'], right_on = ['block_id'])
        blocks_xwalk = blocks_xwalk[['block_id','gadm_code','country_code','geometry','urban_id','conurbation_id']]
    
        if blocks_xwalk[blocks_xwalk['block_id'].duplicated()].shape[0] > 0:
            print('Duplicated rows')
    
        # Blocks crosswalk
        blocks_xwalk = blocks_xwalk[['block_id','gadm_code','country_code','urban_id','conurbation_id']]
        blocks_xwalk['conurbation_urban_id'] = blocks_xwalk['country_code'].astype(str) + '_' + blocks_xwalk['conurbation_id'].astype(str) + '_' + blocks_xwalk['urban_id'].astype(str)

        # Read blocks
        blocks_all = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
        blocks_all = blocks_all.merge(right = blocks_xwalk[['block_id','conurbation_urban_id']], how='inner', on='block_id')
    
        region = gpd.GeoDataFrame({'conurbation_urban_id': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326)

        for i in blocks_all['conurbation_urban_id'].unique():
            print(i)
            if i == str(country_code + '_' + str(np.nan) + '_' + str(np.nan)):
                resid_code = str(country_code + '_' + str(np.nan) + '_' + str(np.nan))
                residual_blocks = blocks_all[blocks_all['conurbation_urban_id'] == resid_code][['conurbation_urban_id','gadm_code','geometry']]
                residual_blocks = dask_geopandas.from_geopandas(residual_blocks, npartitions = 50)
                residual_blocks =  residual_blocks.dissolve(by = "gadm_code").compute()
                residual_blocks.reset_index(inplace = True)
                residual_blocks['group_col'] = residual_blocks['gadm_code'].str.split(".").str[:3].str.join(".")
                residual_blocks = dask_geopandas.from_geopandas(residual_blocks, npartitions = 10)
                residual_blocks =  residual_blocks.dissolve(by = "group_col").compute()
                residual_blocks.reset_index(inplace = True)
                residual_blocks = dask_geopandas.from_geopandas(residual_blocks, npartitions = 5)
                residual_blocks = residual_blocks.dissolve(by = "conurbation_urban_id").compute()
                residual_blocks.reset_index(inplace = True)
                blocks_sub = residual_blocks[['conurbation_urban_id','geometry']]
            else: 
                #blocks_sub = pygeos.from_shapely(blocks_all[blocks_all['conurbation_urban_id'] == i]['geometry'])
                #blocks_sub = gpd.GeoSeries(pygeos.to_shapely(pygeos.union_all(blocks_sub)))
                #blocks_sub = gpd.GeoDataFrame.from_dict({'conurbation_urban_id': i, 'geometry': blocks_sub}).set_crs(epsg=4326)
                blocks_sub = blocks_all[blocks_all['conurbation_urban_id'] == i][['conurbation_urban_id','geometry']]
                blocks_sub = dask_geopandas.from_geopandas(blocks_sub, npartitions = 50)
                blocks_sub =  blocks_sub.dissolve(by = "conurbation_urban_id").compute()
                blocks_sub.reset_index(inplace = True)
                blocks_sub = blocks_sub[['conurbation_urban_id','geometry']]
            region = pd.concat([region, blocks_sub], ignore_index=True)
        
        print('Writing...')
        region.to_parquet(path = Path(output_dir) / f'regions_{country_code}.parquet')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute crosswalk.')
    parser.add_argument('--ghsl_dir', required=True, type=Path, dest="ghsl_dir", help="Path to GHSL directory")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to outputs directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))

