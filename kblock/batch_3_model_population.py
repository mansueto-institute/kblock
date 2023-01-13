
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

from pathlib import Path
from typing import List, Union
from pathlib import Path

import time
import re
import psutil
import os
import argparse
import warnings
import contextlib
import logging

import rasterio as rio
from rasterio import features
import xarray as xr
import rioxarray as rxr

np.set_printoptions(suppress=True)
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')


def rasterize_geometry(gpd_dataframe: gpd.GeoDataFrame, gpd_col: str, xr_array: xr.DataArray, fill=np.nan, encode=False, **kwargs) -> xr.DataArray:
    """Burn a set of geometries into the given raster dataset
    Args:
        gpd_dataframe: geopandas.GeoDataFrame containing geometries to burn into raster
        xr_array: xarray.DataArray of target raster containing dimensions 'y' and 'x'
        encode: boolean, defaults to False, encodes geometry column label as integers
        fill: (optional) int or float, used as fill value for all areas not covered by input geometries
    Returns:
        xarray.DataArray in shape of 'x' and 'y' raster dimensions
    """
    bounds = xr_array.rio.bounds()
    transform_constructor = rio.transform.from_bounds(west = bounds[0], south = bounds[1], east = bounds[2], north = bounds[3], width = len(xr_array['x']), height = len(xr_array['y']))
    dims_remove = [x for x in list(xr_array.dims) if x not in set(['x','y'])]
    xr_dataset = xr_array.to_dataset(name='data').drop_dims(drop_dims = dims_remove)
    shapes = [(shape, n) for n, shape in enumerate(gpd_dataframe.geometry)]
    shapes_map = {n : shape for n, shape in enumerate(gpd_dataframe[gpd_col])}    
    out_shape = (len(xr_dataset['y']), len(xr_dataset['x']))
    raster = features.rasterize(shapes, out_shape=out_shape, fill=fill, transform=transform_constructor, all_touched = True, dtype=float, **kwargs)
    raster_map = np.array(list(map(np.vectorize(shapes_map.get), raster)))
    if encode is True:
        gpd_col_encode = str(gpd_col + '_int')
        xr_dataset[gpd_col_encode] = xr.DataArray(data=raster, coords=xr_dataset, dims=('y', 'x'))
    else: 
        xr_dataset[gpd_col] = xr.DataArray(data=raster_map, coords=xr_dataset, dims=('y', 'x'))
    return xr_dataset


def vectorize_pixels(xr_array: xr.DataArray) ->  xr.Dataset:
    """Vectorize pixel coordinates in a xarray.DataArray raster
    Args:
        xr_array: xarray.DataArray of target raster containing dimensions 'y' and 'x'
    Returns:
        xarray.Dataset with corner coordinates of each pixel as data variables
    """
    bounds = xr_array.rio.bounds()
    transform_constructor = rio.transform.from_bounds(west = bounds[0], south = bounds[1], east = bounds[2], north = bounds[3], width = len(xr_array['x']), height = len(xr_array['y']))
    height = len(xr_array['y'])
    width = len(xr_array['x'])
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    
    xs, ys = rio.transform.xy(transform = transform_constructor, rows = rows, cols = cols, offset = 'center')
    xs_ul, ys_ul = rio.transform.xy(transform = transform_constructor, rows = rows, cols = cols, offset = 'ul')
    xs_ur, ys_ur = rio.transform.xy(transform = transform_constructor, rows = rows, cols = cols, offset = 'ur')
    xs_ll, ys_ll = rio.transform.xy(transform = transform_constructor, rows = rows, cols = cols, offset = 'll')
    xs_lr, ys_lr = rio.transform.xy(transform = transform_constructor, rows = rows, cols = cols, offset = 'lr')
    
    lons = np.array(xs)
    lats = np.array(ys)
    lons_ul = np.array(xs_ul)
    lats_ul = np.array(ys_ul)
    lats_ll = np.array(ys_ll)
    lons_ur = np.array(xs_ur)
    
    dims_remove = [x for x in list(xr_array.dims) if x not in set(['x','y'])]
    xr_data = xr_array.to_dataset(name='data').drop_dims(drop_dims = dims_remove)
    xr_data['ymax'] = xr.DataArray(data=lats_ul, coords = xr_data.coords, dims=('y', 'x'))
    xr_data['ymin'] = xr.DataArray(data=lats_ll, coords = xr_data.coords, dims=('y', 'x'))
    xr_data['xmax'] = xr.DataArray(data=lons_ur, coords = xr_data.coords, dims=('y', 'x'))
    xr_data['xmin'] = xr.DataArray(data=lons_ul, coords = xr_data.coords, dims=('y', 'x'))
    
    return xr_data

def read_tif_dir(dir_path: Union[str, Path], country_code: str, gadm_gpd: gpd.GeoDataFrame, gadm_code_list = []) -> xr.DataArray:
    """ Read GeoTiff or .tif file format from LandScan and Worldpop for a particular country
    Args:
        dir_path: string or Path to the file to open, director must have this structure:
            Path(dir_path) / 'landscan/<filename>.tif' (single global file)
            Path(dir_path) / f'worldpop/{country_code}.tif' (multiple country files)
        country_code: string of 3-digit country codes ('country_code')
        gadm_gpd: geopandas.GeoDataFrame of GADM codes ('gadm_code') and geometries
        gadm_code: list of gadm_codes (optional)
    Returns:
        two xarray.DataArrays with 'x' and 'y 'dimension and single band population value,
        first object is 'landscan' and second object is 'wordpop'
    """
    tif_list = os.listdir(dir_path)
    tif_list = [x for x in ['landscan', 'worldpop'] if x in set(tif_list)]
    for i in tif_list:
        raster_files = [s for s in list(os.listdir(Path(dir_path) / i)) if s.endswith('.tif')]
        if gadm_code_list:
            #input_geom = gpd.GeoSeries(gadm_gpd[gadm_gpd['gadm_code'].isin(gadm_code_list)].make_valid().unary_union)
            input_geom = gadm_gpd[gadm_gpd['gadm_code'].isin(gadm_code_list)]
            assert input_geom['geometry'].crs == '4326', 'gadm_gpd CRS not 4326 in read_tif_dir() function'
            input_geom = gpd.GeoSeries(pygeos.to_shapely(pygeos.union_all(pygeos.make_valid(pygeos.from_shapely(input_geom['geometry']))))).set_crs(4326)
        else:
            #input_geom = gpd.GeoSeries(gadm_gpd[gadm_gpd['country_code'] == country_code].make_valid().unary_union)
            input_geom = gadm_gpd[gadm_gpd['country_code'] == country_code]
            assert input_geom['geometry'].crs == '4326', 'gadm_gpd CRS not 4326 in read_tif_dir() function'
            input_geom = gpd.GeoSeries(pygeos.to_shapely(pygeos.union_all(pygeos.make_valid(pygeos.from_shapely(input_geom['geometry']))))).set_crs(4326)
        if len(raster_files) == 1:
            raster_landscan = rxr.open_rasterio(filename = Path(dir_path) / i / raster_files[0], masked = True).rio.clip(input_geom, all_touched = True, from_disk = True)      
        if len(raster_files) > 1:
            file = list(filter(re.compile(country_code).search, sorted(list(os.listdir(Path(dir_path) / 'worldpop') ))))
            raster_worldpop = rxr.open_rasterio(filename = Path(dir_path) / i / file[0], masked = True).rio.clip(input_geom, all_touched = True, from_disk = True)
    return (raster_landscan, raster_worldpop);

def merge_pixels(xr_array: xr.DataArray, xr_name: str, gpd_dataframe: gpd.GeoDataFrame) -> pd.DataFrame:  
    """Take a xarray Dataset representing a raster, map on column labels from GeoDataFrame, and return a concordance DataFrame associating pixels and labels
    Args:
        xr_array: xarray.DataArray of target raster containing dimensions 'y' and 'x'
        xr_name: str name of the raster value contained in the xr_array
        gpd_dataframe: geopandas.DataFrame containing labels and geometries to left join against the xr_array of pixels
    Returns:
        pandas.DataFrame containing rows for each unique pixel of xr_array with column labels from gpd_dataframe
    """
    raster_geoms = vectorize_pixels(xr_array = xr_array)
    xr_array.name = xr_name
    raster = xr.merge(objects = [xr_array, raster_geoms], join= 'left')
    data = raster.to_dataframe()

    data = data[data[xr_name].notnull()]
    data.reset_index(inplace=True)
    data['pixel_hash'] = data["x"].astype(str) + '_' + data["y"].astype(str)
    coord = [pygeos.box(xmin, ymin, xmax, ymax) for xmin, ymin, xmax, ymax in zip(data['xmin'], data['ymin'], data['xmax'], data['ymax'])]   

    tree = pygeos.STRtree(coord)
    match = tree.query_bulk(pygeos.from_shapely(gpd_dataframe['geometry']), predicate = 'intersects').tolist()

    gpd_columns = list(gpd_dataframe.loc[:, ~gpd_dataframe.columns.isin(['geometry'])].columns)
    match = pd.DataFrame({'index_pixel': match[1], 'index_gadm': match[0]})
    match = match.merge(data[['pixel_hash','x','y',xr_name,'ymax','ymin','xmax','xmin']], how = 'left', left_on='index_pixel', right_index=True)
    match = match.merge(gpd_dataframe[gpd_columns], how = 'left', left_on='index_gadm', right_index=True)
    match['gadm_region'] = match['gadm_code'].str.split('.').str[0] + '.' + match['gadm_code'].str.split('_|\\.').str[1]  
    match['gadm_region'].fillna(value = match['country_code'], inplace=True)
    
    pixel_group_dict = {}
    for p, g in match.groupby(by=['pixel_hash','gadm_region']).agg('count').index:
        pixel_group_dict.setdefault(p,[]).append(g)
    pixel_groups = pd.DataFrame({"pixel_hash": pixel_group_dict.keys(), "intersection_region_list": pixel_group_dict.values()})
    match = match.merge(pixel_groups, how = 'left', left_on = 'pixel_hash', right_on = 'pixel_hash')
    match['intersection_region_hash'] = [';'.join(map(str, l)) for l in match['intersection_region_list']]
    
    pixel_group_dict = {}
    for p, g in match.groupby(by=['pixel_hash','gadm_code']).agg('count').index:
        pixel_group_dict.setdefault(p,[]).append(g)
    pixel_groups = pd.DataFrame({"pixel_hash": pixel_group_dict.keys(), "intersection_gadm_list": pixel_group_dict.values()})
    match = match.merge(pixel_groups, how = 'left', left_on = 'pixel_hash', right_on = 'pixel_hash')
    match['intersection_gadm_hash'] = [';'.join(map(str, l)) for l in match['intersection_gadm_list']]
    
    match_columns = ['pixel_hash', 'x', 'y', xr_name, 'ymax', 'ymin', 'xmax', 'xmin', 'intersection_gadm_list', 'intersection_gadm_hash', 'intersection_region_list', 'intersection_region_hash'] + gpd_columns
    match = match[match_columns]

    return match

def shapelify_pixels(pd_data: pd.DataFrame) -> gpd.GeoDataFrame: 
    """Convert pandas.DataFrame containing bounding box columns to a geopandas.GeoDataFrame
    Args:
        pd_data: pandas.DataFrame returned from merge_pixels() function
    Returns:
        pandas.DataFrame converted to geopandas.GeoDataFrame with pixels formatted as geometries
    """
    coords = [pygeos.box(xmin, ymin, xmax, ymax) for xmin, ymin, xmax, ymax in zip(pd_data['xmin'], pd_data['ymin'], pd_data['xmax'], pd_data['ymax'])]   
    data = gpd.GeoDataFrame(pd_data, geometry = pygeos.to_shapely(coords), crs ="EPSG:4326")
    return data 


def allocate_population(pixel_data: pd.DataFrame, population_col: str, country_code: str, buildings_dir: Union[str, Path], blocks_dir: Union[str, Path], gadm_code_list = [], all_area = False, include_residual = False) -> pd.DataFrame: 
    """Allocate raster population in pandas.DataFrame to block geometries
    Args:
        pixel_data: pandas.DataFrame returned from merge_pixels() function
        population_col: string, name of column containing pixel's population value (e.g.,'landscan','worldpop')
        country_code: string, three-digit country code
        buildings_dir: string or Path, of top level directory containing building files at country level
        blocks_dir: string or Path, of top level directory containing building files at country level
        gadm_code_list: list, list of gadm_codes in target region, defaults to False (optional)
        all_area: boolean, keep block_id and gadm_codes outside of gadm_code_list, defaults to False (optional)
        include_residual: boolean, allocate non-intersecting residual population over entire target area, defaults to False (optional)
    Returns:
        pandas.DataFrame with following columns: 'block_id', 'gadm_code', 'country_code', str(population_col + '_population')
    """
    print(population_col)
    # Read in country buildings and blocks
    if gadm_code_list:
        buildings_full = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_points_{country_code}.parquet', memory_map = True, filters = [('gadm_code', 'in', gadm_code_list)])
        blocks_full = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True, filters = [('gadm_code', 'in', gadm_code_list)])
        pixel_map = pixel_data[~pixel_data['gadm_code'].isin(gadm_code_list)]
    else:
        buildings_full = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_points_{country_code}.parquet', memory_map = True) 
        blocks_full = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True) 
        pixel_map = pixel_data[~pixel_data['country_code'].isin([country_code])]
    
    # Dictionary of gadm_code and country_code in pixel dataframe
    pixel_map = pixel_map[['country_code','gadm_code']].drop_duplicates().reset_index(drop = True).groupby('country_code', group_keys=False)['gadm_code'].apply(list).to_dict()

    # Load buildings and blocks
    for i, j in pixel_map.items():
        # Add buildings from adjacent countries that intersect pixel area
        buildings = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_points_{i}.parquet', memory_map = True, filters = [('gadm_code', 'in', j)])
        buildings_full = pd.concat([buildings_full, buildings], ignore_index=True)
        # Add blocks from adjacent countries that intersect pixel area
        blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{i}.parquet', memory_map = True, filters = [('gadm_code', 'in', j)])
        blocks_full = pd.concat([blocks_full, blocks], ignore_index=True)

    # Check for gadm_codes that were not read but were intersected (usually b/c file not in directory)
    gadm_check_list = [item for sublist in list(pixel_map.values()) for item in sublist]
    missing_blocks_gadm_codes = [x for x in gadm_check_list if x not in set(blocks_full['gadm_code'].unique())] 
    missing_buildings_gadm_codes = [x for x in gadm_check_list if x not in set(buildings_full['gadm_code'].unique())] 
    print(f'Missing gadm_codes in block data: {missing_blocks_gadm_codes}')
    print(f'Missing gadm_codes in buildings data: {missing_buildings_gadm_codes}')

    # Spatial join block_id codes to building points
    # buildings_full = gpd.sjoin(left_df = buildings_full, right_df = blocks_full[['block_id','geometry']], how='left', predicate='intersects')
    # FIX NOT NEEDED with BLOCK ID in BUILDINGS

    # Convert pixel dataframe with corner coordinates to vector geometries
    pixel_data = pixel_data.loc[:, ~pixel_data.columns.isin(['gadm_code','country_code', 'gadm_region','geometry'])]
    pixel_data = pixel_data.loc[pixel_data.astype(str).drop_duplicates().index]
    pixel_data = shapelify_pixels(pd_data = pixel_data)

    population_total = sum(pixel_data[population_col])
    print_total_pop = str(round(population_total,2))
    print('Population in full pixel area: ', print_total_pop)

    block_population = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'building_population': pd.Series(dtype= 'float64')})

    # Check if pixel population is above 0
    if population_total > 0:
        # Spatial join building points to pixel geodataframe
        pixel_buildings_data = gpd.sjoin(left_df = pixel_data, right_df = buildings_full[['building_geohash','country_code','gadm_code','block_id','building_area','geometry']], how='left', predicate='intersects')
        # Remove pixels with null buildings and blocks
        population_buildings = pixel_buildings_data[(pixel_buildings_data['building_geohash'].notnull()) & (pixel_buildings_data['block_id'].notnull())].copy() 
    
        print_building_pixels = str(len(population_buildings['pixel_hash'].unique())) + ' out of ' + str(pixel_data.shape[0])
        print('Pixels with building intersections: ',print_building_pixels)
    
        # Allocate pixel population down to building points based on footprint area
        population_buildings['building_area_pixel'] = population_buildings.groupby(['pixel_hash'])['building_area'].transform('sum')
        population_buildings['building_population'] = population_buildings[population_col] * (population_buildings['building_area']/population_buildings['building_area_pixel'])
        # Aggregate building population to block_id level
        block_population = population_buildings[['building_geohash','country_code','gadm_code','block_id','building_area','building_area_pixel','building_population']]
        block_population = block_population.groupby(['block_id','gadm_code','country_code']).agg({'building_area': 'sum','building_population':'sum'}).reset_index()
    
        print_building_population = str(round(sum(block_population['building_population']),2)) + ' ' + str(round((sum(block_population['building_population'])/population_total)*100,3)) + '%'
        print('Population allocated by building area: ', print_building_population)
    
        # Subset residual pixels that do not intersect buildings points
        population_residual = pixel_data[~pixel_data['pixel_hash'].isin(population_buildings['pixel_hash'].unique())].copy()
        pop_residual = sum(population_residual[population_col])
        print_block_population = str(round(pop_residual,2))
        print('Population not intersecting buildings to reallocate: ',print_block_population)

        # Allocate residual of population that did not intersect buildings by block area
        if pop_residual > 0: 
            population_residual['pixel_area'] = population_residual.to_crs(3395).area
            population_residual = population_residual[['pixel_hash','x','y',population_col,'pixel_area','geometry']]
    
            print_residual_pixels = str(population_residual.shape[0]) + ' out of ' + str(pixel_data.shape[0])
            print('Residual pixels to allocate by block area: ',print_residual_pixels)
    
            # Form new geometries based on overlap between residual pixels and block_id
            #population_residual.geometry = population_residual.buffer(0)
            #blocks_full.geometry = blocks_full.buffer(0)
            block_overlay = gpd.overlay(population_residual, blocks_full, how='intersection', keep_geom_type=True) 
            block_overlay = block_overlay[block_overlay[population_col].notnull()]
            block_overlay['pixel_area_split'] = block_overlay.to_crs(3395).area
    
            # Allocate residual population by block areas that overlaps with pixel area
            block_overlay = block_overlay[block_overlay['block_id'].notnull()]
            block_overlay['pixel_area_split_total_notnull'] = block_overlay.groupby(['pixel_hash'])['pixel_area_split'].transform('sum')
            block_overlay['building_population'] = block_overlay[population_col] * (block_overlay['pixel_area_split'] / block_overlay['pixel_area_split_total_notnull'])
            
            # Aggregate block residual population to block_id level
            block_population_residual = block_overlay.groupby(['block_id','gadm_code','country_code']).agg({'building_population': 'sum'}).reset_index()
            block_population_residual['building_area'] = 0
    
            # Combine block_id population totals allocated to buildings and residual population allocated to blocks
            block_population = pd.concat([block_population, block_population_residual], ignore_index=True)
            block_population = block_population.groupby(['block_id','gadm_code','country_code']).agg({'building_area': 'sum','building_population':'sum'}).reset_index()
    
            # Filter out mis-aligned geometries (occurs rarely due to minute differences between block geometries and GADM geometries) 
            block_population['gadm_code_length'] = block_population['gadm_code'].str.len()
            block_population['extracted_gadm_code'] = block_population.apply(lambda x: str(x['block_id'])[:x['gadm_code_length']], axis=1)
            block_population[block_population['extracted_gadm_code'] == block_population['gadm_code']]
            # block_population = block_population[block_population['block_id'].str.rsplit(pat='_', n=1, expand=True)[0] == block_population['gadm_code']] 
            block_population = block_population[block_population['block_id'].str.slice(start=0, stop=3) == block_population['country_code']]
    
            # Zero out residual population allocated to block_id level with no building area and redistribute to gadm_code level
            pop_adjust = block_population[block_population['building_area'] == 0].groupby(['gadm_code']).agg({'building_population': 'sum'}).reset_index().rename(columns={"building_population":"building_population_adjust"})
            pop_residual = sum(pop_adjust['building_population_adjust'])
            # Remove blocks with no building area
            block_population = block_population.loc[(block_population['building_area'] > 0) | (block_population['building_area'].notnull())]
    
            print_allocated_population = str(round(sum(block_population['building_population']),2)) + ' ' + str(round((sum(block_population['building_population'])/population_total)*100,3)) + '%'
            print('Running sum of population allocated: ',print_allocated_population)
            print_blocks_no_buildings_population = str(round(pop_residual,2))
            print('Population allocated to blocks with no buildings: ', print_blocks_no_buildings_population)
    
            # Re-allocate residual population allocated to blocks with no buildings to the GADM level by building area
            if pop_residual > 0: 
                block_population = pd.merge(left = block_population, right = pop_adjust, how='left', on='gadm_code')
                block_population['building_population_adjust'] = block_population['building_population_adjust'].fillna(0)
                block_population['building_area_gadm'] = block_population.groupby(['gadm_code'])['building_area'].transform('sum')
                block_population['building_population'] = (block_population['building_population_adjust'] * (block_population['building_area']/block_population['building_area_gadm'])) + block_population['building_population']
    
                # Re-allocate any remaining residual population based to country scale block_id level building area shares (should be a small number from pixels touching no blocks or buildings)
                pop_residual = (population_total-sum(block_population['building_population']))
    
                print_no_blocks_no_buildings_population = str(round(pop_residual,2))
                print('Population not touching blocks or buildings: ',print_no_blocks_no_buildings_population)
                if pop_residual > 0: 
                    if include_residual is True:
                        block_population['building_population'] = ((population_total-sum(block_population['building_population'])) * (block_population['building_area']/sum(block_population['building_area']))) + block_population['building_population']   
    
    # Append block_id codes with no population
    block_population = block_population[['block_id','gadm_code','country_code','building_population']]
    if gadm_code_list:
        block_universe = blocks_full[blocks_full['gadm_code'].isin(gadm_code_list)][['block_id','gadm_code','country_code']]
    else:
        block_universe = blocks_full[blocks_full['country_code'] == country_code][['block_id','gadm_code','country_code']]        
    block_universe['building_population'] = 0
    block_universe = block_universe[~block_universe['block_id'].isin(block_population['block_id'].unique())]
    block_population = pd.concat([block_population, block_universe], ignore_index=True)
    if population_total > 0:
        print_share_overall = str(round((sum(block_population['building_population'])/population_total)*100,3)) + '%'
        print('Allocated share of population in full pixel area: ', print_share_overall)

    # Create target area flag
    if gadm_code_list:
        block_population.loc[block_population['gadm_code'].isin(gadm_code_list), 'target_area'] = int(1)
        block_population['target_area'] = block_population['target_area'].fillna(0)
        if population_total > 0:
            print_share_target = str(round((sum(block_population[block_population['target_area'] == 1]['building_population'])/population_total*100),3)) + '%'
            print('Allocated share in gadm_code_list target area: ',print_share_target)
        if all_area is False:
            block_population = block_population[block_population['target_area'] == 1]
            block_population = block_population.drop(['target_area'], axis=1)

    block_population = block_population.rename(columns={"building_population":str(population_col + '_population')})

    return block_population


def mem_profile() -> str: 
    """
    Return memory usage, str
    """
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use


def main(log_file: Path, country_chunk: list, gadm_dir: Path, blocks_dir: Path, rasters_dir: Path, targets_file: Path, buildings_dir: Path, output_dir: Path):

    # logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f"Countries to process: {country_chunk}")
    
    population_dir =  str(output_dir) + '/population'
    Path(population_dir).mkdir(parents=True, exist_ok=True)

    combined_dir =  str(output_dir) + '/combined'
    Path(combined_dir).mkdir(parents=True, exist_ok=True)

    #country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}
    #country_list = list(country_dict.keys())

    # Check files and match to countries in each source directory
    gadm_dir_list = list(filter(re.compile("gadm_").match, sorted(list(os.listdir(Path(gadm_dir))))))
    gadm_dir_list = [(re.sub('gadm_', '', re.sub('.parquet', '', i))) for i in gadm_dir_list] 

    blocks_dir_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    blocks_dir_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in blocks_dir_list] 
    
    buildings_dir_list = list(filter(re.compile("buildings_points_").match, sorted(list(os.listdir(Path(buildings_dir))))))
    buildings_dir_list = [(re.sub('buildings_points_', '', re.sub('.parquet', '', i))) for i in buildings_dir_list] 

    worldpop_dir_list = [s for s in list(os.listdir(Path(rasters_dir) / 'worldpop')) if s.endswith('.tif')]
    worldpop_dir_list = [(re.sub('.tif', '', i)) for i in worldpop_dir_list] 

    landscan_dir_list = [s for s in list(os.listdir(Path(rasters_dir) / 'landscan')) if s.endswith('.tif')]
    if len(landscan_dir_list) > 1:
        raise ValueError('More than one .tif file in ./population/landscan directory. Only one expected.')

    # Check completed files 
    output_file_list = list(filter(re.compile("population_").match, sorted(list(os.listdir(Path(population_dir))))))
    output_country_list = [(re.sub('population_', '', re.sub('.parquet', '', i))) for i in output_file_list] 
    logging.info(f"Finished countries: {output_country_list}")

    # Check country_chunk parameter against the country files across directories
    if country_chunk: 

        country_list = [x for x in country_chunk if x not in set(output_country_list)]  
        if not country_list: 
            raise ValueError('All countries in country_chunk argument already finished.')

        intersected_list = [gadm_dir_list, blocks_dir_list, buildings_dir_list, worldpop_dir_list] 
        intersected_list = list(set.intersection(*map(set,intersected_list)))

        chunk_missing = [x for x in country_list if x not in set(intersected_list)] 
        if chunk_missing:
            raise ValueError(f'Elements in country_chunk argument have missing input files: {chunk_missing}')
        
        country_list = [x for x in country_list if x in set(intersected_list)] 
        if not country_list: 
            raise ValueError('Empty country list.')
    else: 
        raise ValueError('Empty country_chunk argument.')

    if targets_file:
        pop_targets = pd.read_csv(Path(targets_file), low_memory=False)
        pop_targets = pop_targets[['ISO3_code',  'Time', 'TPopulation1July']]
        pop_targets = pop_targets[pop_targets['Time'] == 2020]
        pop_targets = pop_targets[pop_targets['ISO3_code'].isin(country_list)]
        pop_targets = pop_targets.rename(columns={'ISO3_code':'country_code','Time':'year','TPopulation1July':'un_population'}).sort_values(by=['country_code']).reset_index(drop=True)
        pop_targets['un_population'] = (pop_targets['un_population']*1000).astype(int)

    logging.info(f"Remaining countries: {country_list}")

    # Consolidate GADM data into one file
    all_gadm_gpd = gpd.GeoDataFrame({'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 
    for country_code in gadm_dir_list: 
        gadm_gpd = gpd.read_parquet(Path(gadm_dir) / f'gadm_{country_code}.parquet')
        gadm_gpd['country_code'] = country_code
        if not all(x in ['Polygon','MultiPolygon'] for x in gadm_gpd['geometry'].geom_type.unique()):
            gadm_gpd = gadm_gpd.explode(index_parts=False)
            gadm_gpd = gadm_gpd[gadm_gpd['geometry'].geom_type == 'Polygon']
            gadm_gpd = gadm_gpd.dissolve(by='gadm_code', as_index=False)
        all_gadm_gpd = pd.concat([all_gadm_gpd, gadm_gpd[['gadm_code', 'country_code', 'geometry']]], ignore_index=True)
    all_gadm_gpd['gadm_region'] = all_gadm_gpd['gadm_code'].str.split('.').str[0] + '.' + all_gadm_gpd['gadm_code'].str.split('_|\\.').str[1] 
    all_gadm_gpd['gadm_region'].fillna(value = all_gadm_gpd['country_code'], inplace=True)
    logging.info('---------')
    logging.info('---------')
    
    # Run population model
    for country_code in country_list: 
        logging.info(f"Running: {country_code}")
        c0 = time.time()
        logging.info(f"Node status: {mem_profile()}")
        country_gadm = all_gadm_gpd[all_gadm_gpd['country_code'] == country_code]
        region_list = country_gadm['gadm_region'].unique()
        blocks_pop_full = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'landscan_population': pd.Series(dtype= 'float64'), 'worldpop_population': pd.Series(dtype= 'float64')})
        
        # Run for each region within each country
        logging.info(f"Regions: {region_list}")
        for region_code in region_list:
            r0 = time.time()

            logging.info(f'Reading raster for {region_code}')
            region_gadm_list = list(country_gadm[country_gadm['gadm_region'] == region_code]['gadm_code'].unique())
            ls, wp = read_tif_dir(dir_path = rasters_dir, country_code = country_code, gadm_code_list = region_gadm_list, gadm_gpd = all_gadm_gpd)

            #logging.info('Merge pixels to GADMs')
            pixels_ls = merge_pixels(xr_array = ls, xr_name = 'landscan', gpd_dataframe = all_gadm_gpd)
            pixels_wp = merge_pixels(xr_array = wp, xr_name = 'worldpop', gpd_dataframe = all_gadm_gpd)

            #logging.info('Allocate population')
            log_file_detail = Path(os.path.dirname(log_file)) / 'deploy_3_model_population_verbose.log'
            with open(log_file_detail, 'a') as f:
                with contextlib.redirect_stdout(f):
                    a0 = time.time()
                    print('---------')
                    print(region_code)
                    blocks_ls = allocate_population(pixel_data = pixels_ls, population_col = 'landscan', country_code = country_code, buildings_dir = buildings_dir, blocks_dir = blocks_dir, gadm_code_list = region_gadm_list, all_area = False, include_residual = True)
                    blocks_wp = allocate_population(pixel_data = pixels_wp, population_col = 'worldpop', country_code = country_code, buildings_dir = buildings_dir, blocks_dir = blocks_dir, gadm_code_list = region_gadm_list, all_area = False, include_residual = True)
                    a1 = time.time()
                    print(f"{region_code}: {round((a1-a0)/60,3)} minutes")
            blocks_pop = blocks_ls.merge(right = blocks_wp, how='outer', on= ['block_id','gadm_code','country_code'])
            blocks_pop['gadm_region'] = blocks_pop['gadm_code'].str.split('.').str[0] + '.' + blocks_pop['gadm_code'].str.split('_|\\.').str[1] 
            blocks_pop['gadm_region'].fillna(value = blocks_pop['country_code'], inplace=True)
            blocks_pop = blocks_pop[blocks_pop['gadm_region'].isin([region_code])]
            blocks_pop = blocks_pop[['block_id','gadm_code','country_code','landscan_population','worldpop_population']]
            
            # Remove mis-aligned geometries (occurs rarely)
            blocks_pop['gadm_code_length'] = blocks_pop['gadm_code'].str.len()
            blocks_pop['extracted_gadm_code'] = blocks_pop.apply(lambda x: str(x['block_id'])[:x['gadm_code_length']], axis=1)
            blocks_pop[blocks_pop['extracted_gadm_code'] == blocks_pop['gadm_code']]
            #blocks_pop = blocks_pop[blocks_pop['block_id'].str.rsplit(pat = '_', n = 1, expand=True)[0] == blocks_pop['gadm_code']] 
            blocks_pop = blocks_pop[blocks_pop['block_id'].str.slice(start=0, stop = 3) == blocks_pop['country_code']]

            blocks_pop['landscan_population'] = blocks_pop['landscan_population'].fillna(0)
            blocks_pop['worldpop_population'] = blocks_pop['worldpop_population'].fillna(0)
            r1 = time.time()
            logging.info(f"{region_code}: {str(round((r1-r0)/60,3))} minutes")
            blocks_pop_full = pd.concat([blocks_pop_full, blocks_pop], ignore_index=True)
            del ls, wp, pixels_ls, pixels_wp, blocks_ls, blocks_wp, blocks_pop
        
        # Check allocation and re-scale if necessary (rare, happens with islands)
        ls_xr, wp_xr = read_tif_dir(dir_path = rasters_dir, country_code = country_code, gadm_gpd = all_gadm_gpd)
        ls_xr.name = 'landscan'
        wp_xr.name = 'worldpop'
        ls_df = ls_xr.to_dataframe()
        wp_df = wp_xr.to_dataframe()
        landscan_share = (sum(blocks_pop_full['landscan_population']) / sum(ls_df[ls_df['landscan'].notnull()]['landscan']))
        worldpop_share = (sum(blocks_pop_full['worldpop_population']) / sum(wp_df[wp_df['worldpop'].notnull()]['worldpop']))
        logging.info(f'Clipped {country_code} raster vs allocated raster:')
        logging.info('landscan: ' + str(round(landscan_share*100,4)) + f' % of ' + str(round(sum(ls_df[ls_df['landscan'].notnull()]['landscan']),3)))
        if landscan_share > 1:
            blocks_pop_full['landscan_population'] = sum(ls_df[ls_df['landscan'].notnull()]['landscan']) * ((blocks_pop_full['landscan_population'])/(sum(blocks_pop_full['landscan_population'])))
            landscan_share = (sum(blocks_pop_full['landscan_population']) / sum(ls_df[ls_df['landscan'].notnull()]['landscan']))
            logging.info('rescaled landscan: ' + str(round(landscan_share*100,4)) + f' % of ' + str(round(sum(ls_df[ls_df['landscan'].notnull()]['landscan']),3)))
        logging.info('worldpop: ' + str(round(worldpop_share*100,4)) + f' % of ' + str(round(sum(wp_df[wp_df['worldpop'].notnull()]['worldpop']),3)))
        if worldpop_share > 1:
            blocks_pop_full['worldpop_population'] = sum(wp_df[wp_df['worldpop'].notnull()]['worldpop']) * ((blocks_pop_full['worldpop_population'])/(sum(blocks_pop_full['worldpop_population'])))
            worldpop_share = (sum(blocks_pop_full['worldpop_population']) / sum(wp_df[wp_df['worldpop'].notnull()]['worldpop']))
            logging.info('rescaled worldpop: ' + str(round(worldpop_share*100,4)) + f' % of ' + str(round(sum(wp_df[wp_df['worldpop'].notnull()]['worldpop']),3)))
        c1 = time.time()
        logging.info(f"{country_code}: {str(round((c1-c0)/60,3))} minutes")

        blocks_pop_full = blocks_pop_full[['block_id', 'gadm_code', 'country_code', 'landscan_population','worldpop_population']]

        # Add population adjustment
        if targets_file:
            blocks_pop_full = pd.merge(left = blocks_pop_full, right = pop_targets, how='left', on='country_code')
            blocks_pop_full['landscan_total'] = blocks_pop_full.groupby(['country_code'])['landscan_population'].transform('sum')
            blocks_pop_full['worldpop_total'] = blocks_pop_full.groupby(['country_code'])['worldpop_population'].transform('sum')
            blocks_pop_full['landscan_population_un'] = (blocks_pop_full['landscan_population']/blocks_pop_full['landscan_total'])*blocks_pop_full['un_population']
            blocks_pop_full['worldpop_population_un'] = (blocks_pop_full['worldpop_population']/blocks_pop_full['worldpop_total'])*blocks_pop_full['un_population']
            blocks_pop_full = blocks_pop_full[['block_id', 'gadm_code', 'country_code', 'landscan_population', 'landscan_population_un', 'worldpop_population', 'worldpop_population_un']]

        # Write population file
        blocks_pop_full.to_parquet(path = Path(population_dir) / f'population_{country_code}.parquet', engine='pyarrow', compression='snappy')
        logging.info('---------')
        del blocks_pop_full, ls_xr, wp_xr, ls_df, wp_df, landscan_share, worldpop_share

    # Consolidate GADM data into one file
    logging.info(f"Consolidating countries.")
    population_output_list = list(filter(re.compile("population_").match, sorted(list(os.listdir(Path(population_dir))))))
    population_output_list = [(re.sub('population_', '', re.sub('.parquet', '', i))) for i in population_output_list] 

    if targets_file:
        all_population = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'landscan_population': pd.Series(dtype= 'float64'), 'landscan_population_un': pd.Series(dtype= 'float64'), 'worldpop_population': pd.Series(dtype= 'float64'), 'worldpop_population_un': pd.Series(dtype= 'float64')})
    else: 
        all_population = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'landscan_population': pd.Series(dtype= 'float64'), 'worldpop_population': pd.Series(dtype= 'float64')})

    for country_code in population_output_list: 
        population_country = pd.read_parquet(Path(population_dir) / f'population_{country_code}.parquet')
        all_population = pd.concat([all_population, population_country], ignore_index=True)   

    # Write the all-country file 
    logging.info(f'Writing all-country files.')
    all_population.to_parquet(Path(combined_dir) / f'all_population.parquet', compression='snappy')
    all_population.to_csv(Path(combined_dir) / f'all_population.csv')

    logging.info(f"Finished job.")
    logging.info(f"------------")
    logging.info(f"------------")

def setup(args=None):    
    parser = argparse.ArgumentParser(description='Estimate population at the street block level.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file.") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format.")
    parser.add_argument('--gadm_dir', required=True, type=Path, dest="gadm_dir", help="Path to GADM parquet directory.")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory.")
    parser.add_argument('--rasters_dir', required=True, type=Path, dest="rasters_dir", help="Path to population raster directory. Must have format of population/landscan/* and population/worldpop/*")
    parser.add_argument('--targets_file', required=False, type=Path, dest="targets_file", help="Path to file with population targets at the country level.")
    parser.add_argument('--buildings_dir', required=True, type=Path, dest="buildings_dir", help="Path to building points directory.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to top level output directory. Creates /population subdirectory.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))

