
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
from typing import Union
gpd.options.use_pygeos = True

import shapely.ops
import itertools

from pathlib import Path
import psutil
import warnings
import logging
import contextlib
import io
import sys
import time
import math
import re
import os
import argparse
import dask_geopandas
import pygeohash
np.set_printoptions(suppress=True)
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')


def transform_crs(geometry: pygeos.Geometry, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395") -> pygeos.Geometry:
    """
    Transform PyGEOS Geometry from one coordinate 
    reference system to another
    Args:
        geometry: PyGEOS Geometry array
        epsg_from: string, anything accepted by pyproj.CRS.from_user_input(), (i.e., "EPSG:4326")
        epsg_to: string, anything accepted by pyproj.CRS.from_user_input(), (i.e., "EPSG:4326")
    Returns:
        PyGEOS Geometry array
    """
    transformer = pyproj.Transformer.from_crs(crs_from = epsg_from, crs_to = epsg_to, always_xy=True)
    coords = pygeos.get_coordinates(geometry)
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])
    data = pygeos.set_coordinates(geometry.copy(), np.array(new_coords).T)
    return data

def from_shapely_srid(geometry: gpd.GeoDataFrame, srid: int) -> pygeos.Geometry:
    """
    Transform GeoDataFrame to a PyGEOS Geometry array 
    and set spatial reference identifier (SRID).
    Args:
        geometry: GeoSeries or GeoDataFrame 
        srid: integer representing spatial reference identifier (i.e., 4326)
    Returns:
        PyGEOS Geometry array
    """
    data = pygeos.set_srid(pygeos.from_shapely(geometry['geometry'].to_crs(epsg=srid)), srid)
    return data

def compute_area(geometry: pygeos.Geometry, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395") -> np.array:
    """
    Transform a PyGEOS Geometry from one coordinate 
    reference system to another and calculate area
    Args:
        geometry: PyGEOS Geometry array
        epsg_from: string, anything accepted by pyproj.CRS.from_user_input(), (i.e., "EPSG:4326")
        epsg_to: string, anything accepted by pyproj.CRS.from_user_input(), (i.e., "EPSG:4326")
    Returns:
        Numpy array of areas in meters square
    """
    transformer = pyproj.Transformer.from_crs(crs_from = epsg_from, crs_to = epsg_to, always_xy=True)
    coords = pygeos.get_coordinates(geometry)
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])
    data = pygeos.area(pygeos.set_coordinates(geometry.copy(), np.array(new_coords).T))
    return data

def remove_overlaps(data: gpd.GeoDataFrame, group_column: str, partition_count: int = 10) -> gpd.GeoDataFrame:
    """ 
    Args:
        data: GeoDataFrame, containing delineations, requires CRS WGS 84 EPSG 4326.
        group_column: str, column name with label uniquely identifying each row geometry.
        partition_count: number of partitions to use with dask_geopandas.sjoin(), defaults to 10.
    Returns:
        GeoDataFrame with all overlapping geometries removed. 
        Assigns overlapping sections to the geometry group_column with most area in common or is closest.
    """
    assert data.crs == 'epsg:4326', "data is not epsg:4326."
    column_list = list(data.columns)

    column_list_no_geo = list(data.columns)
    column_list_no_geo.remove('geometry')

    column_list_id = list(data.columns)
    column_list_id.remove('geometry')
    column_list_id.append('overlap_id')

    data_overlap = dask_geopandas.from_geopandas(data, npartitions = partition_count)
    data_overlap = dask_geopandas.sjoin(left = data_overlap, right = data_overlap, predicate="overlaps")
    data_overlap = data_overlap.compute()

    if data_overlap.shape[0] > 0: 
        print(f'{data_overlap.shape[0]} overlaps.')
        logging.info(f'{data_overlap.shape[0]} overlaps.')

        data_overlap = data_overlap.rename(columns={str(group_column + '_left'): group_column})
        data_overlap = data_overlap[[group_column,'geometry']]
        
        # Filter to overlapping areas
        overlap_list = data_overlap[group_column].unique()
        data_overlap = data[data[group_column].isin(overlap_list)][[group_column,'geometry']]

        # Resolve overlaps via intersection and polygonization
        overlap_array = pygeos.union_all([pygeos.line_merge(pygeos.boundary(pygeos.from_shapely(data_overlap['geometry'])))])
        overlap_array = pygeos.polygonize_full([overlap_array])[0]
        overlap_array = pygeos.get_parts(pygeos.normalize(pygeos.get_parts(overlap_array))) 
        overlap_array = pygeos.make_valid(overlap_array)

        data_overlap = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(pygeos.to_shapely(overlap_array))}).set_crs(4326).reset_index(drop=True)  
        data_overlap['geometry'] = data_overlap['geometry'].make_valid()
        data_overlap = gpd.overlay(df1 = data_overlap, df2 = data[~data[group_column].isin(overlap_list)], how = 'difference', keep_geom_type = True, make_valid = True)
        data_overlap = data_overlap.assign(overlap_id = data_overlap.index.astype(int))
        
        # Use overlay to relabel overlaps defaulting to label with largest area
        data_overlay = gpd.overlay(df1 = data_overlap, df2 = data, how='intersection', keep_geom_type = True, make_valid = True)
        data_overlay = data_overlay.assign(area = round(data_overlay['geometry'].to_crs(3395).area,0))
        data_overlay['area_rank'] = data_overlay.groupby('overlap_id')['area'].rank(method='first', ascending=False)
        data_overlay = data_overlay[data_overlay[group_column].notnull()]
        data_overlay = data_overlay[data_overlay['area_rank'] == 1]
        data_corrected = data_overlay[column_list_id]
                
        # Merge in labels
        data_corrected = pd.merge(left = data_overlap, right = data_corrected, how='left', on='overlap_id')
        data_corrected = pd.concat([data_corrected[column_list], data[~data[group_column].isin(overlap_list)]])
        
        data_corrected['geometry'] = data_corrected['geometry'].make_valid()
        data_corrected = data_corrected.dissolve(by=column_list_no_geo, as_index=False)
        data_corrected['geometry'] = data_corrected['geometry'].make_valid()
        
        check = dask_geopandas.from_geopandas(data_corrected, npartitions = partition_count)
        check = dask_geopandas.sjoin(left = check, right = check, predicate="overlaps")
        check = check.compute()

        if check.shape[0] > 0: 
            print(f'Unable to resolve all overlaps. {check.shape[0]} overlaps remain.')
            logging.warning(f'Unable to resolve all overlaps. {check.shape[0]} overlaps remain.')
            check = check[[group_column + '_left','geometry']].rename(columns={str(group_column + '_left'): group_column})
            unresolved_area = gpd.overlay(df1 = check, df2 = check, how='intersection', keep_geom_type = True, make_valid = True)
            unresolved_area = sum(unresolved_area[unresolved_area[str(group_column + '_1')] != unresolved_area[str(group_column + '_2')]].to_crs(3395).area)*1e-6
            print(f'Unresolvable overlapping area: {unresolved_area} km2')
            logging.warning(f'Unresolvable overlapping area: {unresolved_area} km2')
        else:
            print('All overlaps resolved.')
            logging.info('All overlaps resolved.')
    else:
        print('No overlaps found.')
        logging.info('No overlaps found.')

        data_corrected = data

    return data_corrected 
    

def build_blocks(gadm_data: gpd.GeoDataFrame, osm_data: Union[pygeos.Geometry, gpd.GeoDataFrame], gadm_column: str, gadm_code: str) -> gpd.GeoDataFrame:
    """
    For a given district level GADM in a country, this 
    subdivides the area enclosed by the GADM boundary and generates 
    polygons based on areas enclosed by OSM linestrings. 
    Args:
        gadm_data: GeoDataFrame, containing GADM delineations (downloaded from https://gadm.org/data.html), requires CRS WGS 84 EPSG 4326
        osm_data: PyGEOS Geometry array (or GeoDataFrame with slower performance), of linestrings that delineate blocks within GADMs, requires CRS WGS 84 EPSG 4326
        gadm_column: string, name of column containing gadm_code
        gadm_code: string, unique identification code for each GADM in gadm_column
    Returns:
        GeoDataFrame with enclosed block geometries for a given 'gadm_code', each block is assigned a unique 'block_id'.
        To alphanumeric suffix in the block_id is a geohash that can be decoded into a lat-lon using # http://geohash.co/
        Geometry projected in CRS WGS 84 EPSG 4326.
    """
    assert gadm_data.crs == 'epsg:4326', "gadm_data is not epsg:4326."
    if type(osm_data) == gpd.geodataframe.GeoDataFrame: 
        osm_data = from_shapely_srid(geometry = osm_data, srid = 4326)
    assert len(np.unique(pygeos.get_srid(osm_data))) == 1 and np.unique(pygeos.get_srid(osm_data))[0] == 4326, "osm_data is not epsg:4326."
    print(gadm_code)
    logging.info(f'{gadm_code}')

    gadm_data = gadm_data[gadm_data[gadm_column] == gadm_code]
    gadm_data = gadm_data.explode(index_parts=False)
    gadm_data = gadm_data[gadm_data['geometry'].geom_type == 'Polygon']

    gadm_array = pygeos.from_shapely(gadm_data['geometry'])
    gadm_array = pygeos.multipolygons(gadm_array)
    osm_data = pygeos.line_merge(pygeos.intersection_all([pygeos.multilinestrings(osm_data),gadm_array]))
    gadm_lines = pygeos.line_merge(pygeos.multilinestrings(pygeos.get_exterior_ring(pygeos.get_parts(gadm_array))))

    polys = pygeos.polygonize_full([pygeos.union_all([osm_data,gadm_lines])])[0]
    polys = pygeos.get_parts(pygeos.normalize(pygeos.get_parts(polys))) 
    polys = pygeos.make_valid(polys)

    gadm_blocks = gpd.GeoDataFrame.from_dict({"country_code": gadm_code[0:3],"gadm_code": gadm_code,'geometry': pygeos.to_shapely(polys)}).set_crs(4326)  
    gadm_blocks = gadm_blocks.reset_index(drop=True)

    gadm_blocks = gpd.overlay(df1 = gadm_blocks, df2 = gadm_data[['geometry']], how='intersection', keep_geom_type = True, make_valid = True)
    gadm_blocks['geometry'] = gadm_blocks['geometry'].make_valid()
    gadm_blocks = gadm_blocks.explode(index_parts=False)
    gadm_blocks = gadm_blocks[gadm_blocks.geom_type == "Polygon"]
    gadm_blocks = gadm_blocks[round(gadm_blocks['geometry'].to_crs(3395).area,0) > 0]
    gadm_blocks = gadm_blocks.assign(block_id = [gadm_code + '_' + str(x) for x in list(gadm_blocks.index)])
    
    if math.ceil(gadm_blocks.shape[0]/10000) > 1: num_partitions = math.ceil(gadm_blocks.shape[0]/10000)
    else: num_partitions = 1
    gadm_blocks = remove_overlaps(data = gadm_blocks, group_column = 'block_id', partition_count = num_partitions)
    gadm_blocks['geometry'] = gadm_blocks['geometry'].make_valid()
    gadm_blocks = gadm_blocks.explode(index_parts=False)
    gadm_blocks = gadm_blocks[gadm_blocks.geom_type == "Polygon"]
    gadm_blocks = gadm_blocks[round(gadm_blocks['geometry'].to_crs(3395).area,0) > 0]

    input_area_check = sum(gadm_data['geometry'].to_crs(3395).area)*0.0001
    output_area_check = sum(gadm_blocks['geometry'].to_crs(3395).area)*0.0001
    residual_area_check = input_area_check - output_area_check 

    if residual_area_check >= 1: 
        residual_area = gadm_data['geometry'].unary_union.difference(gadm_blocks['geometry'].unary_union)
        residual_area = gpd.GeoDataFrame.from_dict({"country_code": gadm_code[0:3],"gadm_code": gadm_code,'geometry': gpd.GeoSeries(residual_area)}).set_crs(4326)  
        residual_area = residual_area.explode(index_parts = False)
        residual_area = residual_area[residual_area.geom_type == "Polygon"]
        residual_area = residual_area[residual_area.to_crs(3395).area*0.0001 >= 1]
        if residual_area.shape[0] > 0:
            print(f'Hectares to add {round(sum(residual_area.to_crs(3395).area)*0.0001,2)}.')
            logging.warning(f'Hectares to add {round(sum(residual_area.to_crs(3395).area)*0.0001,2)}.')
            gadm_blocks = pd.concat([gadm_blocks[['country_code','gadm_code','geometry']], residual_area[['country_code','gadm_code','geometry']]], ignore_index=True)
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS.")
        representative_pt = gadm_blocks.geometry.representative_point().to_list()
    precision_level = 18
    gadm_blocks['block_geohash'] = list(map(lambda x: pygeohash.encode(x.x,x.y, precision=precision_level), representative_pt))

    gadm_blocks = gadm_blocks.sort_values(by='block_geohash', ascending=False).reset_index(drop=True)
    gadm_blocks = gadm_blocks.assign(block_id = [gadm_code + '_' + str(x) for x in list(gadm_blocks.index)])
    gadm_blocks = gadm_blocks[['block_id','block_geohash','gadm_code','country_code','geometry']].reset_index(drop = True).to_crs(epsg=4326)
    
    return gadm_blocks


def mem_profile() -> str: 
    """
    Return memory usage, str
    """
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use


def main(log_file: Path, country_chunk: list, osm_dir: Path, gadm_dir: Path, output_dir: Path):
    
    #logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    # Make directories
    block_dir =  str(output_dir) + '/blocks'
    Path(block_dir).mkdir(parents=True, exist_ok=True)

    block_gpkg_dir =  str(output_dir) + '/blocks/gpkg'
    Path(block_gpkg_dir).mkdir(parents=True, exist_ok=True)
    
    #block_overlaps_dir =  str(output_dir) + '/blocks/gpkg/overlaps'
    #Path(block_overlaps_dir).mkdir(parents=True, exist_ok=True)

    street_dir =  str(output_dir) + '/streets'
    Path(street_dir).mkdir(parents=True, exist_ok=True)

    combined_dir =  str(output_dir) + '/combined'
    Path(combined_dir).mkdir(parents=True, exist_ok=True)

    logging.info(f"block_dir: {block_dir}")
    logging.info(f"street_dir: {street_dir}")

    # OSM file check
    osm_inputs_lists = list(filter(re.compile("-linestring.parquet").search, sorted(list(os.listdir(Path(osm_dir))))))
    osm_inputs_lists = [re.sub('-linestring.parquet', '', i) for i in osm_inputs_lists]
    in_chunk_not_in_osm_inputs = [x for x in country_chunk if x not in set(osm_inputs_lists)]
    if len(in_chunk_not_in_osm_inputs) > 0:
        raise ValueError(f'OSM input data does not exist for {in_chunk_not_in_osm_inputs} in country_chunk argument.')

    # GADM file check
    gadm_inputs_list = list(filter(re.compile("gadm_").match, sorted(list(os.listdir(Path(gadm_dir))))))
    gadm_inputs_list = [(re.sub('gadm_', '', re.sub('.parquet', '', i))) for i in gadm_inputs_list] 
    in_chunk_not_in_gadm_inputs = [x for x in country_chunk if x not in set(gadm_inputs_list)]
    if len(in_chunk_not_in_gadm_inputs) > 0:
        raise ValueError(f'GADM input data does not exist for {in_chunk_not_in_gadm_inputs} in country_chunk argument.')

    logging.info(f"------------")
    logging.info(f"------------")

    # Check for completed countries in output directory
    output_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(block_dir))))))
    output_country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in output_file_list] 
    logging.info(f"Finished countries: {output_country_list}")
    
    # Remove completed countries
    if country_chunk: 
        country_list = [x for x in country_chunk if x not in set(output_country_list)]
        if not country_list: 
            raise ValueError('All countries in country_chunk argument already finished.')
    else: 
        raise ValueError('Empty country_chunk arg.')

    logging.info(f"Countries to process: {country_list}")
    logging.info(f"Generate block geometries")
    logging.info(f'-------------')

    # Iterate through country codes
    for country_code in country_list: 
        logging.info(f"Processing: {country_code}")
        print(country_code)

        # Read OSM files
        osm_gpd = gpd.read_parquet(Path(osm_dir) / f'{country_code}-linestring.parquet')
        osm_gpd = osm_gpd.explode(ignore_index = True)
        osm_gpd = osm_gpd[~osm_gpd['highway'].isin(['footway', 'bridleway', 'steps', 'corridor', 'path', 'cycleway'])]
        osm_pygeos = from_shapely_srid(geometry = osm_gpd, srid = 4326) 

        # Write street (non-footpath) geometries
        osm_streets = osm_gpd[osm_gpd['highway'].notnull()]
        osm_streets.to_parquet(Path(street_dir) / f'streets_{country_code}.parquet', compression='snappy')
        del osm_gpd, osm_streets

        # Read GADM files
        gadm_gpd = gpd.read_parquet(Path(gadm_dir) / f'gadm_{country_code}.parquet')
        #gadm_gpd = gpd.read_parquet(Path(gadm_dir) / 'all_gadm.parquet', memory_map = True, filters = [('country_code', 'in', [country_code])])
        
        # GADM list
        gadm_list = list(gadm_gpd['gadm_code'].unique())
        
        # Initialize GeoDataFrame
        block_bulk = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'block_geohash': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': gpd.GeoSeries(dtype='geometry')}).set_crs(epsg=4326) 

        # Build blocks for each GADM
        logging.info(f"Build blocks")
        t0 = time.time()
        output_map = list(map(lambda x: build_blocks(gadm_data = gadm_gpd, osm_data = osm_pygeos, gadm_column = 'gadm_code', gadm_code = x), gadm_list))
        for i,j in enumerate(output_map): block_bulk = pd.concat([block_bulk, output_map[i]], ignore_index=True)
        t1 = time.time()
        logging.info(f"Blocking succeeded: {block_bulk.shape}, {mem_profile()}, {round((t1-t0)/60,2)} minutes")

        # Make valid 
        block_bulk['geometry'] = block_bulk['geometry'].make_valid()

        # Area, perimeter 
        block_bulk['block_area'] = block_bulk['geometry'].to_crs(3395).area
        block_bulk['block_perimeter'] = block_bulk['geometry'].to_crs(3395).length
        block_bulk = block_bulk[['block_id','block_geohash','gadm_code','country_code','block_area','block_perimeter','geometry']]

        # Check area of input and output geometries
        input_area = round(sum(gadm_gpd.to_crs(3395).area)*1e-6,2)
        output_area = round(sum(block_bulk.to_crs(3395).area)*1e-6,2)
        logging.info(f'Input area: {input_area} km2')
        logging.info(f'Output area: {output_area} km2')
        logging.info(f'Difference: {round((output_area - input_area),2)} km2')
        logging.info(f'Pct diff: {(output_area - input_area) / input_area}')

        # # Perform a final overlap correction 
        # if (((output_area - input_area) / input_area) > 0.001) or ((output_area - input_area) > 1):
        #     if math.ceil(block_bulk.shape[0]/10000) > 1: num_partitions = math.ceil(block_bulk.shape[0]/10000)
        #     else: num_partitions = 1
        #     f = io.StringIO()
        #     with contextlib.redirect_stdout(f):
        #         block_bulk = remove_overlaps(data = block_bulk, group_column = 'block_id', partition_count = num_partitions)
        #         overlap_log = f.getvalue().replace('\n', ' ')
        #         logging.info(f'Overlap correction: {overlap_log}')

            # # Check area of input and output geometries
            # input_area = round(sum(gadm_gpd.to_crs(3395).area)*1e-6,2)
            # output_area = round(sum(block_bulk.to_crs(3395).area)*1e-6,2)
            # logging.info(f'Input area: {input_area} km2')
            # logging.info(f'Output area: {output_area} km2')
            # logging.info(f'Difference: {round((output_area - input_area),2)} km2')
            # logging.info(f'Pct diff: {(output_area - input_area) / input_area}')

            # # Report overlaps if they exist
            # check = dask_geopandas.from_geopandas(block_bulk, npartitions = num_partitions)
            # check = dask_geopandas.sjoin(left = check, right = check, predicate="overlaps")
            # check = check.compute()
            # if check.shape[0] > 0: 
            #     print(check.shape[0])
            #     logging.info(f'Number of unresolvable countrywide overlaps: {check.shape[0]}')
            #     check.to_file(Path(block_overlaps_dir) / f'blocks_overlaps_{country_code}.gpkg', driver="GPKG")

        # Write block geometries
        block_bulk.to_parquet(Path(block_dir) / f'blocks_{country_code}.parquet', compression='snappy')
        block_bulk.to_file(Path(block_gpkg_dir) / f'blocks_{country_code}.gpkg', driver="GPKG")
        logging.info(f"Finished {country_code}.")
        logging.info(f'-------------')

    # Consolidate GADM data into one file
    logging.info(f"Consolidating countries.")
    blocks_output_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(block_dir))))))
    blocks_output_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in blocks_output_list] 
    
    all_blocks = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'block_geohash': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'block_area': pd.Series(dtype='float64'), 'block_perimeter': pd.Series(dtype='float64'), 'geometry': gpd.GeoSeries(dtype='geometry')}).set_crs(epsg=4326) 
    for country_code in blocks_output_list: 
        block_country = gpd.read_parquet(Path(block_dir) / f'blocks_{country_code}.parquet')
        all_blocks = pd.concat([all_blocks, block_country], ignore_index=True)   

    # Write the all-country file 
    logging.info(f'Writing all-country files.')
    all_blocks.to_parquet(Path(combined_dir) / f'all_blocks.parquet', compression='snappy')
    all_blocks.to_file(Path(combined_dir) / f'all_blocks.gpkg', driver="GPKG")

    logging.info(f"Finished job.")
    logging.info(f"------------")

def setup(args=None):    
    parser = argparse.ArgumentParser(description='Build blocks geometries.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file.") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format.")
    parser.add_argument('--osm_dir', required=True, type=Path, dest="osm_dir", help="Path to OSM parquet directory.")
    parser.add_argument('--gadm_dir', required=True, type=Path, dest="gadm_dir", help="Path to GADM parquet directory.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to top level output directory. Creates /blocks and /streets subdirectories.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))


