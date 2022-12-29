
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

        data_overlap = data_overlap.rename(columns={str(group_column + '_left'): group_column})
        data_overlap = data_overlap[[group_column,'geometry']]
        
        overlap_list = data_overlap[group_column].unique()
        data_overlap = data[data[group_column].isin(overlap_list)]
        
        all_intersections = [a.intersection(b) for a, b in list(itertools.combinations(data_overlap['geometry'], 2))]
        data_overlap = pd.concat([data_overlap['geometry'], gpd.GeoSeries(all_intersections).set_crs(4326)])
        data_overlap = list(shapely.ops.polygonize(data_overlap.boundary.unary_union))
        data_overlap = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(data_overlap)}).set_crs(4326).reset_index(drop=True)  
        data_overlap['geometry'] = data_overlap['geometry'].make_valid()
        data_overlap = data_overlap.assign(overlap_id = data_overlap.index.astype(int))
        
        data_overlay = gpd.overlay(df1 = data_overlap, df2 = data, how='intersection', keep_geom_type = True, make_valid = True)
        data_overlay = data_overlay.assign(area = round(data_overlay['geometry'].to_crs(3395).area,0))
        data_overlay['area_rank'] = data_overlay.groupby('overlap_id')['area'].rank(method='first', ascending=False)
        data_overlay = data_overlay[data_overlay[group_column].notnull()]
        data_overlay = data_overlay[data_overlay['area_rank'] == 1]
        data_overlay = data_overlay[column_list_id]
        
        data_sjoin = data_overlap[~data_overlap['overlap_id'].isin(data_overlay['overlap_id'].unique())]
        data_sjoin = gpd.sjoin_nearest(left_df = data_sjoin.to_crs(3395), right_df = data.to_crs(3395), how = 'left').to_crs(4326)
        data_sjoin = data_sjoin.drop_duplicates(subset=['overlap_id'], keep='first')
        data_sjoin = data_sjoin[column_list_id]
        
        data_corrected = pd.concat([data_overlay, data_sjoin], ignore_index=True)
        data_corrected = pd.merge(left = data_overlap, right = data_corrected, how='left', on='overlap_id')
        data_corrected = pd.concat([data_corrected[column_list], data[~data[group_column].isin(overlap_list)]])
        
        data_corrected['geometry'] = data_corrected['geometry'].make_valid()
        data_corrected = data_corrected.dissolve(by=column_list_no_geo, as_index=False)
        data_corrected['geometry'] = data_corrected['geometry'].make_valid()
        
        check = dask_geopandas.from_geopandas(data_corrected, npartitions = partition_count)
        check = dask_geopandas.sjoin(left = check, right = check, predicate="overlaps")
        check = check.compute()
        if check.shape[0] > 0: 
            warnings.warn(f'Unable to remove all overlaps. {check.shape[0]} overlaps remain.')
        else:
            print('All overlaps resolved.')
    else:
        print('No overlaps found.')
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

    gadm_data = gadm_data[gadm_data[gadm_column] == gadm_code]
    gadm_data = gadm_data.explode(index_parts=False)
    gadm_data = gadm_data[gadm_data['geometry'].geom_type == 'Polygon']

    gadm_array = pygeos.from_shapely(gadm_data['geometry'])
    gadm_array = pygeos.multipolygons(gadm_array)
    osm_data = pygeos.line_merge(pygeos.intersection(pygeos.multilinestrings(osm_data),gadm_array))
    gadm_lines = pygeos.line_merge(pygeos.multilinestrings(pygeos.get_exterior_ring(pygeos.get_parts(gadm_array))))

    polys = pygeos.polygonize_full([pygeos.union(osm_data,gadm_lines)])[0]
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
    
    if math.ceil(block_bulk.shape[0]/5000) > 10: num_partitions = math.ceil(block_bulk.shape[0]/5000)
    else: num_partitions = 10
    gadm_blocks = remove_overlaps(data = gadm_blocks, group_column = 'block_id', npartitions = num_partitions)
    gadm_blocks['geometry'] = gadm_blocks['geometry'].make_valid()
    gadm_blocks = gadm_blocks.explode(index_parts=False)
    gadm_blocks = gadm_blocks[gadm_blocks.geom_type == "Polygon"]
    gadm_blocks = gadm_blocks[round(gadm_blocks['geometry'].to_crs(3395).area,0) > 0]
    gadm_blocks = gadm_blocks.assign(block_id = [gadm_code + '_' + str(x) for x in list(gadm_blocks.index)])
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS.")
        representative_pt = gadm_blocks.geometry.representative_point().to_list()
    precision_level = 18
    gadm_blocks['block_geohash'] = list(map(lambda x: pygeohash.encode(x.x,x.y, precision=precision_level), representative_pt))
    gadm_blocks = gadm_blocks[['block_id','block_geohash','gadm_code','country_code','geometry']].reset_index(drop = True).to_crs(epsg=4326)
    
    return gadm_blocks

def mem_profile() -> str: 
    """
    Return memory usage, str
    """
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use


def main(log_file: Path, country_chunk: list, osm_dir: Path, gadm_dir: Path, output_dir: Path):
    
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    # Make directories
    block_dir =  str(output_dir) + '/blocks'
    Path(block_dir).mkdir(parents=True, exist_ok=True)

    block_gpkg_dir =  str(output_dir) + '/blocks/gpkg'
    Path(block_gpkg_dir).mkdir(parents=True, exist_ok=True)
    
    block_overlaps_dir =  str(output_dir) + '/blocks/gpkg/overlaps'
    Path(block_overlaps_dir).mkdir(parents=True, exist_ok=True)

    street_dir =  str(output_dir) + '/streets'
    Path(street_dir).mkdir(parents=True, exist_ok=True)

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

    # Iterate through country codes
    for country_code in country_list: 
        logging.info(f"Processing: {country_code}")
        print(country_code)

        # Read OSM files
        osm_gpd = gpd.read_parquet(Path(osm_dir) / f'{country_code}-linestring.parquet')
        osm_gpd = osm_gpd.explode(ignore_index = True)
        osm_pygeos = from_shapely_srid(geometry = osm_gpd, srid = 4326) 
    
        # Read GADM files
        gadm_gpd = gpd.read_parquet(Path(gadm_dir) / f'gadm_{country_code}.parquet')
        #gadm_gpd = gpd.read_parquet(Path(gadm_dir) / 'all_gadm.parquet', memory_map = True, filters = [('country_code', 'in', [country_code])])
        
        # Write street (non-footpath) geometries
        osm_streets = osm_gpd[osm_gpd['highway'].notnull()]
        osm_streets = osm_streets[~osm_streets['highway'].isin(['footway', 'bridleway', 'steps', 'corridor', 'path', 'cycleway'])]
        osm_streets.to_parquet(Path(street_dir) / f'streets_{country_code}.parquet', compression='snappy')
        del osm_gpd
        
        # GADM list
        gadm_list = list(gadm_gpd['gadm_code'].unique())
        
        # Initialize GeoDataFrame
        block_bulk = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'block_geohash': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 

        # Build blocks for each GADM
        logging.info(f"Build blocks")
        t0 = time.time()
        output_map = list(map(lambda x: build_blocks(gadm_data = gadm_gpd, osm_data = osm_pygeos, gadm_column = 'gadm_code', gadm_code = x), gadm_list))
        for i,j in enumerate(output_map): block_bulk = pd.concat([block_bulk, output_map[i]], ignore_index=True)
        t1 = time.time()
        logging.info(f"build_blocks() finished: {block_bulk.shape}, {mem_profile()}, {str(round(t1-t0,3)/60)} minutes")
        
        # Make valid
        block_bulk['geometry'] = block_bulk['geometry'].make_valid()
        
        # Report overlaps
        if math.ceil(block_bulk.shape[0]/5000) > 10: num_partitions = math.ceil(block_bulk.shape[0]/5000)
        else: num_partitions = 10
        check = dask_geopandas.from_geopandas(block_bulk, npartitions = num_partitions)
        check = dask_geopandas.sjoin(left = check, right = check, predicate="overlaps")
        check = check.compute()
        if check.shape[0] > 0: 
            logging.info(f"Number of overlaps: {check[0].shape}")
            check.to_file(Path(block_overlaps_dir) / f'blocks_overlaps_{country_code}.gpkg', driver="GPKG")

        # Write block geometries
        block_bulk.to_parquet(Path(block_dir) / f'blocks_{country_code}.parquet', compression='snappy')
        block_bulk.to_file(Path(block_gpkg_dir) / f'blocks_{country_code}.gpkg', driver="GPKG")
        logging.info(f"Finished {country_code}.")

    logging.info(f"Finished.")

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


