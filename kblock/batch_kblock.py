
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
from typing import List, Union
gpd.options.use_pygeos = True

from pandas._libs.lib import is_integer
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
import functools
import argparse
from itertools import chain
import math

import pyarrow
import multiprocessing
import dask 
import dask_geopandas

from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

def index_buildings(gadm_block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Joins GeoDataFrame of enclosed block geometries to a GeoDataFrame containing building footprints, for footprints that 
    overlap multiple with block geometries this function allocates the footprint to the block containing its centroid.
    Args:
        gadm_block_data: GeoDataFrame output returned from build_blocks() function, requires CRS WGS 84 EPSG 4326 or 3395
        bldg_data: GeoDataFrame containing building geometries, requires CRS WGS 84 EPSG 4326 or 3395 (accepts 'Polygon' or 'Point' geometries)
    Returns:
        GeoDataFrame with building geometries mapped to 'block_id','gadm_code','country_code'.
        Geometry projected in CRS WGS 84 EPSG 4326. 
    """
    assert gadm_block_data.crs in ['epsg:3395','epsg:4326'], "gadm_block_data is not epsg:4326 or epsg:3395."
    assert bldg_data.crs in ['epsg:3395','epsg:4326'], "bldg_data is not epsg:4326 or epsg:3395."

    if gadm_block_data.crs == 'epsg:4326': gadm_block_data = gadm_block_data.to_crs(epsg=3395)
    if bldg_data.crs == 'epsg:4326': bldg_data = bldg_data.to_crs(epsg=3395)

    gtypes = bldg_data['geometry'].geom_type.unique()
    if len(gtypes) > 1 or gtypes[0] != 'Point':
        bldg_points = bldg_data.copy()
        bldg_points['geometry'] = bldg_points.centroid
    else:
        bldg_points = bldg_data.copy()

    bldg_index = bldg_points.sindex

    index_bulk = bldg_index.query_bulk(gadm_block_data['geometry'], predicate="intersects")  ##
    blocks_buildings_map = pd.DataFrame({'index_blocks': index_bulk[0], 'index_buildings': index_bulk[1]})
    blocks_buildings_map = blocks_buildings_map.merge(bldg_data[['geometry']], how = 'left', left_on='index_buildings', right_index=True)
    blocks_buildings_map = blocks_buildings_map.merge(gadm_block_data[['block_id','gadm_code','country_code']], left_on='index_blocks', right_index=True)

    data = gpd.GeoDataFrame(blocks_buildings_map[['block_id','gadm_code','country_code','geometry']]).set_crs(epsg=3395)
    data = data.to_crs(epsg=4326)
    return data

def index_building_points(gadm_block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Joins GeoDataFrame of enclosed block geometries to a GeoDataFrame containing building footprint centroids
    Args:
        gadm_block_data: GeoDataFrame output returned from build_blocks() function, requires CRS WGS 84 EPSG 4326 or 3395
        bldg_data: GeoDataFrame containing building point geometries with building_area column, requires CRS WGS 84 EPSG 4326 or 3395 (accepts 'Polygon' or 'Point' geometries)
    Returns:
        GeoDataFrame with building geometries mapped to 'building_id','building_area','block_id','gadm_code','country_code'.
        Geometry projected in CRS WGS 84 EPSG 4326. 
    """
    assert gadm_block_data.crs in ['epsg:3395','epsg:4326'], "gadm_block_data is not epsg:4326 or epsg:3395."
    assert bldg_data.crs in ['epsg:3395','epsg:4326'], "bldg_data is not epsg:4326 or epsg:3395."

    if gadm_block_data.crs == 'epsg:4326': gadm_block_data = gadm_block_data.to_crs(epsg=3395)
    if bldg_data.crs == 'epsg:4326': bldg_data = bldg_data.to_crs(epsg=3395)

    gtypes = bldg_data['geometry'].geom_type.unique()
    if len(gtypes) > 1 or gtypes[0] != 'Point':
        bldg_points = bldg_data.copy()
        bldg_points['geometry'] = bldg_points.centroid
    else:
        bldg_points = bldg_data.copy()

    bldg_index = bldg_points.sindex

    index_bulk = bldg_index.query_bulk(gadm_block_data['geometry'], predicate="intersects")  ##
    blocks_buildings_map = pd.DataFrame({'index_blocks': index_bulk[0], 'index_buildings': index_bulk[1]})
    blocks_buildings_map = blocks_buildings_map.merge(bldg_data[['building_id','building_area','geometry']], how = 'left', left_on='index_buildings', right_index=True)
    blocks_buildings_map = blocks_buildings_map.merge(gadm_block_data[['block_id','gadm_code','country_code']], left_on='index_blocks', right_index=True)

    data = gpd.GeoDataFrame(blocks_buildings_map[['building_id','building_area','block_id','gadm_code','country_code','geometry']]).set_crs(epsg=3395)
    data = data.to_crs(epsg=4326)
    return data

def compute_k(block_id: str, block_col: str, block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame, street_linestrings: Union[pygeos.Geometry, gpd.GeoDataFrame] = None, buffer_radius: float = 100, include_geometry: bool = True) -> Union[gpd.GeoDataFrame, pd.DataFrame]:
    """
    Computes the k-complexity value for each block as well as block level metadata taking a GeoDataFrame of enclosed blocks, 
    the indexed building GeoDataFrame, and street linestrings in the form of a PyGEOS Geometry (or GeoDataFrame with slower performance).
    Args:
        block_id: string, unique identification code in block_col
        block_col: string, column name that contains the block_id codes (present in block_data and bldg_data)
        block_data: GeoDataFrame output returned from build_blocks() function, requires CRS 4326 or 3395 for faster performance
        bldg_data: GeoDataFrame output returned from index_buildings() or index_building_points() function, requires CRS 4326 or 3395 for faster performance
        street_linestrings: GeoDataFrame or PyGEOS Geometry array (optional), linestrings representing street networks, requires CRS 4326 or 3395 for faster performance
        buffer_radius: float (optional), extends disconnected street_linestrings to street network that are within a set tolerance in meters
        include_geometry: bool (optional), default True, if False the block geometry will not be included and function will return a Pandas DataFrame
    Returns:
        GeoDataFrame with block geometries in row and the columns: 'block_id','block_area','building_area','building_count','building_layers','k_complexity'
        Geometry projected in CRS 4326.
    """
    assert block_data.crs in ['epsg:3395','epsg:4326'], "block_data is not epsg:4326 or epsg:3395."
    assert bldg_data.crs in ['epsg:3395','epsg:4326'], "bldg_data is not epsg:4326 or epsg:3395."
    if block_data.crs == 'epsg:4326': block_data = block_data.to_crs(epsg=3395)
    if bldg_data.crs == 'epsg:4326': bldg_data = bldg_data.to_crs(epsg=3395)

    block_layers = [] 
    
    bldg_data = bldg_data.loc[bldg_data[block_col] == block_id].copy() 
    block_data = block_data.loc[block_data[block_col] == block_id].copy() 
    block = pygeos.from_shapely(block_data['geometry'])

    country_code = block_data['country_code'].unique()[0] 
    gadm_code = block_data['gadm_code'].unique()[0] 
    block_area = round(pygeos.area(pygeos.multipolygons(pygeos.from_shapely(block_data['geometry']))),2)
    bldg_array = pygeos.from_shapely(bldg_data['geometry'])
    bldg_count = np.sum(pygeos.get_num_geometries(bldg_array))
    is_connected = False

    if 'building_area' in bldg_data.columns:
        bldg_area = round(sum(bldg_data['building_area']),2)
    else:     
        bldg_area_list = pygeos.area(bldg_array)
        bldg_area = round(sum(bldg_area_list),2)

    if bldg_count > 0:
        gtypes = bldg_data['geometry'].geom_type.unique()
        if len(gtypes) > 1 or gtypes[0] != 'Point':
            bldg_data["geometry"] = bldg_data.centroid
        building_points = pygeos.multipoints(pygeos.from_shapely(bldg_data["geometry"]))
    else: 
        building_points = pygeos.centroid(pygeos.union_all(block))

    if street_linestrings is not None:
        if type(street_linestrings) == gpd.geodataframe.GeoDataFrame:    
            street_linestrings = from_shapely_srid(geometry = street_linestrings, srid = 3395)
        assert type(street_linestrings) == np.ndarray and any(pygeos.is_geometry(street_linestrings)), 'street_linestrings is an invalid geometry type.'
        street_linestrings = street_linestrings[np.isin(pygeos.get_type_id(street_linestrings), [1,2,5])]
        if np.unique(pygeos.get_srid(street_linestrings))[0] == 4326:
            street_linestrings = transform_crs(geometry = street_linestrings, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395")
        assert len(np.unique(pygeos.get_srid(street_linestrings))) == 1 and np.unique(pygeos.get_srid(street_linestrings))[0] == 3395, 'street_linestrings is not epsg:4326 or epsg:3395.' 
        street_multilinestrings = pygeos.multilinestrings(street_linestrings)
        street_multilinestrings = pygeos.intersection(street_multilinestrings, block)

        if pygeos.intersects(block, street_multilinestrings): 
            nearest_external_street = 0
        else: 
            nearest_external_street = pygeos.distance(pygeos.centroid(building_points), pygeos.extract_unique_points(pygeos.multilinestrings(street_linestrings)))
        
        if not pygeos.is_empty(street_multilinestrings[0]):
            street_multilinestrings = pygeos.line_merge(street_multilinestrings)
            internal_buffer = pygeos.buffer(street_multilinestrings, radius=buffer_radius/2)
            external_buffer = pygeos.buffer(pygeos.get_exterior_ring(block), radius=buffer_radius)
            complete_buffer = pygeos.get_parts(pygeos.union(internal_buffer, external_buffer))
            exterior_access = street_linestrings[pygeos.intersects(street_linestrings, external_buffer)]
            complete_buffer = complete_buffer[pygeos.intersects(complete_buffer, pygeos.multilinestrings(exterior_access))]
            complete_buffer = pygeos.union_all(complete_buffer)
            off_network_geometry = pygeos.difference(street_multilinestrings, complete_buffer)
            off_network_length = pygeos.length(off_network_geometry)[0]
            on_network_geometry = pygeos.intersection(street_multilinestrings, complete_buffer)
            on_network_length = pygeos.length(on_network_geometry)[0]
            on_network_buffer = pygeos.buffer(on_network_geometry, radius = 1)
            if on_network_length > 0: 
                is_connected = True
        else: 
            on_network_length = 0
            off_network_length = 0
    else:
        nearest_external_street = np.nan
        on_network_length = np.nan
        off_network_length = np.nan

    if bldg_count not in [1,0]:
        voronoi = pygeos.get_parts(pygeos.voronoi_polygons(geometry=building_points, extend_to=block))
        if not pygeos.is_valid(block): block = pygeos.make_valid(block)
        block_parcels = pygeos.intersection(block, voronoi)    
        bldg_count = np.sum(pygeos.get_num_geometries(block_parcels))

        if is_connected is True:
            block_intersect = pygeos.intersects(on_network_buffer, block_parcels) 
            block_parcels_outer = block_parcels[block_intersect]
            block_layers.append(str(np.sum(pygeos.get_num_geometries(block_parcels_outer))))
            block_depth = len(block_layers)
            block_parcels = block_parcels[~block_intersect]
            if len(block_parcels_outer) == 0: 
                block_parcels_outer = block_parcels.copy()
        else: 
            block_parcels_outer = block_parcels.copy()

        parcel_residual = np.sum(pygeos.get_num_geometries(block_parcels))

        while parcel_residual > 0:
            try: block_reduced = pygeos.coverage_union_all(block_parcels_outer)
            except: block_reduced = pygeos.union_all(pygeos.make_valid(block_parcels_outer))
            block_parcels_inner = block_parcels[~pygeos.touches(block_parcels, block_reduced)]      
            block_parcels_outer = block_parcels[pygeos.touches(block_parcels, block_reduced)]  
            if np.sum(pygeos.get_num_geometries(block_parcels_outer)) == 0 and np.sum(pygeos.get_num_geometries(block_parcels_inner)) > 0:
                try: block_reduced = pygeos.coverage_union_all(block_parcels_inner)
                except: block_reduced = pygeos.union_all(pygeos.make_valid(block_parcels_inner))
                center = pygeos.get_coordinates(pygeos.centroid(block_reduced)).tolist()
                block_interior = pygeos.apply(block_reduced, lambda x: ((x - center)*.9999 + center) )
                block_exterior = pygeos.difference(block_reduced, block_interior)
                block_ring = pygeos.intersects(block_exterior, block_parcels_inner)
                block_parcels_outer = block_parcels_inner[block_ring]
                block_parcels_inner = block_parcels_inner[~block_ring]
            block_parcels = block_parcels_inner.copy()
            parcel_residual = np.sum(pygeos.get_num_geometries(block_parcels_inner))
            if parcel_residual == 0: 
                block_layers.append( str(np.sum(pygeos.get_num_geometries(block_parcels_outer))))
                block_depth = len(block_layers)
                break
            block_layers.append( str(np.sum(pygeos.get_num_geometries(block_parcels_outer))))
            block_depth = len(block_layers) 
        else:
            block_depth = len(block_layers)
    else:
        block_layers.append(str(bldg_count))
        block_depth = len(block_layers)

    if include_geometry is True:
        data = gpd.GeoDataFrame.from_dict({'block_id': [block_id], 'gadm_code': gadm_code, 'country_code': country_code, 'block_area': float(block_area), 'on_network_street_length': float(on_network_length), 'off_network_street_length': float(off_network_length), 'nearest_external_street': float(nearest_external_street), 'building_area': float(bldg_area), 'building_count': int(bldg_count), 'building_layers': ','.join(block_layers), 'k_complexity': int(block_depth), 'geometry': block_data['geometry']}).set_crs(epsg=3395)
        data = data.to_crs(4326)
    else:
        data = pd.DataFrame.from_dict({'block_id': [block_id], 'gadm_code': gadm_code, 'country_code': country_code, 'block_area': float(block_area), 'on_network_street_length': float(on_network_length), 'off_network_street_length': float(off_network_length), 'nearest_external_street': float(nearest_external_street), 'building_area': float(bldg_area), 'building_count': int(bldg_count), 'building_layers': ','.join(block_layers), 'k_complexity': int(block_depth)})
    return data


def compute_layers(block_id: str, block_col: str, block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame, street_linestrings: Union[pygeos.Geometry, gpd.GeoDataFrame] = None, buffer_radius: float = 100) -> gpd.GeoDataFrame:
    """
    Computes the k-complexity value for each layer of buildings within a block taking a GeoDataFrame of enclosed blocks, 
    the indexed building GeoDataFrame, and street linestrings in the form of a PyGEOS Geometry (or GeoDataFrame with slower performance). 
    Differs from compute_k() because provides a detailed rendering of the internal layers of building access.
    Args:
        block_id: string, unique identification code of block_id string in block_col
        block_col: string, column that contains the block_id codes (present in block_data and bldg_data)
        block_data: GeoDataFrame, output returned from build_blocks() function, requires CRS 4326 or 3395 for faster performance
        bldg_data: GeoDataFrame, output returned from index_buildings() function, requires CRS 4326 or 3395 for faster performance
        street_linestrings: GeoDataFrame or PyGEOS Geometry array (optional), linestrings representing street accesses, requires CRS 4326 or 3395
        buffer_radius: float (optional), extends disconnected street_linestrings to street network that are within a set tolerance in meters
    Returns:
        GeoDataFrame with block-layer geometries and the columns: 'block_id', 'gadm_code', 'country_code', 'block_property', 'building_count', 'k_complexity', 'geometry'
        Geometry projected in CRS 4326.        
    """
    assert block_data.crs in ['epsg:3395','epsg:4326'], "block_data is not epsg:4326 or epsg:3395."
    assert bldg_data.crs in ['epsg:3395','epsg:4326'], "bldg_data  is not epsg:4326 or epsg:3395."
    if block_data.crs == 'epsg:4326': block_data = block_data.to_crs(epsg=3395)
    if bldg_data.crs == 'epsg:4326': bldg_data = bldg_data.to_crs(epsg=3395)
    
    block_layers = [] 
    block_layers_geometry = []
    k_complexity = []
    
    bldg_data = bldg_data.loc[bldg_data[block_col] == block_id].copy() 
    block_data = block_data.loc[block_data[block_col] == block_id].copy() 
    block = pygeos.from_shapely(block_data['geometry'])
    
    country_code = block_data['country_code'].unique()[0] 
    gadm_code = block_data['gadm_code'].unique()[0] 
    bldg_array = pygeos.from_shapely(bldg_data['geometry'])
    bldg_count = np.sum(pygeos.get_num_geometries(bldg_array))
    is_connected = False
    off_network_length = 0 
    on_network_length = 0 

    if bldg_count > 0:
        gtypes = bldg_data['geometry'].geom_type.unique()
        if len(gtypes) > 1 or gtypes[0] != 'Point':
            bldg_data["geometry"] = bldg_data.centroid
        building_points = pygeos.multipoints(pygeos.from_shapely(bldg_data["geometry"]))
    else: 
        building_points = pygeos.centroid(pygeos.union_all(block))

    if street_linestrings is not None:
        if type(street_linestrings) == gpd.geodataframe.GeoDataFrame:    
            street_linestrings = from_shapely_srid(geometry = street_linestrings, srid = 3395)
        assert type(street_linestrings) == np.ndarray and any(pygeos.is_geometry(street_linestrings)), 'street_linestrings is an invalid geometry type.'
        street_linestrings = street_linestrings[np.isin(pygeos.get_type_id(street_linestrings), [1,2,5])]
        if np.unique(pygeos.get_srid(street_linestrings))[0] == 4326:
            street_linestrings = transform_crs(geometry = street_linestrings, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395")
        assert len(np.unique(pygeos.get_srid(street_linestrings))) == 1 and np.unique(pygeos.get_srid(street_linestrings))[0] == 3395, 'street_linestrings is not epsg:4326 or epsg:3395.'     
        street_multilinestrings = pygeos.multilinestrings(street_linestrings)
        street_multilinestrings = pygeos.intersection(street_multilinestrings, block)


        if not pygeos.is_empty(street_multilinestrings[0]):
            street_multilinestrings = pygeos.line_merge(street_multilinestrings)
            internal_buffer = pygeos.buffer(street_multilinestrings, radius=(buffer_radius/2))
            external_buffer = pygeos.buffer(pygeos.get_exterior_ring(block), radius=buffer_radius)
            complete_buffer = pygeos.get_parts(pygeos.union(internal_buffer, external_buffer))
            exterior_access = street_linestrings[pygeos.intersects(street_linestrings, external_buffer)]
            complete_buffer = complete_buffer[pygeos.intersects(complete_buffer, pygeos.multilinestrings(exterior_access))]
            complete_buffer = pygeos.union_all(complete_buffer)
            off_network_geometry = pygeos.difference(street_multilinestrings,complete_buffer)
            on_network_geometry = pygeos.intersection(street_multilinestrings,complete_buffer)
            on_network_buffer = pygeos.buffer(on_network_geometry, radius = 1)
            off_network_length = pygeos.length(off_network_geometry)[0]
            on_network_length = pygeos.length(on_network_geometry)[0]
            if on_network_length > 0: 
                is_connected = True
    
    if bldg_count not in [1,0]:
        voronoi = pygeos.get_parts(pygeos.voronoi_polygons(geometry=building_points, extend_to=block))
        if not pygeos.is_valid(block): block = pygeos.make_valid(block)
        block_parcels = pygeos.intersection(block, voronoi)    
        bldg_count = np.sum(pygeos.get_num_geometries(block_parcels))

        if is_connected is True:
            block_intersect = pygeos.intersects(on_network_buffer, block_parcels)
            block_parcels_outer = block_parcels[block_intersect]
            block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
            block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
            block_depth = len(block_layers)
            k_complexity.append(block_depth)
            block_parcels = block_parcels[~block_intersect]
            if len(block_parcels_outer) == 0: 
                block_parcels_outer = block_parcels.copy()
        else: 
            block_parcels_outer = block_parcels.copy()
 
        parcel_residual = np.sum(pygeos.get_num_geometries(block_parcels))
        
        while parcel_residual > 0:
            try: block_reduced = pygeos.coverage_union_all(block_parcels_outer)
            except: block_reduced = pygeos.union_all(pygeos.make_valid(block_parcels_outer))
            block_parcels_inner = block_parcels[~pygeos.touches(block_parcels,block_reduced)]      
            block_parcels_outer = block_parcels[pygeos.touches(block_parcels,block_reduced)] 
            if np.sum(pygeos.get_num_geometries(block_parcels_outer)) == 0 and np.sum(pygeos.get_num_geometries(block_parcels_inner)) > 0:
                try: block_reduced = pygeos.coverage_union_all(block_parcels_inner)
                except: block_reduced = pygeos.coverage_union_all(pygeos.make_valid(block_parcels_inner))
                center = pygeos.get_coordinates(pygeos.centroid(block_reduced)).tolist()
                block_interior = pygeos.apply(block_reduced, lambda x: ((x - center)*.9999 + center) )
                block_exterior = pygeos.difference(block_reduced, block_interior)
                block_ring = pygeos.intersects(block_exterior, block_parcels_inner)
                block_parcels_outer = block_parcels_inner[block_ring]
                block_parcels_inner = block_parcels_inner[~block_ring]
            block_parcels = block_parcels_inner.copy()
            parcel_residual = np.sum(pygeos.get_num_geometries(block_parcels_inner))
            if parcel_residual == 0: 
                block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
                block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
                block_depth = len(block_layers)
                k_complexity.append(block_depth)
                break
            block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
            block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
            k_complexity.append(len(block_layers))
        #else:
            #block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
            #block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
            #k_complexity.append(len(block_layers))
    else:
        block_layers.append(bldg_count)
        block_layers_geometry.append(pygeos.union_all(block))
        k_complexity.append(len(block_layers))

    data = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'gadm_code': gadm_code, 'country_code': country_code, 'block_property': 'building-parcels', 'building_count': block_layers, 'k_complexity': k_complexity, 'geometry': pygeos.to_shapely(block_layers_geometry)}).set_crs(epsg=3395)
    if on_network_length > 0:
        lines = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'gadm_code': gadm_code, 'country_code': country_code, 'block_property': 'on-network-streets', 'building_count': np.nan, 'k_complexity': np.nan, 'geometry': pygeos.to_shapely(on_network_geometry)}).set_crs(epsg=3395)
        data = pd.concat([data, lines])
    if off_network_length > 0:
        lines = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'gadm_code': gadm_code, 'country_code': country_code, 'block_property': 'off-network-streets', 'building_count': np.nan, 'k_complexity': np.nan, 'geometry': pygeos.to_shapely(off_network_geometry)}).set_crs(epsg=3395)
        data = pd.concat([data, lines])
    data = data.to_crs(4326)
    return data


def from_shapely_srid(geometry: gpd.GeoDataFrame, srid: int) -> pygeos.Geometry:
    """
    Transform GeoDataFrame to a PyGEOS Geometry array and set spatial reference identifier (SRID).
    Args:
        geometry: GeoSeries or GeoDataFrame 
        srid: integer representing spatial reference identifier (i.e., 4326)
    Returns:
        PyGEOS Geometry array
    """
    data = pygeos.set_srid(pygeos.from_shapely(geometry['geometry'].to_crs(epsg=srid)), srid)
    return data

def weighted_qcut(values, weights, q, **kwargs):
    """
    Return weighted quantile cuts from a given series, values.
    """
    if is_integer(q):
        quantiles = np.linspace(0, 1, q + 1)
    else:
        quantiles = q
    order = weights.iloc[values.argsort()].cumsum()
    bins = pd.cut(order / order.iloc[-1], quantiles, **kwargs)
    return bins.sort_index()

def mem_profile() -> str: 
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use

country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}    
    

def main(log_file: Path, country_chunk: list, chunk_size: int, core_count: int, blocks_dir: Path, streets_dir: Path, buildings_dir: Path, population_dir: Path, output_dir: Path):

    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info('---------')
    logging.info(f"Countries to process: {country_chunk}")

    buildings_per_partition=int(chunk_size)
    number_of_cores=int(core_count)

    # Make directory
    complexity_dir =  str(output_dir) + '/complexity'
    Path(complexity_dir).mkdir(parents=True, exist_ok=True)
    dask_dir =  str(output_dir) + '/complexity' + '/dask'
    Path(dask_dir).mkdir(parents=True, exist_ok=True)
    logging.info(f"complexity_dir: {complexity_dir}")
    
    # Check if country is completed
    output_file_list = list(filter(re.compile("complexity_").match, sorted(list(os.listdir(Path(complexity_dir))))))
    output_country_list = [(re.sub('complexity_', '', re.sub('.parquet', '', i))) for i in output_file_list] 
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
        print(country_code)
        logging.info(f"Processing: {country_code}")
        t0_country = time.time()

        # Build GADM list
        country_blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)    
        gadm_list = list(country_blocks['gadm_code'].unique())
        del country_blocks

        country_buildings = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_{country_code}.parquet')
        dask_folder_exists = os.path.isdir(Path(dask_dir) / f'{country_code}.parquet')
        if dask_folder_exists: 
            if len(os.listdir(Path(dask_dir) / f'{country_code}.parquet')) >= 1:            
                completed_blocks = dask.dataframe.read_parquet(path = Path(dask_dir) / f'{country_code}.parquet').compute()
                completed_gadm_list = list(completed_blocks['gadm_code'].unique())
                logging.info(f"Completed GADMs: {completed_gadm_list}")
                country_buildings = country_buildings[~country_buildings['gadm_code'].isin(completed_gadm_list)]
                #gadm_list = [x for x in gadm_list if x in set(country_buildings['gadm_code'].unique())] 
                del completed_blocks
        logging.info(f"GADMs to process: {gadm_list}")

        # Reconcile building and block GADM lists 
        building_gadm_list = country_buildings['gadm_code'].unique()
        building_check = [x for x in building_gadm_list if x not in set(gadm_list)] 
        logging.info(f"GADMs in building file not in blocks: {building_check}")
        block_check = [x for x in gadm_list if x not in set(building_gadm_list)] 
        logging.info(f"GADMs in blocks file not in buildings: {block_check}")
        intersected_list = [building_gadm_list, gadm_list] 
        gadm_list = list(set.intersection(*map(set,intersected_list)))
        # logging.info(f"GADM list intersected: {gadm_list}")
        # Filter the GADMs in buildings file to list that intersects block file
        country_buildings = country_buildings[country_buildings['gadm_code'].isin(gadm_list)]

        # Partition buildings into evenly sized chunks
        buildings_count = country_buildings.groupby(['gadm_code']).size().reset_index(name='counts')
        cuts = math.ceil(buildings_count.agg({'counts': 'sum'}).reset_index()[0][0]/buildings_per_partition)
        buildings_count['partition_index'] = weighted_qcut(buildings_count['gadm_code'], buildings_count['counts'], cuts, labels=False)
        gadm_dict = buildings_count[['gadm_code','partition_index']].drop_duplicates().reset_index(drop = True).groupby('partition_index')['gadm_code'].apply(list).to_dict()
        del country_buildings, cuts, buildings_count

        # Read in streets
        streets = gpd.read_parquet(path = Path(streets_dir) / f'streets_{country_code}.parquet', memory_map = True).to_crs(3395)
        streets = streets[streets['geometry'].notnull()].to_crs(3395)
        streets = from_shapely_srid(geometry = streets, srid = 3395)
        
        # Loop over the chunks
        for i in gadm_dict.items():

            gadm_chunk = i[1]
            logging.info(f"Node status: {mem_profile()}")
            logging.info(f"Processing chunk: {gadm_chunk}")

            try:
                blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True, filters = [('gadm_code', 'in', gadm_chunk)]).to_crs(3395)
            except Warning:
                print(f"Warning: No block data available for: {gadm_chunk}")
            try:
                street_zone = pygeos.buffer(pygeos.convex_hull(pygeos.union_all(pygeos.envelope(pygeos.from_shapely(blocks['geometry'])))), radius = 1000)
                street_network = streets[pygeos.intersects(street_zone, streets)]
                if len(street_network) == 0:
                    street_network = streets.copy()
            except Warning:
                print(f"Warning: No street data available for: {gadm_chunk}")
            try:
                buildings = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_{country_code}.parquet', memory_map = True, filters = [('gadm_code', 'in', gadm_chunk)]).to_crs(3395)
            except Warning:
                print(f"Warning: No building data available for: {gadm_chunk}")
                
            # Spatial join
            buildings = index_building_points(gadm_block_data = blocks, bldg_data = buildings)
            buildings = buildings.to_crs(3395)
            block_list = list(blocks['block_id'].unique())
            
            # Parallelize block computation
            k_gadm = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'block_area': pd.Series(dtype='float'), 'on_network_street_length': pd.Series(dtype='float'),  'off_network_street_length': pd.Series(dtype='float'), 'nearest_external_street': pd.Series(dtype='float'), 'building_area': pd.Series(dtype='float'), 'building_count': pd.Series(dtype='int'), 'building_layers': pd.Series(dtype='object'), 'k_complexity': pd.Series(dtype='int')})    
            pool = multiprocessing.Pool(processes= number_of_cores) 
            output = pool.map(functools.partial(compute_k, block_col = 'block_id', block_data = blocks, bldg_data = buildings, street_linestrings = street_network, buffer_radius = 60, include_geometry = False), block_list)
            pool.close() 
            pool.join()
            for i,j in enumerate(output): k_gadm = pd.concat([k_gadm, output[i]], ignore_index=True)

            # Incremental file build
            k_gadm = dask.dataframe.from_pandas(data = k_gadm, npartitions = 1) 
            dask.dataframe.to_parquet(df = k_gadm, path = Path(dask_dir) / f'{country_code}.parquet', engine='pyarrow', compression='snappy', append=True, ignore_divisions=True)

            ## Parallelize layer computation
            #k_layers = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'block_property': pd.Series(dtype='str'), 'building_count': pd.Series(dtype='int'), 'k_complexity': pd.Series(dtype='int'), 'geometry': gpd.GeoSeries()})
            #pool = multiprocessing.Pool(processes= number_of_cores) 
            #output = pool.map(functools.partial(compute_layers, block_col = 'block_id', block_data = blocks, bldg_data = buildings, street_linestrings = street_network, buffer_radius = 60), block_list)
            #pool.close() 
            #pool.join()
            #for i,j in enumerate(output): k_layers = pd.concat([k_layers, output[i]], ignore_index=True)

        # Combine partitioned dataframes into single parquet
        k_bulk = dask.dataframe.read_parquet(path = Path(dask_dir) / f'{country_code}.parquet').compute()
        k_bulk.to_parquet(path = Path(complexity_dir) / f'complexity_{country_code}.parquet')

        #k_layers.to_parquet(path = Path(complexity_dir) / f'complexity_{country_code}_layers.parquet')

        t1_country = time.time()
        logging.info(f"Finished: {country_code}, {str(round((t1_country-t0_country)/60,3))} minutes")
        logging.info('---------')
        logging.info('---------')


def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute K.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to log file") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="Array of country codes in ISO 3166-1 alpha-3 format")
    parser.add_argument('--chunk_size', required=False, type=int, dest="chunk_size", help="Number of building points per chunk")
    parser.add_argument('--core_count', required=False, type=int, dest="core_count", help="Number of core processes to run in parallel")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory")
    parser.add_argument('--streets_dir', required=True, type=Path, dest="streets_dir", help="Path to streets directory")
    parser.add_argument('--buildings_dir', required=True, type=Path, dest="buildings_dir", help="Path to buildings directory")
    parser.add_argument('--population_dir', required=True, type=Path, dest="population_dir", help="Path to population directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to outputs directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))


