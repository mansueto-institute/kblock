
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

import pyproj
from typing import List, Union
from pandas._libs.lib import is_integer
from pathlib import Path
import psutil
import gc
import warnings
import logging
import time
import re
import os
import functools
import argparse
import math

import multiprocessing
import dask 
import dask.dataframe
import momepy

from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')
#warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

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

def compute_k(block_id: str, block_col: str, block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame, street_linestrings: Union[pygeos.Geometry, gpd.GeoDataFrame] = None, buffer_radius: float = 100, include_geometry: bool = True) -> Union[gpd.GeoDataFrame, pd.DataFrame]:
    """
    Computes the k-complexity value for each block as well as block level metadata taking a GeoDataFrame of enclosed blocks, 
    the indexed building GeoDataFrame, and street linestrings in the form of a PyGEOS Geometry (or GeoDataFrame with slower performance).
    Args:
        block_id: string, unique identification code in block_col
        block_col: string, column name that contains the block_id codes (present in block_data and bldg_data)
        block_data: GeoDataFrame output returned from build_blocks() function, requires CRS 4326 or 3395
        bldg_data: GeoDataFrame, requires CRS 4326 or 3395, must contain block_col containing block_id
        street_linestrings: GeoDataFrame or PyGEOS Geometry array (optional), linestrings representing street networks, requires CRS 4326 or 3395
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
    assert block_data.shape[0] == 1, 'block_id is not unique.'
    block = pygeos.from_shapely(block_data['geometry'])
    if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)

    if pygeos.get_num_geometries(block)[0] > 1:
        block = pygeos.get_parts(block)
        block = np.array([x for x in block if (pygeos.get_type_id(x) == 3 or pygeos.get_type_id(x) == 6)])
        block = np.array([pygeos.buffer(pygeos.multipolygons(pygeos.get_parts(block)),1)])
        #block = pygeos.union_all(np.array([pygeos.buffer(pygeos.multipolygons(pygeos.get_parts(block)),1)]))
        
    country_code = block_data['country_code'].unique()[0] 
    gadm_code = block_data['gadm_code'].unique()[0] 
    block_area = round(pygeos.area(pygeos.multipolygons(pygeos.from_shapely(block_data['geometry']))),2)
    bldg_array = pygeos.from_shapely(bldg_data['geometry'])
    if not pygeos.is_valid(bldg_array).all(): bldg_array = pygeos.make_valid(bldg_array)
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
    if not pygeos.is_valid(building_points).all(): building_points = pygeos.make_valid(building_points)

    if street_linestrings is not None:
        if type(street_linestrings) == gpd.geodataframe.GeoDataFrame:    
            street_linestrings = from_shapely_srid(geometry = street_linestrings, srid = 3395)
        assert type(street_linestrings) == np.ndarray and any(pygeos.is_geometry(street_linestrings)), 'street_linestrings is an invalid geometry type.'
        street_linestrings = street_linestrings[np.isin(pygeos.get_type_id(street_linestrings), [1,2,5])]
        if np.unique(pygeos.get_srid(street_linestrings))[0] == 4326:
            street_linestrings = transform_crs(geometry = street_linestrings, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395")
        assert len(np.unique(pygeos.get_srid(street_linestrings))) == 1 and np.unique(pygeos.get_srid(street_linestrings))[0] == 3395, 'street_linestrings is not epsg:4326 or epsg:3395.' 
        street_multilinestrings = pygeos.multilinestrings(street_linestrings)
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
        street_multilinestrings = pygeos.intersection(street_multilinestrings, block)
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)

        if pygeos.intersects(block, street_multilinestrings): 
            nearest_external_street = 0
        else: 
            nearest_external_street = pygeos.distance(pygeos.centroid(building_points), pygeos.extract_unique_points(pygeos.multilinestrings(street_linestrings)))

        assert len(street_multilinestrings) == 1, 'street_linestrings not formatted as single geometry.'
        if not pygeos.is_empty(street_multilinestrings[0]):
            street_multilinestrings = pygeos.line_merge(street_multilinestrings)
            if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
            internal_buffer = pygeos.buffer(street_multilinestrings, radius=buffer_radius/2)
            external_buffer = pygeos.buffer(pygeos.get_exterior_ring(block), radius=buffer_radius)
            if not pygeos.is_valid(internal_buffer).all(): internal_buffer= pygeos.make_valid(internal_buffer)
            if not pygeos.is_valid(external_buffer).all(): external_buffer= pygeos.make_valid(external_buffer)
            complete_buffer = pygeos.get_parts(pygeos.union(internal_buffer, external_buffer))
            if not pygeos.is_valid(street_linestrings).all(): street_linestrings = pygeos.make_valid(street_linestrings)
            exterior_access = street_linestrings[pygeos.intersects(street_linestrings, external_buffer)]
            complete_buffer = complete_buffer[pygeos.intersects(complete_buffer, pygeos.multilinestrings(exterior_access))]
            complete_buffer = pygeos.union_all(complete_buffer)
            if not pygeos.is_valid(complete_buffer).all(): complete_buffer = pygeos.make_valid(complete_buffer)
            off_network_geometry = pygeos.difference(street_multilinestrings, complete_buffer)
            # assert len(off_network_geometry) == 1, 'off_network_geometry more than .'
            off_network_length = pygeos.length(off_network_geometry)[0]
            on_network_geometry = pygeos.intersection(street_multilinestrings, complete_buffer)
            # assert len(on_network_geometry) == 1, 'on_network_geometry unable to format as MultiLineString.'
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
        if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)
        if not pygeos.is_valid(voronoi).all(): voronoi = pygeos.make_valid(voronoi)
        block_parcels = pygeos.intersection(block, voronoi)  
        if not pygeos.is_valid(block_parcels).all(): block_parcels = pygeos.make_valid(block_parcels)  
        bldg_count = np.sum(pygeos.get_num_geometries(block_parcels))

        if is_connected is True:
            block_intersect = pygeos.intersects(on_network_buffer, block_parcels) 
            block_parcels_outer = block_parcels[block_intersect]
            block_layers.append(str(np.sum(pygeos.get_num_geometries(block_parcels_outer))))
            block_depth = len(block_layers)
            block_parcels = block_parcels[~block_intersect]
            if not pygeos.is_valid(block_parcels).all(): block_parcels = pygeos.make_valid(block_parcels)
            if len(block_parcels_outer) == 0: 
                block_parcels_outer = block_parcels.copy()
        else: 
            block_parcels_outer = block_parcels.copy()

        parcel_residual = np.sum(pygeos.get_num_geometries(block_parcels))

        while parcel_residual > 0:
            try: block_reduced = pygeos.coverage_union_all(block_parcels_outer)
            except: block_reduced = pygeos.union_all(pygeos.make_valid(block_parcels_outer))
            if not pygeos.is_valid(block_parcels).all(): block_parcels = pygeos.make_valid(block_parcels)
            if not pygeos.is_valid(block_reduced).all(): block_reduced = pygeos.make_valid(block_reduced)
            block_parcels_inner = block_parcels[~pygeos.touches(block_parcels, block_reduced)]      
            if not pygeos.is_valid(block_parcels_inner).all(): block_parcels_inner = pygeos.make_valid(block_parcels_inner)
            block_parcels_outer = block_parcels[pygeos.touches(block_parcels, block_reduced)]  
            if not pygeos.is_valid(block_parcels_outer).all(): block_parcels_outer = pygeos.make_valid(block_parcels_outer)
            if np.sum(pygeos.get_num_geometries(block_parcels_outer)) == 0 and np.sum(pygeos.get_num_geometries(block_parcels_inner)) > 0:
                try: block_reduced = pygeos.coverage_union_all(block_parcels_inner)
                except: block_reduced = pygeos.union_all(pygeos.make_valid(block_parcels_inner))
                center = pygeos.get_coordinates(pygeos.centroid(block_reduced)).tolist()
                block_interior = pygeos.apply(block_reduced, lambda x: ((x - center)*.9999 + center) )
                if not pygeos.is_valid(block_reduced).all(): block_reduced = pygeos.make_valid(block_reduced) 
                if not pygeos.is_valid(block_interior).all(): block_interior = pygeos.make_valid(block_interior) 
                block_exterior = pygeos.difference(block_reduced, block_interior)
                if not pygeos.is_valid(block_exterior).all(): block_exterior = pygeos.make_valid(block_exterior) 
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


def compute_layers(block_id: str, block_col: str, block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame, street_linestrings: Union[pygeos.Geometry, gpd.GeoDataFrame] = None, use_building_polygons: bool = True, buffer_radius: float = 100, shrink_value: float = 0.1, segment_value: float = 0.5, threshold_value: float = 0.01) -> gpd.GeoDataFrame:
    """
    Computes the k-complexity value for each layer of buildings within a block taking a GeoDataFrame of enclosed blocks, 
    the indexed building GeoDataFrame, and street linestrings in the form of a PyGEOS Geometry (or GeoDataFrame with slower performance). 
    Differs from compute_k() because provides a detailed rendering of the internal layers of building access.
    Args:
        block_id: string, unique identification code of block_id string in block_col
        block_col: string, column that contains the block_id codes (present in block_data and bldg_data)
        block_data: GeoDataFrame, output returned from build_blocks() function, requires CRS 4326 or 3395 
        bldg_data: GeoDataFrame, GeoDataFrame, requires CRS 4326 or 3395, must contain block_col containing block_id
        street_linestrings: GeoDataFrame or PyGEOS Geometry array (optional), linestrings representing street accesses, requires CRS 4326 or 3395
        use_building_polygons: bool (optional), contour parcels around input polygon bldg_data, otherwise defaults to points
        buffer_radius: float (optional), extends disconnected street_linestrings to street network that are within a set tolerance in meters
        For info about shrink_value, segment_value, threshold_value see: http://docs.momepy.org/en/stable/generated/momepy.Tessellation.html
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
    assert block_data.shape[0] == 1, 'block_id is not unique.'
    block = pygeos.from_shapely(block_data['geometry'])
    if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)
    
    if pygeos.get_num_geometries(block)[0] > 1:
        block = pygeos.get_parts(block)
        block = np.array([x for x in block if (pygeos.get_type_id(x) == 3 or pygeos.get_type_id(x) == 6)])
        block = np.array([pygeos.buffer(pygeos.multipolygons(pygeos.get_parts(block)),1)])

    country_code = block_data['country_code'].unique()[0] 
    gadm_code = block_data['gadm_code'].unique()[0] 
    bldg_array = pygeos.from_shapely(bldg_data['geometry'])
    if not pygeos.is_valid(bldg_array).all(): bldg_array = pygeos.make_valid(bldg_array)
    bldg_count = np.sum(pygeos.get_num_geometries(bldg_array))
    is_connected = False
    off_network_length = 0 
    on_network_length = 0 

    if use_building_polygons == False:
        if bldg_count > 0:
            gtypes = bldg_data['geometry'].geom_type.unique()
            if len(gtypes) > 1 or gtypes[0] != 'Point':
                bldg_data["geometry"] = bldg_data.centroid
            building_points = pygeos.multipoints(pygeos.from_shapely(bldg_data["geometry"]))
        else: 
            building_points = pygeos.centroid(pygeos.union_all(block))
        if not pygeos.is_valid(building_points).all(): building_points = pygeos.make_valid(building_points)
    else:
        if bldg_count > 0:
            building_polygons = bldg_data.copy()
        else: 
            building_polygons = pygeos.centroid(pygeos.union_all(block))
        if not pygeos.is_valid(building_polygons).all(): building_polygons = pygeos.make_valid(building_polygons)

    if street_linestrings is not None:
        if type(street_linestrings) == gpd.geodataframe.GeoDataFrame:    
            street_linestrings = from_shapely_srid(geometry = street_linestrings, srid = 3395)
        assert type(street_linestrings) == np.ndarray and any(pygeos.is_geometry(street_linestrings)), 'street_linestrings is an invalid geometry type.'
        street_linestrings = street_linestrings[np.isin(pygeos.get_type_id(street_linestrings), [1,2,5])]
        if np.unique(pygeos.get_srid(street_linestrings))[0] == 4326:
            street_linestrings = transform_crs(geometry = street_linestrings, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395")
        assert len(np.unique(pygeos.get_srid(street_linestrings))) == 1 and np.unique(pygeos.get_srid(street_linestrings))[0] == 3395, 'street_linestrings is not epsg:4326 or epsg:3395.'     
        street_multilinestrings = pygeos.multilinestrings(street_linestrings)
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
        street_multilinestrings = pygeos.intersection(street_multilinestrings, block)
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)

        assert len(street_multilinestrings) == 1, 'street_linestrings not formatted as single geometry.'
        if not pygeos.is_empty(street_multilinestrings[0]):
            street_multilinestrings = pygeos.line_merge(street_multilinestrings)
            if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
            internal_buffer = pygeos.buffer(street_multilinestrings, radius=(buffer_radius/2))
            external_buffer = pygeos.buffer(pygeos.get_exterior_ring(block), radius=buffer_radius)
            if not pygeos.is_valid(internal_buffer).all(): internal_buffer= pygeos.make_valid(internal_buffer)
            if not pygeos.is_valid(external_buffer).all(): external_buffer= pygeos.make_valid(external_buffer)
            complete_buffer = pygeos.get_parts(pygeos.union(internal_buffer, external_buffer))
            if not pygeos.is_valid(street_linestrings).all(): street_linestrings = pygeos.make_valid(street_linestrings)
            exterior_access = street_linestrings[pygeos.intersects(street_linestrings, external_buffer)]
            complete_buffer = complete_buffer[pygeos.intersects(complete_buffer, pygeos.multilinestrings(exterior_access))]
            complete_buffer = pygeos.union_all(complete_buffer)
            if not pygeos.is_valid(complete_buffer).all(): complete_buffer = pygeos.make_valid(complete_buffer)
            off_network_geometry = pygeos.difference(street_multilinestrings,complete_buffer)            
            on_network_geometry = pygeos.intersection(street_multilinestrings,complete_buffer)
            on_network_buffer = pygeos.buffer(on_network_geometry, radius = 1)
            off_network_length = pygeos.length(off_network_geometry)[0]
            on_network_length = pygeos.length(on_network_geometry)[0]
            if on_network_length > 0: 
                is_connected = True
    
    if bldg_count not in [1,0]:

        if use_building_polygons == True:
            building_polygons['bldg_id'] = momepy.unique_id(building_polygons)
            enclosed_tess = momepy.Tessellation(gdf=building_polygons, unique_id='bldg_id', enclosures=block_data, enclosure_id='block_id', shrink = shrink_value, segment = segment_value, threshold = threshold_value).tessellation
            voronoi = pygeos.get_parts(pygeos.from_shapely(enclosed_tess['geometry']))
        else:
            voronoi = pygeos.get_parts(pygeos.voronoi_polygons(geometry=building_points, extend_to=block))

        if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)
        if not pygeos.is_valid(voronoi).all(): voronoi = pygeos.make_valid(voronoi)
        block_parcels = pygeos.intersection(block, voronoi)    
        if not pygeos.is_valid(block_parcels).all(): block_parcels = pygeos.make_valid(block_parcels)  
        bldg_count = np.sum(pygeos.get_num_geometries(block_parcels))

        if is_connected is True:
            block_intersect = pygeos.intersects(on_network_buffer, block_parcels)
            block_parcels_outer = block_parcels[block_intersect]
            block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
            block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
            block_depth = len(block_layers)
            k_complexity.append(block_depth)
            block_parcels = block_parcels[~block_intersect]
            if not pygeos.is_valid(block_parcels).all(): block_parcels = pygeos.make_valid(block_parcels)
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
                except: block_reduced = pygeos.union_all(pygeos.make_valid(block_parcels_inner))
                if not pygeos.is_valid(block_reduced).all(): block_reduced = pygeos.make_valid(block_reduced)
                center = pygeos.get_coordinates(pygeos.centroid(block_reduced)).tolist()
                block_interior = pygeos.apply(block_reduced, lambda x: ((x - center)*.9999 + center) )
                block_exterior = pygeos.difference(block_reduced, block_interior)
                block_ring = pygeos.intersects(block_exterior, block_parcels_inner)
                block_parcels_outer = block_parcels_inner[block_ring]
                if not pygeos.is_valid(block_parcels_outer).all(): block_parcels_outer = pygeos.make_valid(block_parcels_outer)
                block_parcels_inner = block_parcels_inner[~block_ring]
                if not pygeos.is_valid(block_parcels_inner).all(): block_parcels_inner = pygeos.make_valid(block_parcels_inner)
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

def mem_profile_detail() -> str:
    mem = psutil.virtual_memory()
    mem_report = f'''
        Percent in use: {mem.percent}
        Available: {round(mem.available / mem.total, 3)}
        Used: {round(mem.used / mem.total, 3)}
        Free: {round(mem.free / mem.total, 3)}
        Active: {round(mem.active / mem.total, 3)}
        Inactive: {round(mem.inactive / mem.total, 3)}'''
    return mem_report

def main(log_file: Path, country_chunk: list, chunk_size: int, core_count: int, blocks_dir: Path, streets_dir: Path, buildings_dir: Path, output_dir: Path):

    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info('---------')
    logging.info(f"Countries to process: {country_chunk}")

    buildings_per_partition=int(chunk_size)
    number_of_cores=int(core_count)
    assert buildings_per_partition > 0 and type(buildings_per_partition) is int,"chunk_size argument must be an integer and greater than 0."
    assert number_of_cores > 0 and type(number_of_cores) is int,"core_count argument must be an integer and greater than 0."
    if number_of_cores > multiprocessing.cpu_count():
        warnings.warn(f"core_count argument exceeds CPUs on machine. Available CPUs: {multiprocessing.cpu_count()}", UserWarning)

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
    
    # Check files and match to countries in each source directory
    blocks_dir_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    blocks_dir_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in blocks_dir_list] 

    buildings_dir_list = list(filter(re.compile("buildings_points_").match, sorted(list(os.listdir(Path(buildings_dir))))))
    buildings_dir_list = [(re.sub('buildings_points_', '', re.sub('.parquet', '', i))) for i in buildings_dir_list] 

    streets_dir_list = list(filter(re.compile("streets_").match, sorted(list(os.listdir(Path(streets_dir))))))
    streets_dir_list = [(re.sub('streets_', '', re.sub('.parquet', '', i))) for i in streets_dir_list] 

    # Check country_chunk parameter against the country files across directories
    if country_chunk: 
        country_list = [x for x in country_chunk if x not in set(output_country_list)] 
        if not country_list: 
            raise ValueError('All countries in country_chunk argument already finished.')

        intersected_list = [blocks_dir_list, buildings_dir_list, streets_dir_list] 
        intersected_list = list(set.intersection(*map(set,intersected_list)))
        chunk_missing = [x for x in country_list if x not in set(intersected_list)] 
        if chunk_missing:
            raise ValueError(f'Elements in country_chunk argument have missing input files: {chunk_missing}')
        
        country_list = [x for x in country_list if x in set(intersected_list)] 
        if not country_list: 
            raise ValueError('Empty country list.')
    else: 
        raise ValueError('Empty country_chunk argument.')
            
    logging.info(f"Remaining countries: {country_list}")

    # Process country list
    for country_code in country_list:
        print(country_code)
        logging.info(f"Processing: {country_code}")
        t0_country = time.time()

        # Build block_id list
        country_blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)    
        # gadm_list = list(country_blocks['gadm_code'].unique())
        block_list = list(country_blocks['block_id'].unique())
        del country_blocks

        # Read in bulilding points
        country_buildings = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_points_{country_code}.parquet')
        dask_folder_exists = os.path.isdir(Path(dask_dir) / f'{country_code}.parquet')
        if dask_folder_exists: 
            if len(os.listdir(Path(dask_dir) / f'{country_code}.parquet')) >= 1:            
                completed_blocks = dask.dataframe.read_parquet(path = Path(dask_dir) / f'{country_code}.parquet').compute()
                # completed_gadm_list = list(completed_blocks['gadm_code'].unique())
                # country_buildings = country_buildings[~country_buildings['gadm_code'].isin(completed_gadm_list)]
                completed_block_list = list(completed_blocks['block_id'].unique())
                country_buildings = country_buildings[~country_buildings['block_id'].isin(completed_block_list)]
                logging.info(f"Completed blocks: {len(completed_block_list)} {round(len(completed_block_list)/len(block_list),2)*100} %")
                del completed_blocks

        # Reconcile building and block GADM lists 
        building_block_list = country_buildings['block_id'].unique()
        building_check = [x for x in building_block_list if x not in set(block_list)] 
        logging.info(f"Blocks in building file not in blocks: {len(building_check)}")
        block_check = [x for x in block_list if x not in set(building_block_list)] 
        logging.info(f"Blocks in blocks file not in buildings: {len(block_check)}")
        intersected_list = [building_block_list, block_list] 
        block_list = list(set.intersection(*map(set,intersected_list)))
        logging.info(f"Blocks to process: {len(block_list)}")
        
        # Partition buildings into evenly sized chunks
        country_buildings = country_buildings[country_buildings['block_id'].isin(block_list)]
        country_buildings['building_geogroup'] = country_buildings['building_geohash'].str.slice(start=2, stop=6)
        buildings_count = country_buildings.groupby(['building_geogroup', 'block_id']).size().reset_index(name='counts')
        buildings_count = buildings_count.sort_values(by='building_geogroup', ascending=False).reset_index(drop=True)
        cuts = math.ceil(buildings_count.agg({'counts': 'sum'}).reset_index()[0][0]/buildings_per_partition)
        buildings_count['partition_index'] = weighted_qcut(buildings_count['block_id'], buildings_count['counts'], cuts, labels=False)
        chunk_dict = buildings_count[['block_id','partition_index']].drop_duplicates().reset_index(drop = True).groupby('partition_index')['block_id'].apply(list).to_dict()
        logging.info(f'Number of partitions: {len(chunk_dict)}')
        del country_buildings, cuts, buildings_count
        
        # Read in streets
        streets = gpd.read_parquet(path = Path(streets_dir) / f'streets_{country_code}.parquet', memory_map = True).to_crs(3395)
        streets = streets[streets['geometry'].notnull()].to_crs(3395)
        streets = from_shapely_srid(geometry = streets, srid = 3395)
        
        # Loop over the chunks
        for i in chunk_dict.items():

            chunk_list = i[1]
            
            logging.info(f"Node status: {mem_profile_detail()}")
            #logging.info(f"Uncollectable garbage: {gc.garbage}")
            #logging.info(f"Garbage stats: {gc.get_stats()}")
            logging.info(f"Blocks in chunk: {len(chunk_list)}")

            # Read in block data within chunk
            try:
                blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True, filters = [('block_id', 'in', chunk_list)]).to_crs(3395)
            except Warning:
                logging.info(f"Warning: No block data available for: {chunk_list}")
            logging.info(f"GADMs in chunk: {blocks['gadm_code'].unique()}")
            logging.info(f"Geohashes in chunk: {blocks['block_geohash'].str.slice(start=0, stop=4).unique()}")
            
            # Read in streets and limit to streets that intersect chunk
            try:
                street_zone = pygeos.buffer(pygeos.convex_hull(pygeos.union_all(pygeos.envelope(pygeos.from_shapely(blocks['geometry'])))), radius = 1000)
                street_network = streets[pygeos.intersects(street_zone, streets)]
                del street_zone
                if len(street_network) == 0:
                    street_network = streets.copy()
            except Warning:
                logging.info(f"Warning: No street data available for: {chunk_list}")

            # Read in buildings within chunk
            try:
                buildings = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_points_{country_code}.parquet', memory_map = True, filters = [('block_id', 'in', chunk_list)]).to_crs(3395)
            except Warning:
                logging.info(f"Warning: No building data available for: {chunk_list}")
            logging.info(f"Buildings in chunk: {buildings.shape[0]}")
                
            # Spatial join of buildings to blocks            
            buildings = buildings.to_crs(3395)
            parallel_block_list = list(blocks['block_id'].unique())
            
            k_output = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'block_area': pd.Series(dtype='float'), 'on_network_street_length': pd.Series(dtype='float'),  'off_network_street_length': pd.Series(dtype='float'), 'nearest_external_street': pd.Series(dtype='float'), 'building_area': pd.Series(dtype='float'), 'building_count': pd.Series(dtype='int'), 'building_layers': pd.Series(dtype='object'), 'k_complexity': pd.Series(dtype='int')})    

            # Parallelize block computation of k and related metrics
            if number_of_cores > 1: 
                pool = multiprocessing.Pool(processes = number_of_cores, maxtasksperchild = 100) 
                k_chunk = pool.map(functools.partial(compute_k, block_col = 'block_id', block_data = blocks, bldg_data = buildings, street_linestrings = street_network, buffer_radius = 60, include_geometry = False), parallel_block_list)
                pool.close() 
                pool.join()
                for i,j in enumerate(k_chunk): k_output = pd.concat([k_output, k_chunk[i]], ignore_index=True)
            else: 
                k_chunk = list(map(lambda x: compute_k(block_id = x, block_col = 'block_id', block_data = blocks, bldg_data = buildings, street_linestrings = street_network, buffer_radius = 60, include_geometry = False), parallel_block_list))
                for i,j in enumerate(k_chunk): k_output = pd.concat([k_output, k_chunk[i]], ignore_index=True)

            # Incremental file build
            k_output = dask.dataframe.from_pandas(data = k_output, npartitions = 1) 
            dask.dataframe.to_parquet(df = k_output, path = Path(dask_dir) / f'{country_code}.parquet', engine='pyarrow', compression='snappy', append=True, ignore_divisions=True)
            del blocks, street_network, buildings, k_output, k_chunk 
            gc.collect()

            ## Parallelize layer computation
            #k_layers = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'block_property': pd.Series(dtype='str'), 'building_count': pd.Series(dtype='int'), 'k_complexity': pd.Series(dtype='int'), 'geometry': gpd.GeoSeries()})
            #pool = multiprocessing.Pool(processes= number_of_cores) 
            #output = pool.map(functools.partial(compute_layers, block_col = 'block_id', block_data = blocks, bldg_data = buildings, street_linestrings = street_network, buffer_radius = 60), parallel_block_list)
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
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to log file.") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="Array of country codes in ISO 3166-1 alpha-3 format.")
    parser.add_argument('--chunk_size', required=False, type=int, dest="chunk_size", help="Number of building per chunk.")
    parser.add_argument('--core_count', required=False, type=int, dest="core_count", help="Number of core processes to run in parallel.")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory.")
    parser.add_argument('--streets_dir', required=True, type=Path, dest="streets_dir", help="Path to streets directory.")
    parser.add_argument('--buildings_dir', required=True, type=Path, dest="buildings_dir", help="Path to /buildings/points directory. It is recommended to use building centroids for faster computation.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to top level output directory. Creates /complexity subdirectory.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))




