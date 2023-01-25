
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

import pyproj
from pandas._libs.lib import is_integer
from pathlib import Path
import psutil
import gc
import warnings
import logging
import time
import re
import os
import collections
import argparse
import math
import typing

import multiprocessing
import dask 
import dask.dataframe
import momepy

from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')
#warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

def inputs_generator(block_id_col: str, blocks: gpd.GeoDataFrame, buildings: gpd.GeoDataFrame, streets: gpd.GeoDataFrame) -> typing.Generator[typing.Tuple[str, gpd.GeoSeries, gpd.GeoSeries, gpd.GeoSeries], None, None]:
    for block_id in blocks[block_id_col].unique():
        yield (
            block_id,
            blocks.loc[(blocks[block_id_col] == block_id), "geometry"],
            buildings.loc[(buildings[block_id_col] == block_id), "geometry"],
            streets.loc[(streets[block_id_col] == block_id), "geometry"],
        )    

def compute_k(block_id: str, block: gpd.GeoSeries, buildings: gpd.GeoSeries, streets: gpd.GeoSeries, buffer_radius: float = 100, include_geometry: bool = False, srid: int = 3395) -> typing.Union[gpd.GeoDataFrame, pd.DataFrame]:
    """
    Computes the k-complexity value for an enclosed block, buildings, and street linestrings.
    Args:
        block_id: string, unique identification code
        block: GeoSeries of a single Polygon representing enclosed block
        buildings: GeoSeries of Points or Polygons representing buildings
        streets: GeoSeries of Linestrings or MultiLinestrings representing street networks
        buffer_radius: float (optional), defaults to 100 meters, extends disconnected street_linestrings to street network that are within a set tolerance in meters
        include_geometry: bool (optional), defaults to True, if False the block geometry will not be included and function will return a Pandas DataFrame
        srid: integer (optional), defaults to 3395, representing spatial reference identifier, for meter conversion
    Returns:
        DataFrame or GeoDataFrame for individual block geometry. Geometry projected in CRS 4326.
    """
    block_layers = [] 

    # Check block geometry
    assert block.geom_type.unique().all() == 'Polygon', 'block geometry is not a Polygon.'
    block = pygeos.set_srid(pygeos.from_shapely(block.to_crs(epsg=srid)), srid)
    if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)
    if pygeos.get_num_geometries(block)[0] > 1:
        block = pygeos.get_parts(block)
        block = np.array([x for x in block if (pygeos.get_type_id(x) == 3 or pygeos.get_type_id(x) == 6)])
        block = np.array([pygeos.buffer(pygeos.multipolygons(pygeos.get_parts(block)),1)])
        if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)

    # Check buildings geometry
    if buildings.geom_type.unique().any() in ['MultiPoint','MultiPolygon']: buildings = buildings.explode(ignore_index=True)
    if buildings.geom_type.unique().any() in ['Polygon']: buildings = buildings.centroid
    assert buildings.geom_type.unique().all() in [None, 'Point'], 'building geometry is not a Point.'
    buildings = pygeos.set_srid(pygeos.from_shapely(buildings.to_crs(epsg=srid)), srid)
    if not pygeos.is_valid(buildings).all(): buildings = pygeos.make_valid(buildings)
    buildings = buildings[np.isin(pygeos.get_type_id(buildings), [0,4])]
    buildings = pygeos.multipoints(buildings)

    is_connected = False

    # Check street networks
    if streets.geom_type.unique().any() in ['MultiLineString', 'LineString']:
        streets = pygeos.set_srid(pygeos.from_shapely(streets.to_crs(epsg=srid)), srid)
        if not pygeos.is_valid(streets).all(): streets = pygeos.make_valid(streets)
        streets = streets[np.isin(pygeos.get_type_id(streets), [1,2,5])]

        street_multilinestrings = pygeos.multilinestrings(pygeos.get_parts(streets))
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
        street_multilinestrings = pygeos.intersection(street_multilinestrings, block)
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)

        if pygeos.intersects(block, street_multilinestrings).any(): 
            nearest_external_street = 0
        else: 
            nearest_external_street = pygeos.distance(pygeos.centroid(buildings), pygeos.extract_unique_points(pygeos.multilinestrings(pygeos.get_parts(streets))))

        if not pygeos.is_empty(street_multilinestrings).all():
            
            street_multilinestrings = pygeos.line_merge(street_multilinestrings)
            if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
            internal_buffer = pygeos.buffer(street_multilinestrings, radius=buffer_radius/2)
            if not pygeos.is_valid(internal_buffer).all(): internal_buffer= pygeos.make_valid(internal_buffer)
            external_buffer = pygeos.buffer(pygeos.get_exterior_ring(block), radius=buffer_radius)
            if not pygeos.is_valid(external_buffer).all(): external_buffer= pygeos.make_valid(external_buffer)

            complete_buffer = pygeos.get_parts(pygeos.union(internal_buffer, external_buffer))
            exterior_access = streets[pygeos.intersects(streets, external_buffer)]
            complete_buffer = complete_buffer[pygeos.intersects(complete_buffer, pygeos.multilinestrings(pygeos.get_parts(exterior_access)))]
            
            # Complete buffer is geometry covering streets connected to exterior
            complete_buffer = pygeos.union_all(complete_buffer) 
            if not pygeos.is_valid(complete_buffer).all(): complete_buffer = pygeos.make_valid(complete_buffer)

            off_network_geometry = pygeos.difference(street_multilinestrings, complete_buffer)
            off_network_length = sum(pygeos.length(off_network_geometry))
            on_network_geometry = pygeos.intersection(street_multilinestrings, complete_buffer)
            on_network_length = sum(pygeos.length(on_network_geometry))
            on_network_buffer = pygeos.buffer(on_network_geometry, radius = 1)
            if on_network_length > 0: 
                is_connected = True
        else: 
            on_network_length = 0
            off_network_length = 0
    else:
        nearest_external_street = np.nan
        on_network_length = 0
        off_network_length = 0

    # K-complexity calculation
    bldg_count = np.sum(pygeos.get_num_geometries(buildings)) 
    if bldg_count not in [1,0]:

        voronoi = pygeos.get_parts(pygeos.voronoi_polygons(geometry=buildings, extend_to=block))
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
                block_layers.append(str(np.sum(pygeos.get_num_geometries(block_parcels_outer))))
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
        data = gpd.GeoDataFrame.from_dict({'block_id': [block_id], 'on_network_street_length': float(on_network_length), 'off_network_street_length': float(off_network_length), 'nearest_external_street': float(nearest_external_street), 'building_count': int(bldg_count), 'building_layers': ','.join(block_layers), 'k_complexity': int(block_depth), 'geometry': gpd.GeoSeries(pygeos.to_shapely(block))}).set_crs(epsg=srid)
        data = data.to_crs(4326)
    else:
        data = pd.DataFrame.from_dict({'block_id': [block_id], 'on_network_street_length': float(on_network_length), 'off_network_street_length': float(off_network_length), 'nearest_external_street': float(nearest_external_street), 'building_count': int(bldg_count), 'building_layers': ','.join(block_layers), 'k_complexity': int(block_depth)})
    return data

def compute_layers(block_id: str, block: gpd.GeoSeries, buildings: gpd.GeoSeries, streets: gpd.GeoSeries, use_building_polygons: bool = False, buffer_radius: float = 100, srid: int = 3395, shrink_value: float = 0.4, segment_value: float = 0.5, threshold_value: float = 0.01) -> gpd.GeoDataFrame:
    """
    Computes the k-complexity value for each layer of buildings within an enclosed block.
    Differs from compute_k() because provides a detailed rendering of the internal layers of building access.
    Args:
        block_id: string, unique identification code
        block: GeoSeries of a single Polygon representing enclosed block
        buildings: GeoSeries of Points or Polygons representing buildings
        streets: GeoSeries of Linestrings or MultiLinestrings representing street networks
        use_building_polygons: bool (optional), defaults to False to use default triangular voronoi. Set to True to contour parcels around building polygons (slower performance).
        buffer_radius: float (optional), extends disconnected street_linestrings to street network that are within a set tolerance in meters
        srid: integer (optional), defaults to 3395, representing spatial reference identifier, for meter conversion
        For info about shrink_value, segment_value, threshold_value see: http://docs.momepy.org/en/stable/generated/momepy.Tessellation.html
    Returns:
        GeoDataFrame with building-parcel and street geometries. Geometry output projected in CRS 4326.
    """
    block_layers = [] 
    block_layers_geometry = []
    on_network_geometry = pygeos.Geometry("LINESTRING EMPTY")
    off_network_geometry = pygeos.Geometry("LINESTRING EMPTY")

    data = gpd.GeoDataFrame.from_dict({
        'block_id': pd.Series(dtype='str'), 
        'block_property':  pd.Series(dtype='str'), 
        'building_count': pd.Series(dtype='int'), 
        'k_complexity': pd.Series(dtype='int'), 
        'geometry': gpd.GeoSeries(dtype='geometry')
    }).set_crs(epsg=srid)

    # Check block geometry
    assert block.geom_type.unique().all() == 'Polygon', 'block geometry is not a Polygon.'
    block = pygeos.set_srid(pygeos.from_shapely(block.to_crs(epsg=srid)), srid)
    if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)
    if pygeos.get_num_geometries(block)[0] > 1:
        block = pygeos.get_parts(block)
        block = np.array([x for x in block if (pygeos.get_type_id(x) == 3 or pygeos.get_type_id(x) == 6)])
        block = np.array([pygeos.buffer(pygeos.multipolygons(pygeos.get_parts(block)),1)])
        if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)

    # Check buildings geometry
    buildings = pygeos.set_srid(pygeos.from_shapely(buildings.to_crs(epsg=srid)), srid)
    buildings = buildings[np.isin(pygeos.get_type_id(buildings), [0,3,4,6])]
    buildings = pygeos.get_parts(buildings)

    if np.sum(pygeos.get_num_geometries(buildings)) > 0:
        if use_building_polygons == False:
            if np.isin(pygeos.get_type_id(buildings), [3]).any(): 
                warnings.warn(f"Converting building Polygons to Points. Specify use_building_polygons = True to use Polygons.")
                buildings = buildings.centroid
            type_count = dict(collections.Counter(pygeos.get_type_id(buildings)))
            assert type_count.get(0) > 0, 'No Points present in buildings geopandas.GeoSeries.'
            buildings = buildings[np.isin(pygeos.get_type_id(buildings), [0])]
            buildings = pygeos.multipoints(buildings)
        else: 
            type_count = dict(collections.Counter(pygeos.get_type_id(buildings)))
            assert type_count.get(3) > 0, 'No Polygons present in buildings geopandas.GeoSeries.'
            buildings = buildings[np.isin(pygeos.get_type_id(buildings), [3])]
            buildings = pygeos.multipolygons(buildings)
        
    if not pygeos.is_valid(buildings).all(): buildings = pygeos.make_valid(buildings)

    bldg_count = np.sum(pygeos.get_num_geometries(buildings))
    is_connected = False

    # Check street networks
    if streets.geom_type.unique().any() in ['MultiLineString', 'LineString']:
        streets = pygeos.set_srid(pygeos.from_shapely(streets.to_crs(epsg=srid)), srid)
        if not pygeos.is_valid(streets).all(): streets = pygeos.make_valid(streets)
        streets = streets[np.isin(pygeos.get_type_id(streets), [1,2,5])]

        street_multilinestrings = pygeos.multilinestrings(pygeos.get_parts(streets))
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
        street_multilinestrings = pygeos.intersection(street_multilinestrings, block)
        if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)

        if not pygeos.is_empty(street_multilinestrings).all():
            
            street_multilinestrings = pygeos.line_merge(street_multilinestrings)
            if not pygeos.is_valid(street_multilinestrings).all(): street_multilinestrings = pygeos.make_valid(street_multilinestrings)
            internal_buffer = pygeos.buffer(street_multilinestrings, radius=buffer_radius/2)
            if not pygeos.is_valid(internal_buffer).all(): internal_buffer= pygeos.make_valid(internal_buffer)
            external_buffer = pygeos.buffer(pygeos.get_exterior_ring(block), radius=buffer_radius)
            if not pygeos.is_valid(external_buffer).all(): external_buffer= pygeos.make_valid(external_buffer)

            complete_buffer = pygeos.get_parts(pygeos.union(internal_buffer, external_buffer))
            exterior_access = streets[pygeos.intersects(streets, external_buffer)]
            complete_buffer = complete_buffer[pygeos.intersects(complete_buffer, pygeos.multilinestrings(pygeos.get_parts(exterior_access)))]
            
            # Complete buffer is geometry covering streets connected to exterior
            complete_buffer = pygeos.union_all(complete_buffer) 
            if not pygeos.is_valid(complete_buffer).all(): complete_buffer = pygeos.make_valid(complete_buffer)

            off_network_geometry = pygeos.difference(street_multilinestrings, complete_buffer)
            on_network_geometry = pygeos.intersection(street_multilinestrings, complete_buffer)
            on_network_buffer = pygeos.buffer(on_network_geometry, radius = 1)
            if sum(pygeos.length(on_network_geometry)) > 0: 
                is_connected = True
    
    if bldg_count not in [1,0]:

        if use_building_polygons == True:
            buildings_gdf = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(pygeos.to_shapely(buildings)) }).set_crs(epsg=srid)
            buildings_gdf = buildings_gdf.assign(index_id = buildings_gdf.index.astype(int))
            buildings_gdf['index_id'] = momepy.unique_id(buildings_gdf)
            block_gdf = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(pygeos.to_shapely(block)) }).set_crs(epsg=srid)
            buildings_gdf = buildings_gdf.assign(block_id = block_id.astype(str))
            enclosed_tess = momepy.Tessellation(gdf=buildings_gdf, unique_id='index_id', enclosures=block_gdf, enclosure_id='block_id', shrink = shrink_value, segment = segment_value, threshold = threshold_value).tessellation
            voronoi = pygeos.get_parts(pygeos.from_shapely(enclosed_tess['geometry']))
        else:
            voronoi = pygeos.get_parts(pygeos.voronoi_polygons(geometry=buildings, extend_to=block))

        if not pygeos.is_valid(block).all(): block = pygeos.make_valid(block)
        if not pygeos.is_valid(voronoi).all(): voronoi = pygeos.make_valid(voronoi)
        block_parcels = pygeos.intersection(block, voronoi)    
        if not pygeos.is_valid(block_parcels).all(): block_parcels = pygeos.make_valid(block_parcels)  
        bldg_count = np.sum(pygeos.get_num_geometries(block_parcels))

        if is_connected is True:
            block_intersect = pygeos.intersects(on_network_buffer, block_parcels)
            block_parcels_outer = block_parcels[block_intersect]
            block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
            #block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
            block_depth = len(block_layers)
            block_layers_geometry = pygeos.multipolygons(pygeos.get_parts(block_parcels_outer))
            building_count = np.sum(pygeos.get_num_geometries(block_parcels_outer))
            data_layer = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'block_property': 'building-parcels', 'building_count': building_count, 'k_complexity': block_depth, 'geometry':  gpd.GeoSeries(pygeos.to_shapely(block_layers_geometry))}).set_crs(epsg=srid)
            data = pd.concat([data, data_layer])

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
                block_depth = len(block_layers)
                block_layers_geometry = pygeos.multipolygons(pygeos.get_parts(block_parcels_outer))
                building_count = np.sum(pygeos.get_num_geometries(block_parcels_outer))
                data_layer = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'block_property': 'building-parcels', 'building_count': building_count, 'k_complexity': block_depth, 'geometry': gpd.GeoSeries(pygeos.to_shapely(block_layers_geometry)) }).set_crs(epsg=srid)
                data = pd.concat([data, data_layer])
                break
            block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
            block_depth = len(block_layers)
            block_layers_geometry = pygeos.multipolygons(pygeos.get_parts(block_parcels_outer))
            building_count = np.sum(pygeos.get_num_geometries(block_parcels_outer))
            data_layer = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'block_property': 'building-parcels', 'building_count': building_count, 'k_complexity': block_depth, 'geometry': gpd.GeoSeries(pygeos.to_shapely(block_layers_geometry)) }).set_crs(epsg=srid)
            data = pd.concat([data, data_layer])
    else:
        block_layers.append(bldg_count)
        block_depth = len(block_layers)
        block_layers_geometry = pygeos.multipolygons(pygeos.get_parts(block))
        data_layer = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'block_property': 'building-parcels', 'building_count': bldg_count, 'k_complexity': block_depth, 'geometry':  gpd.GeoSeries(pygeos.to_shapely(block_layers_geometry)) }).set_crs(epsg=srid)
        data = pd.concat([data, data_layer])

    if not pygeos.is_empty(on_network_geometry).all():
        lines = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'block_property': 'on-network-streets', 'building_count': np.nan, 'k_complexity': np.nan, 'geometry':  gpd.GeoSeries(pygeos.to_shapely(pygeos.multilinestrings(pygeos.get_parts(on_network_geometry)))) }).set_crs(epsg=srid)
        data = pd.concat([data, lines])
    if not pygeos.is_empty(off_network_geometry).all():
        lines = gpd.GeoDataFrame.from_dict({'block_id': block_id, 'block_property': 'off-network-streets', 'building_count': np.nan, 'k_complexity': np.nan, 'geometry':  gpd.GeoSeries(pygeos.to_shapely(pygeos.multilinestrings(pygeos.get_parts(off_network_geometry)))) }).set_crs(epsg=srid)
        data = pd.concat([data, lines])
    data = data.to_crs(4326)
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
        block_id_universe = country_blocks[['block_id']]
        del country_blocks

        # Read in building points
        country_buildings = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_points_{country_code}.parquet')
        block_id_universe_buildings = country_buildings[['block_id']]   
        block_id_universe_buildings = block_id_universe_buildings.drop_duplicates()
        block_id_universe_buildings = block_id_universe_buildings.assign(buildings_present=int(1))
        block_id_universe = pd.merge(left = block_id_universe, right = block_id_universe_buildings, how='left', on='block_id')
        block_id_universe['buildings_present'] = block_id_universe['buildings_present'].fillna(0)
        del block_id_universe_buildings
        
        # Remove completed blocks
        dask_folder_exists = os.path.isdir(Path(dask_dir) / f'{country_code}.parquet')
        if dask_folder_exists: 
            if len(os.listdir(Path(dask_dir) / f'{country_code}.parquet')) >= 1:            
                completed_blocks = dask.dataframe.read_parquet(path = Path(dask_dir) / f'{country_code}.parquet').compute()
                completed_block_list = list(completed_blocks['block_id'].unique())
                block_id_universe_completed = completed_blocks[['block_id']]
                block_id_universe_completed = block_id_universe_completed.drop_duplicates()
                block_id_universe_completed = block_id_universe_completed.assign(completed=int(1)) 
                block_id_universe = pd.merge(left = block_id_universe, right = block_id_universe_completed, how='left', on='block_id')
                block_id_universe['completed'] = block_id_universe['completed'].fillna(0)
                country_buildings = country_buildings[~country_buildings['block_id'].isin(completed_block_list)]
                logging.info(f"Completed blocks: {len(completed_block_list)} {round(len(completed_block_list)/block_id_universe.shape[0],2)*100} %")
                del completed_blocks, completed_block_list, block_id_universe_completed
            else:
                block_id_universe = block_id_universe.assign(completed=int(0)) 
        else:
            block_id_universe = block_id_universe.assign(completed=int(0)) 

        # Reconcile buildings and completed blocks with full universe
        logging.info(f"Blocks total: {block_id_universe.shape[0]}")
        logging.info(f"Blocks with buildings present: {block_id_universe.loc[(block_id_universe['buildings_present'] == 0)].shape[0]}")
        logging.info(f"Blocks without buildings present: {block_id_universe.loc[(block_id_universe['buildings_present'] == 1)].shape[0]}")
        logging.info(f"Blocks completed: {block_id_universe.loc[(block_id_universe['completed'] == 1)].shape[0]}")
        logging.info(f"Blocks not completed: {block_id_universe.loc[(block_id_universe['completed'] == 0)].shape[0]}")
        block_list = block_id_universe.loc[(block_id_universe['completed'] == 0) & (block_id_universe['buildings_present'] == 1)]['block_id'].unique()
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
        del country_buildings, cuts, buildings_count, block_id_universe, block_list
        gc.collect()
        
        # Read in streets
        streets = gpd.read_parquet(path = Path(streets_dir) / f'streets_{country_code}.parquet', memory_map = True).to_crs(4326)
        streets = streets[streets['geometry'].notnull()]

        # Loop over the chunks
        for i in chunk_dict.items():

            chunk_list = i[1]
            
            logging.info(f"Node status: {mem_profile_detail()}")
            logging.info(f"Partition {i[0]} of {len(chunk_dict)}")
            logging.info(f"Blocks in chunk: {len(chunk_list)}")

            # Read in block data within chunk
            try:
                blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True, filters = [('block_id', 'in', chunk_list)]).to_crs(4326)
                blocks = blocks[['block_id', 'block_geohash', 'gadm_code', 'country_code', 'geometry']]
            except Warning:
                logging.warning(f"Warning: No block data available for: {chunk_list}")

            logging.info(f"GADMs in chunk: {blocks['gadm_code'].unique()}")
            logging.info(f"Geohashes in chunk: {blocks['block_geohash'].str.slice(start=0, stop=4).unique()}")
            
            # Read in streets and limit to streets that intersect chunk
            try:
                blocks_buffered = blocks[['block_id','geometry']].copy()
                # Buffer set to 100 - is upper bound for nearest_external_street value in compute_k()
                blocks_buffered['geometry'] = blocks_buffered['geometry'].to_crs(3395).buffer(100).to_crs(4326)
                street_network = gpd.overlay(df1 = streets, df2 = blocks_buffered, how='intersection', keep_geom_type=True, make_valid=True)
                street_network = street_network[['block_id', 'geometry']]
                del blocks_buffered
            except Warning:
                logging.warning(f"Warning: No street data available for: {chunk_list}")

            # Read in buildings within chunk
            try:
                buildings = gpd.read_parquet(path = Path(buildings_dir) / f'buildings_points_{country_code}.parquet', memory_map = True, filters = [('block_id', 'in', chunk_list)]).to_crs(4326)
            except Warning:
                logging.warning(f"Warning: No building data available for: {chunk_list}")
            logging.info(f"Buildings in chunk: {buildings.shape[0]}")

            # Initialize dataframe  
            k_output = pd.DataFrame({
                'block_id': pd.Series(dtype='str'), 
                'on_network_street_length': pd.Series(dtype='float'),  
                'off_network_street_length': pd.Series(dtype='float'), 
                'nearest_external_street': pd.Series(dtype='float'), 
                'building_count': pd.Series(dtype='int'), 
                'building_layers': pd.Series(dtype='object'), 
                'k_complexity': pd.Series(dtype='int')
            })   

            # Generate inputs to prevent copying data across workers
            inputs = inputs_generator(block_id_col = 'block_id', blocks = blocks, buildings = buildings, streets = street_network)

            # Parallelize computation
            with multiprocessing.Pool(processes = number_of_cores, maxtasksperchild = 100) as pool:
                results = pool.starmap(func = compute_k, iterable = inputs, chunksize = 10)
            for i in results: k_output = pd.concat([k_output, i])

            # Incremental file build
            k_output = dask.dataframe.from_pandas(data = k_output, npartitions = 1) 
            dask.dataframe.to_parquet(df = k_output, path = Path(dask_dir) / f'{country_code}.parquet', engine='pyarrow', compression='snappy', append=True, ignore_divisions=True)
            del blocks, street_network, buildings, k_output, results
            #gc.collect()

        # Combine partitioned dataframes into single parquet
        k_bulk = dask.dataframe.read_parquet(path = Path(dask_dir) / f'{country_code}.parquet').compute()
        k_bulk.to_parquet(path = Path(complexity_dir) / f'complexity_{country_code}.parquet')

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
