
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
import re 
import os 
from typing import List, Union
from pathlib import Path

gpd.options.use_pygeos = True
pd.options.mode.chained_assignment = None 

# print to logs as well as terminal nbETV0m169rs

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

def trim_coastline(gadm_data: gpd.GeoDataFrame, osm_data: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Trim the geometries along coastline which can have 
    pixelated boundaries that does not match contours of 
    natural features. 
    Args:
        gadm_data: GeoDataFrame containing GADM delineations, requires CRS WGS 84 EPSG 4326
        osm_data: GeoDataFrame containing linestrings coastline features, requires CRS WGS 84 EPSG 4326
    Returns:
        GeoDataFrame with polygon(s) of the input geometry with trimmed coastline.
        Geometry projected in CRS WGS 84 EPSG 4326. 
    """
    assert gadm_data.crs == 'epsg:4326', "gadm_data is not epsg:4326."
    assert osm_data.crs == 'epsg:4326', "osm_data is not epsg:4326."

    bbox = pygeos.envelope(pygeos.union_all(pygeos.from_shapely(gadm_data['geometry'])))
    bbox_center = pygeos.get_coordinates(pygeos.centroid(bbox)).tolist()
    bbox_enlarge = pygeos.apply(bbox, lambda x: ((x - bbox_center)*1.1 + bbox_center) )

    land = pygeos.from_shapely(gadm_data['geometry'])
    land_center = pygeos.buffer(pygeos.centroid(land), radius=.00001, quadsegs=50) 
    coast = pygeos.union_all(pygeos.from_shapely(osm_data[osm_data['natural'].isin(['coastline','water'])]['geometry']))
    coast_buffer = pygeos.buffer(coast, radius=.00001, quadsegs=50) 
    bbox_poly = pygeos.difference(bbox_enlarge,coast_buffer)
    bbox_poly_parts = pygeos.get_parts(bbox_poly)

    positive_space = bbox_poly_parts[pygeos.intersects(pygeos.multipolygons(land_center), bbox_poly_parts)]
    land_intersect = pygeos.get_parts(pygeos.intersection(pygeos.union_all(land), positive_space))
    land_intersect = land_intersect[pygeos.intersects(pygeos.multipolygons(land_center), land_intersect)]
     
    data = gpd.GeoDataFrame({"geometry": gpd.GeoSeries(pygeos.to_shapely(land_intersect))}).set_crs(epsg=4326) 
    return data

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
        Geometry projected in CRS WGS 84 EPSG 4326.
    """
    assert gadm_data.crs == 'epsg:4326', "gadm_data is not epsg:4326."
    if type(osm_data) == gpd.geodataframe.GeoDataFrame: 
        osm_data = from_shapely_srid(geometry = osm_data, srid = 4326)
    assert len(np.unique(pygeos.get_srid(osm_data))) == 1 and np.unique(pygeos.get_srid(osm_data))[0] == 4326, "osm_data is not epsg:4326."

    gadm_data = pygeos.from_shapely(gadm_data[gadm_data[gadm_column] == gadm_code]['geometry'])
    osm_data = pygeos.intersection(pygeos.multilinestrings(osm_data),gadm_data)
    gadm_outline = pygeos.get_exterior_ring(pygeos.get_parts(gadm_data)[compute_area(pygeos.get_parts(gadm_data)) == max(compute_area(pygeos.get_parts(gadm_data)))])
    gadm_lines = pygeos.line_merge(pygeos.multilinestrings(pygeos.get_exterior_ring(pygeos.get_parts(gadm_data))))
    gadm_lines = pygeos.union(osm_data[0],gadm_lines)
    gadm_blocks = pygeos.union(osm_data[0],gadm_outline)
    gadm_blocks = pygeos.polygonize_full([gadm_blocks])
    gadm_blocks = gpd.GeoDataFrame.from_dict({"country_code": gadm_code[0:3],
                                              "gadm_code": gadm_code,
                                              'geometry': pygeos.to_shapely(pygeos.get_parts(pygeos.normalize(gadm_blocks[0])))}, crs=4326) 
    gadm_blocks = gadm_blocks.reset_index(drop=True)
    gadm_blocks = gadm_blocks.assign(block_id = [gadm_code + '_' + str(x) for x in list(gadm_blocks.index)])
    gadm_blocks = gadm_blocks[['block_id','gadm_code','country_code','geometry']].to_crs(epsg=4326)
    return gadm_blocks

def index_buildings(gadm_block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Joins GeoDataFrame of enclosed block geometries to a 
    GeoDataFrame containing building footprints, for footprints that 
    overlap multiple with block geometries this function allocates 
    the footprint to the block containing its centroid
    Args:
        gadm_block_data: GeoDataFrame output returned from build_blocks() function, requires CRS WGS 84 EPSG 4326
        bldg_data: GeoDataFrame containing building geometries, requires CRS WGS 84 EPSG 4326
    Returns:
        GeoDataFrame with building geometries mapped to 'block_id','gadm_code','country_code'.
        Geometry projected in CRS WGS 84 EPSG 4326. 
    """
    assert gadm_block_data.crs == 'epsg:4326', "gadm_block_data is not epsg:4326."
    assert bldg_data.crs == 'epsg:4326', "bldg_data is not epsg:4326."

    gadm_block_data = gadm_block_data.to_crs(epsg=3395)
    bldg_data = bldg_data.to_crs(epsg=3395)
    bldg_points = bldg_data.copy()
    bldg_points['geometry'] = bldg_points.centroid
    bldg_index = bldg_points.sindex

    index_bulk = bldg_index.query_bulk(gadm_block_data['geometry'], predicate="intersects")
    blocks_buildings_map = pd.DataFrame({'index_blocks': index_bulk[0], 'index_buildings': index_bulk[1]})
    blocks_buildings_map = blocks_buildings_map.merge(bldg_data[['geometry']], how = 'left', left_on='index_buildings', right_index=True)
    blocks_buildings_map = blocks_buildings_map.merge(gadm_block_data[['block_id','gadm_code','country_code']], left_on='index_blocks', right_index=True)

    data = gpd.GeoDataFrame(blocks_buildings_map[['block_id','gadm_code','country_code','geometry']]).set_crs(epsg=3395)
    data = data.to_crs(epsg=4326)
    return data

def compute_k(block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame, block_col: str, block_id: str, street_linestrings: Union[pygeos.Geometry, gpd.GeoDataFrame] = None, buffer_streets: bool = False) -> gpd.GeoDataFrame:
    """
    Computes the k-complexity value for each block as well as 
    block level metadata taking a GeoDataFrame of enclosed blocks, 
    the indexed building GeoDataFrame, and street linestrings in 
    the form of a PyGEOS Geometry (or GeoDataFrame with slower performance).
    Args:
        block_data: GeoDataFrame output returned from build_blocks() function, requires CRS WGS 84 EPSG 4326 (or EPSG 3395 for faster performance)
        bldg_data: GeoDataFrame output returned from index_buildings() function, requires CRS WGS 84 EPSG 4326 (or EPSG 3395 for faster performance)
        block_col: string, column name that contains the block_id codes (present in block_data and bldg_data)
        block_id: string, unique identification code in block_col
        street_linestrings: (optional) GeoDataFrame or PyGEOS Geometry array, linestrings representing street networks, requires CRS WGS 84 EPSG 4326 or 3395. For faster performance use EPSG 3395 PyGEOS Geometry array.
        buffer_streets: bool, default False, if True removes all disconnected streets in a buffer radius of 30
    Returns:
        GeoDataFrame with block geometries in row and the columns: 'block_id','block_area','building_area','building_count','building_layers','k_complexity'
        Geometry projected in CRS WGS 84 EPSG 4326.
    """
    assert block_data.crs in ['epsg:3395','epsg:4326'], "block_data is not epsg:4326 or epsg:3395."
    assert bldg_data.crs in ['epsg:3395','epsg:4326'], "bldg_data is not epsg:4326 or epsg:3395."
    
    if block_data.crs == 'epsg:4326': block_data = block_data.to_crs(epsg=3395)
    if bldg_data.crs == 'epsg:4326': bldg_data = bldg_data.to_crs(epsg=3395)

    block_layers = [] 
    bldg_data = bldg_data[bldg_data[block_col] == block_id]
    block_data = block_data[block_data[block_col] == block_id]

    country_code = block_data['country_code'].unique()[0] 
    gadm_code = block_data['gadm_code'].unique()[0] 
    
    block_area = round(pygeos.area(pygeos.multipolygons(pygeos.from_shapely(block_data['geometry']))),2)
    bldg_array = pygeos.from_shapely(bldg_data['geometry'])
    bldg_area_list = pygeos.area(bldg_array)
    bldg_area = round(sum(bldg_area_list),2)
    bldg_count = np.sum(pygeos.get_num_geometries(bldg_array))
        
    if bldg_count not in [1,0]:
        bldg_data["geometry"] = bldg_data.centroid
        building_points = pygeos.multipoints(pygeos.from_shapely(bldg_data["geometry"]))
        points = pygeos.get_parts(building_points)
        block = pygeos.from_shapely(block_data['geometry'])
        voronoi = pygeos.get_parts(pygeos.voronoi_polygons(geometry=building_points, extend_to=block))
        if not pygeos.is_valid(block): block = pygeos.make_valid(block)
        block_parcels = pygeos.intersection(block, voronoi)    
        bldg_count = np.sum(pygeos.get_num_geometries(block_parcels))

        if street_linestrings is not None:
            if type(street_linestrings) == gpd.geodataframe.GeoDataFrame:    
                street_linestrings = from_shapely_srid(geometry = street_linestrings, srid = 3395)
            assert type(street_linestrings) == np.ndarray and all(pygeos.is_geometry(street_linestrings)), 'street_linestrings is an invalid geometry type.'
            if np.unique(pygeos.get_srid(street_linestrings))[0] == 4326:
                street_linestrings = transform_crs(geometry = street_linestrings, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395")
            assert len(np.unique(pygeos.get_srid(street_linestrings))) == 1 and np.unique(pygeos.get_srid(street_linestrings))[0] == 3395, 'street_linestrings is not epsg:4326 or epsg:3395.' 
            street_linestrings = pygeos.intersection(pygeos.multilinestrings(street_linestrings),block)
            if buffer_streets == True:
                connected_streets = pygeos.get_parts(pygeos.buffer(street_linestrings, radius=30, quadsegs=5)) 
                connected_streets = connected_streets[pygeos.length(connected_streets).argsort()[::-1] == 0]
                exterior_boundary = pygeos.buffer(pygeos.get_exterior_ring(block), radius=60, quadsegs=5)
                complete_buffer = pygeos.union(connected_streets, exterior_boundary)
                street_length = pygeos.length(street_linestrings)
                street_linestrings = pygeos.intersection(street_linestrings,complete_buffer)
                street_length_trim = pygeos.length(street_linestrings)
                street_length_delta = round(street_length[0] - street_length_trim[0],3)
                if street_length_delta > 0:
                    print(f'Buffering: {block_id}, {street_length_delta} meters, {round((street_length_delta/street_length[0])*100,5)} %')    
            if pygeos.is_empty(street_linestrings[0]) == False:
                block_intersect = pygeos.intersects(street_linestrings,block_parcels)
                block_parcels_outer = block_parcels[block_intersect]
                block_layers.append( str(np.sum(pygeos.get_num_geometries(block_parcels_outer))))
                block_depth = len(block_layers)
                block_parcels = block_parcels[~block_intersect]
                if len(block_parcels_outer) == 0: 
                    block_parcels_outer = block_parcels.copy()
            else: block_parcels_outer = block_parcels.copy()

        parcel_residual = np.sum(pygeos.get_num_geometries(block_parcels))

        while parcel_residual > 0:

            block_reduced = pygeos.coverage_union_all(block_parcels_outer)
            block_parcels_inner = block_parcels[~pygeos.touches(block_parcels,block_reduced)]      
            block_parcels_outer = block_parcels[pygeos.touches(block_parcels,block_reduced)]  
            if np.sum(pygeos.get_num_geometries(block_parcels_outer)) == 0 and np.sum(pygeos.get_num_geometries(block_parcels_inner)) > 0:
                block_reduced = pygeos.coverage_union_all(block_parcels_inner)
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
        block_layers.append( str(bldg_count))
        block_depth = len(block_layers)
        
    data = gpd.GeoDataFrame.from_dict({'block_id': [block_id], 
                                       'gadm_code': gadm_code,
                                       'country_code': country_code,
                                       'block_area': float(block_area),
                                       'building_area': float(bldg_area),
                                       'building_count': int(bldg_count),
                                       'building_layers': ','.join(block_layers), 
                                       'k_complexity': int(block_depth), 
                                       'geometry': block_data['geometry']}).set_crs(epsg=3395)
    data = data.to_crs(4326)
    return data

def compute_layers(block_data: gpd.GeoDataFrame, bldg_data: gpd.GeoDataFrame, block_col: str, block_id: str, street_linestrings: Union[pygeos.Geometry, gpd.GeoDataFrame] = None, buffer_streets: bool = False) -> gpd.GeoDataFrame:
    """
    Computes the k-complexity value for each layer of buildings 
    within a block taking a GeoDataFrame of enclosed blocks, 
    the indexed building GeoDataFrame, and street linestrings in 
    the form of a PyGEOS Geometry (or GeoDataFrame with slower performance). 
    Differs from compute_k() because provides a detailed rendering 
    of the internal layers of building access.
    Args:
        block_data: GeoDataFrame, output returned from build_blocks() function, requires CRS WGS 84 EPSG 4326 or 3395 (for faster performance)
        bldg_data: GeoDataFrame, output returned from index_buildings() function, requires CRS WGS 84 EPSG 4326 or 3395 (for faster performance)
        block_col: string, column that contains the block_id codes (present in block_data and bldg_data)
        block_id: string, unique identification code of block_id string in block_col
        street_linestrings: GeoDataFrame or PyGEOS Geometry array, linestrings representing street accesses (optional), requires CRS WGS 84 4326 or 3395. For faster performance use EPSG 3395 PyGEOS Geometry array.
        buffer_streets: bool, default False, if True removes all disconnected streets in a buffer radius of 30
    Returns:
        GeoDataFrame with block-layer geometries and the columns: 'block_id','gadm_code','country_code','building_count','k_complexity'
        Geometry projected in CRS WGS 84 EPSG 4326.
    """
    assert block_data.crs in ['epsg:3395','epsg:4326'], "block_data is not epsg:4326 or epsg:3395."
    assert bldg_data.crs in ['epsg:3395','epsg:4326'], "bldg_data  is not epsg:4326 or epsg:3395."
        
    if block_data.crs == 'epsg:4326': block_data = block_data.to_crs(epsg=3395)
    if bldg_data.crs == 'epsg:4326': bldg_data = bldg_data.to_crs(epsg=3395)
    
    block_layers = [] 
    block_layers_geometry = []
    k_complexity = []
    bldg_data = bldg_data[bldg_data[block_col] == block_id]
    block_data = block_data[block_data[block_col] == block_id]
    country_code = block_data['country_code'].unique()[0] 
    gadm_code = block_data['gadm_code'].unique()[0] 
    bldg_array = pygeos.from_shapely(bldg_data['geometry'])
    bldg_count = np.sum(pygeos.get_num_geometries(bldg_array))
    
    if bldg_count not in [1,0]:
        bldg_data["geometry"] = bldg_data.centroid
        building_points = pygeos.multipoints(pygeos.from_shapely(bldg_data["geometry"]))
        points = pygeos.get_parts(building_points)
        block = pygeos.from_shapely(block_data['geometry'])
        voronoi = pygeos.get_parts(pygeos.voronoi_polygons(geometry=building_points, extend_to=block))
        if not pygeos.is_valid(block): block = pygeos.make_valid(block)
        block_parcels = pygeos.intersection(block, voronoi)    
        bldg_count = np.sum(pygeos.get_num_geometries(block_parcels))

        if street_linestrings is not None:
            if type(street_linestrings) == gpd.geodataframe.GeoDataFrame:    
                street_linestrings = from_shapely_srid(geometry = street_linestrings, srid = 3395)
            assert type(street_linestrings) == np.ndarray and all(pygeos.is_geometry(street_linestrings)), 'street_linestrings is an invalid geometry type.'
            if np.unique(pygeos.get_srid(street_linestrings))[0] == 4326:
                street_linestrings = transform_crs(geometry = street_linestrings, epsg_from = "EPSG:4326", epsg_to = "EPSG:3395")
            assert len(np.unique(pygeos.get_srid(street_linestrings))) == 1 and np.unique(pygeos.get_srid(street_linestrings))[0] == 3395, 'street_linestrings is not epsg:4326 or epsg:3395.' 
            street_linestrings = pygeos.intersection(pygeos.multilinestrings(street_linestrings),block)
            if buffer_streets == True:
                connected_streets = pygeos.get_parts(pygeos.buffer(street_linestrings, radius=30, quadsegs=5)) 
                connected_streets = connected_streets[pygeos.length(connected_streets).argsort()[::-1] == 0]
                exterior_boundary = pygeos.buffer(pygeos.get_exterior_ring(block), radius=60, quadsegs=5)
                complete_buffer = pygeos.union(connected_streets, exterior_boundary)
                street_length = pygeos.length(street_linestrings)
                street_linestrings = pygeos.intersection(street_linestrings,complete_buffer)
                street_length_trim = pygeos.length(street_linestrings)
                street_length_delta = round(street_length[0] - street_length_trim[0],3)
                if street_length_delta > 0:
                    print(f'Buffering: {block_id}, {street_length_delta} meters, {round((street_length_delta/street_length[0])*100,5)} %')    
            if pygeos.is_empty(street_linestrings[0]) == False:
                block_intersect = pygeos.intersects(street_linestrings,block_parcels)
                block_parcels_outer = block_parcels[block_intersect]
                block_layers.append(np.sum(pygeos.get_num_geometries(block_parcels_outer)))
                block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
                block_depth = len(block_layers)
                k_complexity.append(block_depth)
                block_parcels = block_parcels[~block_intersect]
                if len(block_parcels_outer) == 0: 
                    block_parcels_outer = block_parcels.copy()
            else: block_parcels_outer = block_parcels.copy()
 
        parcel_residual = np.sum(pygeos.get_num_geometries(block_parcels))
        
        while parcel_residual > 0:
            block_reduced = pygeos.coverage_union_all(block_parcels_outer)
            block_parcels_inner = block_parcels[~pygeos.touches(block_parcels,block_reduced)]      
            block_parcels_outer = block_parcels[pygeos.touches(block_parcels,block_reduced)] 
            if np.sum(pygeos.get_num_geometries(block_parcels_outer)) == 0 and np.sum(pygeos.get_num_geometries(block_parcels_inner)) > 0:
                block_reduced = pygeos.coverage_union_all(block_parcels_inner)
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
            block_depth = len(block_layers) 
            k_complexity.append(block_depth)
        else:
            block_depth = len(block_layers)
            k_complexity.append(block_depth)
    else:
        block_layers.append(bldg_count)
        block_layers_geometry.append(pygeos.geometrycollections(block_parcels_outer))
        block_depth = len(block_layers)
        k_complexity.append(block_depth)

    data = gpd.GeoDataFrame.from_dict({'block_id': block_id, 
                                       'gadm_code': gadm_code,
                                       'country_code': country_code,
                                       'building_count': block_layers,
                                       'k_complexity': k_complexity,
                                       'geometry': pygeos.to_shapely(block_layers_geometry)}).set_crs(epsg=3395)
    data = data.to_crs(4326)
    return data



