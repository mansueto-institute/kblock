import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
from typing import List, Union
gpd.options.use_pygeos = True

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
    coast = pygeos.line_merge(pygeos.union_all(pygeos.from_shapely(osm_data[osm_data['natural'].isin(['coastline'])]['geometry'])))
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
    gadm_lines = pygeos.line_merge(pygeos.multilinestrings(pygeos.get_exterior_ring(pygeos.get_parts(gadm_data))))
    gadm_blocks = pygeos.polygonize_full([pygeos.union(osm_data[0],gadm_lines)])[0]
    gadm_blocks = gpd.GeoDataFrame.from_dict({"country_code": gadm_code[0:3],"gadm_code": gadm_code,'geometry': pygeos.to_shapely(pygeos.get_parts(pygeos.normalize(pygeos.get_parts(gadm_blocks))))}).set_crs(4326)  
    gadm_blocks = gadm_blocks.reset_index(drop=True)
    gadm_blocks = gadm_blocks.assign(block_id = [gadm_code + '_' + str(x) for x in list(gadm_blocks.index)])
    gadm_blocks = gadm_blocks[['block_id','gadm_code','country_code','geometry']].to_crs(epsg=4326)
    return gadm_blocks


def index_streets(osm_data: Union[pygeos.Geometry, gpd.GeoDataFrame], gadm_data: gpd.GeoDataFrame, gadm_column: str) -> gpd.GeoDataFrame:
    """
    Map the OSM linestrings to GADM boundaries. Includes duplicates 
    for linestrings that cross multiple GADM boundaries.
    Args:
        gadm_data: GeoDataFrame, containing GADM delineations (downloaded from https://gadm.org/data.html), requires CRS WGS 84 EPSG 4326
        osm_data: GeoDataFrame (or PyGEOS Geometry array), of OSM linestrings requires CRS WGS 84 EPSG 4326
        gadm_column: string, name of column containing gadm_code
    Returns:
        GeoDataFrame with linestrings mapped to 'gadm_code'
        Geometry projected in CRS WGS 84 EPSG 4326.
    """
    assert gadm_data.crs == 'epsg:4326', "gadm_data is not epsg:4326."
    assert osm_data.crs == 'epsg:4326', "osm_data is not epsg:4326."

    osm_streets = osm_data[osm_data['highway'].notnull()]
    osm_streets_index = osm_streets.sindex
    
    index_1 = osm_streets_index.query_bulk(gadm_data['geometry'], predicate="contains")
    index_2 = osm_streets_index.query_bulk(gadm_data['geometry'], predicate="intersects")
    streets_map1 = pd.DataFrame({'index_gadm': index_1[0], 'index_streets': index_1[1]})
    streets_map2 = pd.DataFrame({'index_gadm': index_2[0], 'index_streets': index_2[1]})
    
    streets_map = pd.concat([streets_map1, streets_map2], ignore_index=True)
    streets_map = streets_map.drop_duplicates()
    streets_map = streets_map.merge(osm_streets[['id','geometry']], how = 'left', left_on='index_streets', right_index=True)
    streets_map = streets_map.merge(gadm_data[[gadm_column]], how = 'left', left_on='index_gadm', right_index=True)
    
    osm_mapped = gpd.GeoDataFrame(streets_map[['id',gadm_column,'geometry']]).set_crs(epsg=4326)
    osm_mapped = osm_mapped.rename(columns={gadm_column:"gadm_code"})
    return osm_mapped
