import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
from typing import List, Union
gpd.options.use_pygeos = True

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
import argparse
np.set_printoptions(suppress=True)
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

import pyarrow
import dask 
import dask_geopandas


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

def mem_profile() -> str: 
    """
    Return memory usage, str
    """
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use

country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}

def code_to_geofabrik(country_code: str, country_dict: dict = country_dict) -> str:
    """
    Map country_codes to Geofabrik country name
    Args:
        country_code: str, country code in ISO 3166-1 alpha-3 format
        country_dict: dictionary, maps country codes to Geofabrik file names
    Returns:
        Name of Geofabrik country for a corresponding country_code, str
    """
    geofabrik_name = [value[0] for (key, value) in country_dict.items() if key == country_code][0]
    return geofabrik_name

def geofabrik_to_code(geofabrik_name: str, country_dict: dict = country_dict) -> str:
    """
    Map Geofabrik file names to country_codes in ISO 3166-1 alpha-3 format
    Args:
        geofabrik_name: str, Geofabrik country name
        country_dict: dictionary, maps country codes to Geofabrik file names
    Returns:
        File name of country_code for a corresponding Geofabrik country name, str
        Returns more than one country code for the following countries {'israel-and-palestine' : ['PSE', 'ISR'], 'malaysia-singapore-brunei' : ['BRN', 'MYS', 'SGP'], 'morocco' : ['MAR', 'ESH'], 'senegal-and-gambia' : ['GMB', 'SEN']}
    """
    country_code = [key for (key, value) in country_dict.items() if value[0] == geofabrik_name]
    return country_code

def read_osm(country_code: str, directory_path: Union[str, Path]) -> gpd.GeoDataFrame:
    """
    Read in OSM data ingested from Geofabrik using the download_job.sh script and convert to pipeline-ready format
    Args:
        country_code: str, country codes in ISO 3166-1 alpha-3 format
        directory_path: str or Path, parent directory holding OSM linestring data from Geofabrik
    Returns:
        GeoDataFrame containing OSM linestrings for a particular country
    """
    geofabrik_name = code_to_geofabrik(country_code = country_code)
    data = gpd.read_file(Path(directory_path) / f'{geofabrik_name}-latest-linestring.geojson')
    data = data.rename(columns={'@type':'type','@id':'id','@version':'version','@changeset':'changeset','@timestamp':'timestamp'})
    return data

def main(log_file: Path, country_chunk: list, osm_dir: Path, gadm_dir: Path, output_dir: Path):

    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f"Countries to process: {country_chunk}")
    
    # Make directories
    block_dir =  str(output_dir) + '/blocks'
    Path(block_dir).mkdir(parents=True, exist_ok=True)
    street_dir =  str(output_dir) + '/streets'
    Path(street_dir).mkdir(parents=True, exist_ok=True)
    logging.info(f"block_dir: {block_dir}")
    logging.info(f"street_dir: {street_dir}")

    # Check if country is completed
    output_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(block_dir))))))
    output_country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in output_file_list] 
    logging.info(f"Finished countries: {output_country_list}")
    
    # Subset to country list (remove invalid country codes and remove completed countries)
    if country_chunk: 
        country_list = list(country_dict.keys())
        country_list = [x for x in country_chunk if x in set(country_list)]
        country_list = [x for x in country_chunk if x not in set(output_country_list)]
        if not country_list: 
            raise ValueError('Empty country list')

    logging.info(f"Remaining countries: {country_list}")

    # Confirm remaining countries have data in input files
    osm_list_of_lists = list(filter(re.compile("-latest-linestring.geojson").search, sorted(list(os.listdir(Path(osm_dir))))))
    osm_list_of_lists = list(map(geofabrik_to_code, [re.sub('-latest-linestring.geojson', '', i) for i in osm_list_of_lists]))
    osm_inputs_list = [i for flatlist in osm_list_of_lists for i in flatlist]
    gadm_inputs_list = list(filter(re.compile("gadm_").match, sorted(list(os.listdir(Path(gadm_dir))))))
    gadm_inputs_list = [(re.sub('gadm_', '', re.sub('.geojson', '', i))) for i in gadm_inputs_list] 

    in_chunk_not_in_osm_inputs = [x for x in country_list if x not in set(osm_inputs_list)]
    if len(in_chunk_not_in_osm_inputs) > 0:
        raise ValueError(f'OSM input data does not exist for {in_chunk_not_in_osm_inputs} in country_chunk arg.')
    in_chunk_not_in_gadm_inputs = [x for x in country_list if x not in set(gadm_inputs_list)]
    if len(in_chunk_not_in_gadm_inputs) > 0:
        raise ValueError(f'GADM input data does not exist for {in_chunk_not_in_gadm_inputs} in country_chunk arg.')

    logging.info(f"Generate block geometries")

    # Process country list
    for country_code in country_list: 
        logging.info(f"Processing: {country_code}")
        
        # Read OSM files
        t0 = time.time()
        osm_gpd = read_osm(country_code = country_code, directory_path = osm_dir) 
        osm_gpd = osm_gpd.explode(ignore_index = True)
        osm_pygeos = prepare.from_shapely_srid(geometry = osm_gpd, srid = 4326) 
        t1 = time.time()
        logging.info(f"Read OSM: osm_gpd: {osm_gpd.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
    
        # Read GADM files
        t0 = time.time()
        gadm_gpd = gpd.read_file(Path(gadm_dir) / f'gadm_{country_code}.geojson')
        gadm_col = max(list(filter(re.compile("GID_*").match, list(gadm_gpd.columns))))
        t1 = time.time()
        logging.info(f"Read GADM: gadm_gpd: {gadm_gpd.shape}, gadm_col: {gadm_col}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        
        # Trim coastline if applicable
        if osm_gpd[osm_gpd['natural'].isin(['coastline'])].shape[0] > 0:
            logging.info(f"Trim coastline")
            t0 = time.time()
            gadm_gpd_trim = prepare.trim_coastline(gadm_data = gadm_gpd, osm_data = osm_gpd)
            trimmed_area_percent = (sum(gadm_gpd.to_crs(3395).area/10**6)-sum(gadm_gpd_trim.to_crs(3395).area/10**6))/sum(gadm_gpd_trim.to_crs(3395).area/10**6)
            logging.info(f"Removed area with trim_coastline(): {round(trimmed_area_percent,3)} %")
            gadm_gpd_trim = gpd.clip(gdf = gadm_gpd['geometry'], mask = gadm_gpd_trim.unary_union, keep_geom_type=True)
            gadm_gpd_trim = gpd.GeoDataFrame(geometry=gpd.GeoSeries(gadm_gpd_trim)).reset_index(drop=True)
            gadm_gpd_trim = gpd.overlay(df1 = gadm_gpd_trim, df2 = gadm_gpd, keep_geom_type=True)
            gadm_dropped = [x for x in gadm_gpd[gadm_col].unique() if x not in set(gadm_gpd_trim[gadm_col].unique())]
            gadm_gpd_trim = gadm_gpd_trim.dissolve(by=[gadm_col], as_index = False)
            t1 = time.time()
            logging.info(f"Dropped GADMs from trim_coastline(): {gadm_dropped}")
            logging.info(f"Trimmed gadm_gpd: {gadm_gpd.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
            gadm_gpd_union = gadm_gpd_trim.unary_union
            gadm_gpd_explode = gadm_gpd.explode(ignore_index = True)
            gadm_gpd_residual = gadm_gpd_explode.disjoint(gadm_gpd_union, align=True)
            gadm_gpd = pd.concat([gadm_gpd_trim, gadm_gpd_explode[gadm_gpd_residual]], ignore_index=True)
            gadm_gpd = gadm_gpd.dissolve(by=[gadm_col], as_index = False)
            del gadm_gpd_trim, gadm_gpd_union, gadm_gpd_explode, gadm_gpd_residual
            logging.info(f"Trimmed with residuals gadm_gpd: {gadm_gpd.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
            
        # Write street geometries
        t0 = time.time()
        osm_streets = osm_gpd[osm_gpd['highway'].notnull()]
        osm_streets.to_parquet(Path(street_dir) / f'streets_{country_code}.parquet', compression='snappy')
        t1 = time.time()
        logging.info(f"Writing osm_streets ({osm_streets.shape}) from osm_gpd ({osm_gpd.shape}) ./streets_{country_code}.parquet, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        del osm_gpd
        
        gadm_list = list(gadm_gpd[gadm_col].unique())
        
        # Initialize
        block_bulk = gpd.GeoDataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 

        # Build blocks for each GADM
        logging.info(f"Build blocks")

        t0 = time.time()
        output_map = list(map(lambda x: prepare.build_blocks(gadm_data = gadm_gpd, osm_data = osm_pygeos, gadm_column = gadm_col, gadm_code = x), gadm_list))
        for i,j in enumerate(output_map): block_bulk = pd.concat([block_bulk, output_map[i]], ignore_index=True)
        #check_area = np.sum(gadm_blocks['geometry'].to_crs(3395).area)/np.sum(gadm_gpd[gadm_gpd[gadm_col] == i]['geometry'].to_crs(3395).area)
        #if (check_area/1) < .99: logging.info(f"Area proportion: {round(check_area,4)}")
        t1 = time.time()
        logging.info(f"build_blocks() function: block_bulk: {block_bulk.shape}, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        
        # Write block geometries
        t0 = time.time()
        block_bulk.to_parquet(Path(block_dir) / f'blocks_{country_code}.parquet', compression='snappy')
        t1 = time.time()
        logging.info(f"Writing ./blocks_{country_code}.parquet, {mem_profile()}, {str(round(t1-t0,3))} seconds")
        logging.info(f"Finished {country_code}.")
        
    logging.info(f"Finished.")

def setup(args=None):    
    parser = argparse.ArgumentParser(description='Build blocks.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format")
    parser.add_argument('--osm_dir', required=True, type=Path, dest="osm_dir", help="OSM directory")
    parser.add_argument('--gadm_dir', required=True, type=Path, dest="gadm_dir", help="GADM directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Output directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))


