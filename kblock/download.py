import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
gpd.options.use_pygeos = True

import warnings
import os
import sys
import re
import io
from pathlib import Path
from typing import List, Union
from urlpath import URL
import requests
import tempfile
import zipfile

import pyrosm
from pyrosm.data import sources
from pyrosm import OSM, get_data
from pyrosm.config import Conf

def get_osm_lines(country_code: str, download_dir: Union[str, Path]) -> gpd.GeoDataFrame:
    """
    Download *.osm.pbf file from Geofabrik and parse into linestring geometries
    Args:
        country_code: three-letter country code following ISO 3166-1 alpha-3 format
    Returns:
        GeoDataFrame for specific country
    """
    warnings.filterwarnings("ignore", message="")    
    country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}
    geofabrik_name = [value[0] for (key, value) in country_dict.items() if key == country_code][0]
    #temp_dir = tempfile.TemporaryDirectory()
    finished_files = list(filter(re.compile(f"{geofabrik_name}-latest.osm.pbf").match, sorted(list(os.listdir(Path(download_dir))))))
    if len(finished_files) == 0:
        pyrosm.data.get_data(dataset = geofabrik_name, update = True, directory = download_dir)  
    #print(os.listdir(download_dir))
    extract = pyrosm.OSM(str(download_dir + f'/{geofabrik_name}-latest.osm.pbf'))
    
    osm_data = gpd.GeoDataFrame({'id': pd.Series(dtype='int64'), 'timestamp': pd.Series(dtype='int64'), 
                                 'version': pd.Series(dtype='int8'), 'osm_type': pd.Series(dtype='object'), 
                                 'highway': pd.Series(dtype='object'), 'boundary': pd.Series(dtype='object'), 
                                 'natural': pd.Series(dtype='object'), 'waterway': pd.Series(dtype='object'),  
                                 'railway': pd.Series(dtype='object'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 
    
    osm_network = extract.get_network(network_type='all')
    if osm_network is not None:
        osm_network["geometry_types"] = osm_network.geometry.geom_type
        osm_network_lines = osm_network.loc[osm_network["geometry_types"].isin(["LineString", "MultiLineString"])]
        osm_network_lines = osm_network_lines[['id','timestamp','version','osm_type','highway','geometry']].to_crs(4326)
        osm_data = osm_data.append(osm_network_lines, ignore_index=True)

    osm_boundaries = extract.get_boundaries(boundary_type='all')
    if osm_boundaries is not None:
        osm_boundaries["geometry_types"] = osm_boundaries.geometry.geom_type
        osm_boundaries_lines = osm_boundaries.loc[osm_boundaries["geometry_types"].isin(["LineString", "MultiLineString"])]
        osm_boundaries_lines = osm_boundaries_lines[['id','timestamp','version','osm_type','boundary','geometry']].to_crs(4326)
        osm_data = osm_data.append(osm_boundaries_lines, ignore_index=True)

    osm_natural = extract.get_natural(custom_filter={'natural': ['coastline','water']}) #{'natural': ['coastline']})
    if osm_natural is not None:
        osm_natural["geometry_types"] = osm_natural.geometry.geom_type
        osm_natural_lines = osm_natural.loc[osm_natural["geometry_types"].isin(["LineString", "MultiLineString"])]
        osm_natural_lines = osm_natural_lines[['id','timestamp','version','osm_type','natural','geometry']].to_crs(4326)
        osm_data = osm_data.append(osm_natural_lines, ignore_index=True)

    osm_waterways = extract.get_data_by_custom_criteria(custom_filter={'waterway': Conf.tags.waterway})
    if osm_waterways is not None:
        osm_waterways["geometry_types"] = osm_waterways.geometry.geom_type
        osm_waterways_lines = osm_waterways.loc[osm_waterways["geometry_types"].isin(["LineString", "MultiLineString"])]
        osm_waterways_lines = osm_waterways_lines[['id','timestamp','version','osm_type','waterway','geometry']].to_crs(4326)
        osm_data = osm_data.append(osm_waterways_lines, ignore_index=True)
    
    osm_railways = extract.get_data_by_custom_criteria(custom_filter={'railway': Conf.tags.railway})
    if osm_railways is not None:
        osm_railways["geometry_types"] = osm_railways.geometry.geom_type
        osm_railways_lines = osm_railways.loc[osm_railways["geometry_types"].isin(["LineString", "MultiLineString"])]
        osm_railways_lines = osm_railways_lines[['id','timestamp','version','osm_type','railway','geometry']].to_crs(4326)
        osm_data = osm_data.append(osm_railways_lines, ignore_index=True)
    
    #temp_dir.cleanup()
    osm_data.fillna('', inplace=True)
    osm_data = osm_data.reset_index(drop=True)
    osm_data['line_id'] = osm_data.groupby(['id']).cumcount()+1
    osm_data['line_id'] = osm_data['id'].astype(str) + '_' + osm_data['line_id'].astype(str)
    osm_data = osm_data[['line_id','id','timestamp','version','osm_type','highway','boundary','natural','waterway','railway','geometry']]
    return osm_data


def gadm_dir_to_path(gadm_dir: Union[str, Path]) -> str:
    """
    For a given country, the GADM directory contains multiple file levels
    so this convenience function just returns the path to the highest
    resolution GADM file within the directory
    Args:
        gadm_dir: directory containing all gadm files for a country
    Returns:
        Path to specific gadm file
    """
    sort_fn = lambda p: int(p.stem.split("_")[-1])
    gadm_dir = Path(gadm_dir)
    files = [p for p in gadm_dir.iterdir() if ".shp" in p.name]
    files.sort(key=sort_fn)
    return files[-1]

def get_gadm_data(country_code: str, download_dir: Union[str, Path]) -> gpd.GeoDataFrame:
    """
    Download most detailed administrative boundaries from GADM.org
    Args:
        country_code: three-letter country code following ISO 3166-1 alpha-3 format
    Returns:
        GeoDataFrame for a specific country
    """
    dest_exist = list(filter(re.compile(f"{country_code}").match, sorted(list(os.listdir(Path(download_dir))))))
    if dest_exist[0] == country_code:
        gadm_url = f"https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_{country_code}_shp.zip"
        #temp_dir = tempfile.TemporaryDirectory()
        results = requests.get(URL(gadm_url))
        z = zipfile.ZipFile(io.BytesIO(results.content))
        dest_folder = Path(download_dir) / country_code
        dest_folder.mkdir(parents=True, exist_ok=True)
        z.extractall(path=dest_folder)
    gadm_file = gadm_dir_to_path(gadm_dir = Path(download_dir) / country_code )
    data = gpd.read_file(gadm_file).to_crs(4326)
    #print(os.listdir(temp_dir.name))
    #temp_dir.cleanup()
    return data

