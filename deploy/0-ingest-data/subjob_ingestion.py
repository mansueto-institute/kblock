
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj

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
import argparse
import warnings; warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')


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
        download_dir: folder containing downloaded GADM shapefiles
    Returns:
        GeoDataFrame for a specific country
    """
    dest_exist = list(filter(re.compile(f"{country_code}").match, sorted(list(os.listdir(Path(download_dir))))))
    if len(dest_exist) == 0:
        #gadm_url = f"https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_{country_code}_shp.zip"
        gadm_url = f"https://biogeo.ucdavis.edu/data/gadm4.1/shp/gadm41_{country_code}_shp.zip"
        #temp_dir = tempfile.TemporaryDirectory()
        results = requests.get(URL(gadm_url))
        z = zipfile.ZipFile(io.BytesIO(results.content))
        dest_folder = Path(download_dir) / 'shp' / country_code
        dest_folder.mkdir(parents=True, exist_ok=True)
        z.extractall(path=dest_folder)
    gadm_file = gadm_dir_to_path(gadm_dir = Path(download_dir) / 'shp' / country_code)
    data = gpd.read_file(gadm_file).to_crs(4326)
    #print(os.listdir(temp_dir.name))
    #temp_dir.cleanup()
    return data


def get_optimal_gadm_data(country_code: str, optimal_area: float, download_dir: Union[str, Path]) -> gpd.GeoDataFrame:
    """
    Download most detailed administrative boundaries from GADM.org
    Args:
        country_code: three-letter country code following ISO 3166-1 alpha-3 format
        optimal_area: float, average square kilometers for ideal GADM unit level
        download_dir: folder containing downloaded GADM shapefiles
    Returns:
        GeoDataFrame for a specific country
    """
    dest_exist = list(filter(re.compile(f"{country_code}").match, sorted(list(os.listdir(Path(download_dir))))))
    dest_folder = Path(download_dir) / 'shp' / country_code
    if len(dest_exist) == 0:
        #gadm_url = f"https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_{country_code}_shp.zip"
        gadm_url = f"https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_{country_code}_shp.zip"
        #temp_dir = tempfile.TemporaryDirectory()
        results = requests.get(URL(gadm_url))
        z = zipfile.ZipFile(io.BytesIO(results.content))
        dest_folder.mkdir(parents=True, exist_ok=True)
        z.extractall(path=dest_folder)
    gadm_file = gadm_dir_to_path(gadm_dir = dest_folder)
    data = gpd.read_file(gadm_file).to_crs(3395)
    #print(os.listdir(temp_dir.name))
    #temp_dir.cleanup()
    gid_list = list(filter(re.compile("GID_*").match, list(data.columns)))
    data_levels = pd.DataFrame(list(zip(list(data[gid_list].nunique().index), list(data[gid_list].nunique()))), columns =['GID', 'units'])
    data_levels['avg_km2'] = round((data.area.agg({'sum'})*1e-6)[0]/data_levels['units'],1)
    optimal_gid = list(data_levels.iloc[(data_levels['avg_km2']-optimal_area).abs().argsort()[:1]]['GID'])
    shp_list = list(filter(re.compile(".shp").search, sorted(list(os.listdir(dest_folder)))))
    shp = list(filter(re.compile(str('_' + re.findall('\d+', optimal_gid[0])[0])).search, shp_list))[0]
    data = gpd.read_file(dest_folder / shp).to_crs(4326)
    return data


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


def main(gadm_chunk: list, osm_chunk: list, osm_dir: Path, output_dir: Path):

    for country_code in gadm_chunk:
        print(country_code) 

        gadm_output_dir = str(output_dir) + '/gadm'
        osm_output_dir = str(output_dir) + '/osm/parquet'
        Path(gadm_output_dir).mkdir(parents=True, exist_ok=True)
        Path(osm_output_dir).mkdir(parents=True, exist_ok=True)

        gadm_gpd = get_optimal_gadm_data(country_code = country_code, optimal_area = 4000, download_dir = gadm_output_dir) 

        Path(str(gadm_output_dir) + '/geojson').mkdir(parents=True, exist_ok=True) 
        gadm_gpd.to_file(Path(str(gadm_output_dir) + '/geojson') / f'gadm_{country_code}.geojson')

    for country_code in osm_chunk:
        print(country_code) 

        geofabrik_name = code_to_geofabrik(country_code = country_code)
        osm_lines = gpd.read_file(Path(osm_dir) / f'{geofabrik_name}-latest-linestring.geojson')
        osm_lines = osm_lines.rename(columns={'@type':'type','@id':'id','@version':'version','@changeset':'changeset','@timestamp':'timestamp'})
        osm_lines.to_parquet(Path(osm_output_dir) / f'{country_code}-linestring.parquet', compression='snappy')

        osm_polys = gpd.read_file(Path(osm_dir) / f'{geofabrik_name}-latest-polygon.geojson')
        osm_polys = osm_polys.rename(columns={'@type':'type','@id':'id','@version':'version','@changeset':'changeset','@timestamp':'timestamp'})
        osm_polys.to_parquet(Path(osm_output_dir) / f'{country_code}-polygon.parquet', compression='snappy')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Download and build blocks.')
    parser.add_argument('--gadm_chunk', required=False, type=str, dest="gadm_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format")
    parser.add_argument('--osm_chunk', required=False, type=str, dest="osm_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format")
    parser.add_argument('--osm_dir', required=True, type=Path, dest="osm_dir", help="OSM directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Output directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))

