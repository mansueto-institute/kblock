
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
        gadm_url = f"https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_{country_code}_shp.zip"
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


def main(country_chunk: list, output_dir: Path):

    for country_code in country_chunk:
        print(country_code) 
        gadm_gpd = get_gadm_data(country_code = country_code, download_dir = output_dir ) 
        gadm_gpd.to_file(Path(output_dir) / f'gadm_{country_code}.geojson')


def setup(args=None):
    parser = argparse.ArgumentParser(description='Download and build blocks.')
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Output directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))

