
import pygeos
from typing import Union, Tuple
from pathlib import Path
import pyarrow
import dask_geopandas as dgpd
import geopandas as gpd
import pandas as pd
import logging

def centroidize(top_building_dir: Union[str, Path]):
    # For each building file
    #     calculate areas
    #     generate centroids
    #     Write out geojson with columns osm_id, gadm, gadm_code, area, geometry
    # For each written geojson
    #     Read in file
    # Aggregate 
    log_path = Path.cwd() / 'log.txt'
    logging.basicConfig(filename=Path(log_path), level=logging.INFO)
    logging.info('Started')
    if isinstance(top_building_dir, str):
        top_building_dir = Path(top_building_dir)
    building_files_list = list(top_building_dir.rglob("buildings_*"))
    # logging.info(f'Got building_files_list {building_files_list}')
    # for building_file in building_files_list:
        # logging.info('calling centroidize_helper on {building_file}')
        # centroidize_helper(building_file)
    list(map(centroidize_helper, building_files_list))
    top_centroid_dir = top_building_dir.parent / 'building_centroids'
    # country_building_dirs = [d for d in top_centroid_dir.iterdir() if d.is_dir()]
    # list(map(aggregate_helper, country_building_dirs))
    # centroid_files_list = top_building_dir.rglob("centroid_*")
    # for centroid_geojson in centroid_files_list:
    #     if 'parquet' not in centroid_geojson.name:
    #         centroid_geojson.unlink()


    
def centroidize_helper(building_file_path: Path):
    logging.getLogger().setLevel(logging.INFO)
    log_file = Path.cwd() / 'log.txt'
    logging.basicConfig(filename=log_file, format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f'Working on {building_file_path} in centroidize_helper')
    output_file_parent = building_file_path.parent.parent.parent / 'building_centroids' / building_file_path.name.split('_')[1].split('.')[0]
    if not output_file_parent.exists():
        output_file_parent.mkdir(parents=True)
    output_file_name = output_file_parent / ('centroid_' + building_file_path.stem + building_file_path.suffix)
    if output_file_name.exists():
    	logging.info(f'{output_file_name} exists, returning')
    	return
    building_file = gpd.read_file(building_file_path).set_crs(3395, allow_override=True)
    building_file['building_area'] = building_file.area
    building_file.geometry = building_file.geometry.centroid
    building_file = building_file.set_crs(4326, allow_override=True)
    building_file = building_file[['osm_id', 'gadm', 'gadm_code', 'building_area', 'geometry']]
    logging.info(f'Wrote {output_file_name}')
    building_file.to_file(output_file_name, driver='GeoJSON')


def aggregate_helper(country_building_dir: Path):
    logging.getLogger().setLevel(logging.INFO)
    log_file = Path.cwd() / 'log.txt'
    logging.basicConfig(filename=log_file, format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f'Working on {country_building_dir}')
    all_building_files = pd.concat([gpd.read_file(building_file_path.resolve()) for building_file_path in country_building_dir.iterdir() if 'centroid' in building_file_path.name and 'geojson' in building_file_path.name])
    all_building_files = dgpd.from_geopandas(all_building_files, npartitions=10)
    output_dir = country_building_dir.parent.parent / 'building_centroid_parquets' / (country_building_dir.name + '_centroids')
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    all_building_files.to_parquet(output_dir)
    logging.info(f'Done working on {country_building_dir}')
