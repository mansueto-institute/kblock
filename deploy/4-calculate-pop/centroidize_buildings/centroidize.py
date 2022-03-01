
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
    country_building_dirs = [d for d in top_centroid_dir.iterdir() if d.is_dir()]
    top_level_bounds = list(set([file.name.split('_')[1].split('.')[1] for file in top_centroid_dir.rglob("building_centroids*")]))
    top_level_bounds = zip(top_level_bounds, top_centroid_dir*len(top_level_bounds))
    # pass a tuple of the top level building_centroid directory and a string that's of the form countrycode
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


def aggregate_helper(aggregation_tuple):
    logging.getLogger().setLevel(logging.INFO)
    log_file = Path.cwd() / 'log.txt'
    logging.basicConfig(filename=log_file, format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    top_level_bounds = aggregation_tuple[0]
    building_centroid_dir = aggregation_tuple[1]
    logging.info(f'Working on {top_level_bounds}')
    country_centroid_dir = Path(building_centroid_dir) / top_level_bounds.split('.')[0]
    if not country_centroid_parquet_dir.exists():
        country_centroid_parquet_dir.mkdir(parents=True)

    all_building_files = pd.concat([gpd.read_file(centroid_path.resolve()) for centroid_path in country_centroid_dir.rlob("centroid_buildings_"+top_level_bounds+"*")])
    output_path = Path(building_centroid_dir).parent / 'building_centroid_parquets' / top_level_bounds.split('.')[0] / (top_level_bounds.split('.')[1]+'.parquet')
    all_building_files.to_parquet(output_dir)
    logging.info(f'Done working on {country_building_dir}')
