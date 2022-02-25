
import pygeos
from typing import Union, Tuple
from pathlib import Path
import pyarrow
import dask_geopandas as dgpd
import geopandas as gpd
import pandas as pd

def centroidize(top_building_dir: Union[str, Path]):
	# For each building file
	# 	calculate areas
	# 	generate centroids
	# 	Write out geojson with columns osm_id, gadm, gadm_code, area, geometry
	# For each written geojson
	# 	Read in file
	# Aggregate 

	if isinstance(top_building_dir, str):
		top_building_dir = Path(top_building_dir)
	building_files_list = top_building_dir.rglob("buildings*")
	map(centroidize_helper, list(building_files_list))
	country_building_dirs = [d for d in top_building_dir.iterdir() if d.is_dir()]
	for country_building_dir in country_building_dirs:
		aggregate_helper(country_building_dir)
	map(aggregate_helper, country_building_dirs)


	
def centroidize_helper(building_file_path: Path):
	building_file = gpd.read_file(building_file_path).set_crs(3395, allow_override=True)
	building_file = dgpd.from_geopandas(building_file, chunksize=100)
	building_file['building_area'] = building_file.area.compute()
	building_file.geometry = building_file.geometry.centroid
	building_file = building_file.set_crs(4326, allow_override=True)
	building_file = building_file[['osm_id', 'gadm', 'gadm_code', 'building_area', 'geometry']]
	output_file_name = building_file_path.parent / ('centroid_' + building_file_path.stem + building_file_path.suffix)
	building_file.to_file(output_file_name, driver='GeoJSON')


def aggregate_helper(country_building_dir: Path):
	# print(country_building_dir.iterd)
	all_building_files = pd.concat([gpd.read_file(building_file_path.resolve()) for building_file_path in country_building_dir.iterdir() if 'centroid' in building_file_path.name])
	all_building_files = dgpd.from_geopandas(all_building_files, chunksize=100)
	all_building_files.to_parquet(country_building_dir / (country_building_dir.name+'_centroids.parquet'))