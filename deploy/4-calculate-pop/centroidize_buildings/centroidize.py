
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
	map(aggregate_helper, country_building_dirs)
	centroid_files_list = top_building_dir.rglob("centroid_*")
	for centroid_geojson in centroid_files_list:
		if 'parquet' not in centroid_geojson.name:
			centroid_geojson.unlink()


	
def centroidize_helper(building_file_path: Path):
	building_file = gpd.read_file(building_file_path).set_crs(3395, allow_override=True)
	building_file['building_area'] = building_file.area
	building_file.geometry = building_file.geometry.centroid
	building_file = building_file.set_crs(4326, allow_override=True)
	building_file = building_file[['osm_id', 'gadm', 'gadm_code', 'building_area', 'geometry']]
	output_file_name = building_file_path.parent / ('centroid_' + building_file_path.stem + building_file_path.suffix)
	building_file.to_file(output_file_name, driver='GeoJSON')


def aggregate_helper(country_building_dir: Path):
	all_building_files = pd.concat([gpd.read_file(building_file_path.resolve()) for building_file_path in country_building_dir.iterdir() if 'centroid' in building_file_path.name and 'geojson' in building_file_path.name])
	all_building_files = dgpd.from_geopandas(all_building_files, chunksize=100)
	output_dir = country_building_dir.parent.parent / 'building_centroids' / (country_building_dir.name + '_centroids')
	if not output_dir.exists():
		output_dir.mkdir(parents=True)
	all_building_files.to_parquet(output_dir)