import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
gpd.options.use_pygeos = True
import shapely

from typing import Union
from pathlib import Path
import argparse
import pyarrow
import dask_geopandas as dgpd
import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from ghsl import make_buffered_GHSL

def main(log_file: Path, ghsl_file: Path, country_chunk: str, blocks_dir: Path, output_dir: Path, ancillary_dir: Path) -> None:
	ghsl_df = dgpd.read_parquet(ghsl_file, 
								filters=[('Geographical Region', '==', 'Western Africa'), 
										 ('Geographical Region', '==', 'Northern Africa'),
										 ('Geographical Region', '==', 'Eastern Africa'),
										 ('Geographical Region', '==', 'Southern Africa')])
	ghsl_df = ghsl_df.compute()
	# For testing with just freetown area, the cx is larger than Freetown
	# freetown = ghsl_df[ghsl_df['Name of the Urban Centre'] == 'Freetown']
	# freetown = ghsl_df.cx[-1480000:-1430000, 910000:944000]
	# ghsl_df = freetown
	for country in country_chunk:
		blocks = dgpd.read_parquet(blocks_dir, filters=[('country_code', '==', country)]).compute()
		blocks = blocks.to_crs(3395)
		ghsl_df_sindex = ghsl_df.sindex
		index_bulk = ghsl_df_sindex.query_bulk(blocks['geometry'], predicate='within')
		# print(index_bulk)
		ghsl_blocks_map = pd.DataFrame({'index_blocks': index_bulk[0], 'index_ghsl': index_bulk[1]})
		ghsl_blocks_map = ghsl_blocks_map.merge(blocks, how='left', left_on='index_blocks', right_index=True)
		# print(ghsl_blocks_map.head())
		ghsl_blocks_map = ghsl_blocks_map.merge(ghsl_df, how='left', left_on='index_ghsl', right_index=True)
		ghsl_blocks_map = ghsl_blocks_map.rename(columns={'geometry_x': 'blocks_geometry', 'geometry_y': 'ghsl_geometry'})
		data = gpd.GeoDataFrame(ghsl_blocks_map, geometry='ghsl_geometry', crs=3395)
		data = data.drop(columns=['index_blocks', 'index_ghsl'])

		# 10 km buffer
		metro_areas = gpd.GeoDataFrame(data[['ghsl_geometry']], geometry='ghsl_geometry', crs=3395)
		metro_areas.geometry = metro_areas.buffer(10000)

		# To get the metros
		metro_areas = metro_areas.dissolve().explode().reset_index(drop=True)
		metro_areas = metro_areas.rename_geometry('metro_geometry')

		metro_areas_sindex = metro_areas.sindex
		index_bulk = metro_areas_sindex.query_bulk(data['ghsl_geometry'], predicate='intersects')
		metros_blocks_map = pd.DataFrame({'index_blocks': index_bulk[0], 'index_metros': index_bulk[1]})
		metros_blocks_map = metros_blocks_map.merge(data, how='left', left_on='index_blocks', right_index=True)

		metros_blocks_map = metros_blocks_map.merge(metro_areas, how='left', left_on='index_metros', right_index=True)
		metros_blocks_map = gpd.GeoDataFrame(metros_blocks_map, geometry='ghsl_geometry', crs=3395)
		metros_blocks_map['GHSL Area'] = metros_blocks_map.area
		ghsl_area_to_metro_dict = {}
		for metro_area_id in list(set(metros_blocks_map['index_metros'])):
			metro_area_ghsl_areas = metros_blocks_map[metros_blocks_map['index_metros']==metro_area_id]
			metro_area_ghsl_areas = metro_area_ghsl_areas.sort_values(by='GHSL Area', ascending=False)
			biggest_ghsl_area_name_in_metro = metro_area_ghsl_areas.iloc[0]['Name of the Urban Centre']
			ghsl_area_to_metro_dict.update({index_val:biggest_ghsl_area_name_in_metro for index_val in metro_area_ghsl_areas.index.tolist()})
		metros_blocks_map['urban_agglomeration'] = list(map(lambda x: ghsl_area_to_metro_dict[x], metros_blocks_map.index.tolist()))
		metros_blocks_map = metros_blocks_map.rename(columns={'Name of the Urban Centre': 'urban_center'})
		metros_blocks_map = metros_blocks_map.merge(blocks, how='outer', on='block_id')
		metros_blocks_map['urban_agglomeration'] = ["missing" if isinstance(x, float) else x for x in metros_blocks_map['urban_agglomeration']]
		metros_blocks_map['urban_center'] = ["missing" if isinstance(x, float) else x for x in metros_blocks_map['urban_center']]


		data = metros_blocks_map[['block_id', 'gadm_code', 'country_code', 'urban_center', 'urban_agglomeration']]
		assert (len(data['block_id']) == len(list(set(data['block_id']))))
		data_filename = output_dir / (country+'_ghsl.parquet')
		data.to_parquet(data_filename)

		ancillary_data = metros_blocks_map.drop(columns=['index_blocks', 'index_metros', 'block_id', 'country_code', 'blocks_geometry',
			                                     'Unique ID', 'ghsl_geometry', 'Area', 'metro_geometry', 'GHSL Area', 'urban_agglomeration'])
		ancillary_data = ancillary_data.set_index(ancillary_data['urban_center'], drop=True)
		ancillary_data_filename = ancillary_dir / (country+'_ancillary.parquet')
		ancillary_data.to_parquet(ancillary_data_filename)


def setup(args=None):	
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file") 
	parser.add_argument('--ghsl_file', required=True, type=str, dest='ghsl_file', help="Path to GHSL parquet file")
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format")
	parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Blocks directory")
	parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Output directory")
	parser.add_argument('--ancillary_dir', required=True, type=Path, dest='ancillary_dir', help="Directory for ancillary data")
	return parser.parse_args(args)

if __name__ == "__main__":
	main(**vars(setup()))