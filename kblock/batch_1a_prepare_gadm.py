
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

import os
from pathlib import Path
import warnings
import re
import time
import argparse
import logging
import math
import dask_geopandas
np.set_printoptions(suppress=True)
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')
warnings.filterwarnings('ignore', message='.*dropped geometries of different geometry types than.*')

def main(log_file: Path, country_chunk: list, gadm_dir: Path, daylight_dir: Path, osm_dir: Path, output_dir: Path):

    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    
    # Make directory
    gadm_output_dir = str(output_dir) + '/parquet'
    Path(gadm_output_dir).mkdir(parents=True, exist_ok=True)

    gpkg_dir = str(output_dir) + '/gpkg'
    Path(gpkg_dir).mkdir(parents=True, exist_ok=True)

    overlaps_dir = str(output_dir) + '/overlaps'
    Path(overlaps_dir).mkdir(parents=True, exist_ok=True)

    combined_dir = str(output_dir) + '/combined'
    Path(combined_dir).mkdir(parents=True, exist_ok=True)

    # Check GADM files in gadm_dir
    gadm_input_list = list(filter(re.compile("gadm_").match, sorted(list(os.listdir(Path(gadm_dir))))))
    gadm_input_list = [(re.sub('gadm_', '', re.sub('.geojson', '', i))) for i in gadm_input_list] 
    in_chunk_not_in_gadm_inputs = [x for x in country_chunk if x not in set(gadm_input_list)]
    if len(in_chunk_not_in_gadm_inputs) > 0:
        raise ValueError(f'GADM input data does not exist for {in_chunk_not_in_gadm_inputs} in country_chunk argument.')

    # Consolidate GADM data into one file
    all_gadm_gpd = gpd.GeoDataFrame({'GID_0': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 
    for country_code in gadm_input_list: 
        gadm_gpd = gpd.read_file(Path(gadm_dir) / f'gadm_{country_code}.geojson')
        gadm_gpd = gadm_gpd[['GID_0', 'geometry']]
        gadm_gpd = gadm_gpd.dissolve(by='GID_0', as_index=False)
        all_gadm_gpd = pd.concat([all_gadm_gpd, gadm_gpd], ignore_index=True)        

    # Check for completed countries in output_dir
    output_file_list = list(filter(re.compile("gadm_").match, sorted(list(os.listdir(Path(gadm_output_dir))))))
    output_country_list = [(re.sub('gadm_', '', re.sub('.parquet', '', i))) for i in output_file_list] 
    logging.info(f"Finished countries: {output_country_list}")

    # Remove completed countries form country_chunk argument
    if country_chunk: 
        country_list = [x for x in country_chunk if x not in set(output_country_list)]
        if not country_list: 
            raise ValueError('All countries in country_chunk argument already finished.')
    else: 
        raise ValueError('Empty country_chunk arg.')

    logging.info(f"Countries to process: {country_list}")
    logging.info(f"Preprocess GADM geometries")

    # Read in data from daylight_dir
    daylight_inland_land = gpd.read_file(Path(daylight_dir) / 'land_polygons.shp')
    daylight_coastal_water = gpd.read_file(Path(daylight_dir) / 'water_polygons.shp')

    daylight_inland_land = daylight_inland_land.to_crs(4326)
    daylight_coastal_water = daylight_coastal_water.to_crs(4326)

    # Iterate through country codes
    for country_code in country_list: 
        print(country_code)
        t0 = time.time()
        logging.info(f"Country: {country_code}")
        logging.info(f"Removing internal water features")

        # Remove internal water features from administrative boundaries
        osm_inland_water = gpd.read_parquet(Path(osm_dir) / f'{country_code}-polygon.parquet')
        osm_inland_water = osm_inland_water[~osm_inland_water['natural'].isin(['coastline'])]
        gadm_country = gpd.read_file(Path(gadm_dir) / f'gadm_{country_code}.geojson')
        gadm_country['geometry'] = gadm_country['geometry'].make_valid()
        # gadm_country['geometry'] = gadm_country['geometry'].to_crs(3395).buffer(0.001).buffer(-0.001).to_crs(4326)
        gadm_country = gpd.overlay(df1 = gadm_country, df2 = daylight_coastal_water, how = 'difference', keep_geom_type = True, make_valid = True)
        gadm_country = gpd.overlay(df1 = gadm_country, df2 = osm_inland_water, how = 'difference', keep_geom_type = True, make_valid = True)
        gadm_col = max(list(filter(re.compile("GID_*").match, list(gadm_country.columns))))
        gadm_country = gadm_country.rename(columns={gadm_col: "gadm_code"})
        gadm_country['country_code'] = country_code
        gadm_country = gadm_country[['gadm_code', 'country_code', 'geometry']]
        gadm_country['geometry'] = gadm_country['geometry'].make_valid()

        # Fix jagged coastlines
        gadm_buffer = pygeos.union_all(pygeos.buffer(pygeos.get_parts(pygeos.from_shapely(gadm_country['geometry'].to_crs(3395))),3000))
        coastal_buffer = pygeos.union_all(pygeos.intersection(pygeos.from_shapely(daylight_coastal_water['geometry'].to_crs(3395)), pygeos.convex_hull(gadm_buffer))) 
        if pygeos.is_empty(coastal_buffer) == False:
            logging.info(f"Buffering coastline")
            fix_envelope = pygeos.convex_hull(pygeos.intersection(coastal_buffer, pygeos.union_all(pygeos.get_parts(pygeos.from_shapely(gadm_country['geometry'].to_crs(3395))))))
            land_buffer = pygeos.intersection(pygeos.from_shapely(daylight_inland_land['geometry'].to_crs(3395)), pygeos.buffer(pygeos.union_all(coastal_buffer),3000))
            land_buffer = pygeos.intersection(land_buffer, fix_envelope)

            if all(pygeos.is_empty(land_buffer)) == False:

                land_buffer = pygeos.to_shapely(pygeos.get_parts(land_buffer[~pygeos.is_empty(land_buffer)]))
                coast_fixes = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(land_buffer).set_crs(3395).to_crs(4326)}).set_crs(4326)
                coast_fixes = coast_fixes.explode(index_parts=False)
                coast_fixes = coast_fixes[coast_fixes.geom_type == "Polygon"]
                coast_fixes = gpd.overlay(df1 = coast_fixes, df2 = gadm_country, how = 'difference', keep_geom_type = True, make_valid = True)

                coast_fixes['geometry'] = coast_fixes['geometry'].make_valid()
                coast_fixes = coast_fixes.explode(index_parts=False)
                coast_fixes = coast_fixes[coast_fixes.geom_type == "Polygon"]

                osm_inland_water['geometry'] = osm_inland_water['geometry'].make_valid()
                osm_inland_water = osm_inland_water.explode(index_parts=False)
                osm_inland_water = osm_inland_water[osm_inland_water.geom_type == "Polygon"]
                coast_fixes = gpd.overlay(df1 = coast_fixes, df2 = osm_inland_water, how = 'difference', keep_geom_type = True, make_valid = True)
                
                coast_fixes['geometry'] = coast_fixes['geometry'].make_valid()
                coast_fixes = coast_fixes.explode(index_parts=False)
                coast_fixes = coast_fixes[coast_fixes.geom_type == "Polygon"]

                daylight_coastal_water['geometry'] = daylight_coastal_water['geometry'].make_valid()
                daylight_coastal_water = daylight_coastal_water.explode(index_parts=False)
                daylight_coastal_water = daylight_coastal_water[daylight_coastal_water.geom_type == "Polygon"]
                coast_fixes = gpd.overlay(df1 = coast_fixes, df2 = daylight_coastal_water, how = 'difference', keep_geom_type = True, make_valid = True)

                other_countries = all_gadm_gpd[~all_gadm_gpd['GID_0'].isin([country_code])].to_crs(4326)
                coast_fixes = gpd.overlay(df1 = coast_fixes, df2 = other_countries, how = 'difference', keep_geom_type = True, make_valid = True)

                coast_fixes = coast_fixes.explode(ignore_index = True)
                
                area_selection = pygeos.from_shapely(all_gadm_gpd[all_gadm_gpd['GID_0'].isin([country_code])].to_crs(3395).buffer(500).to_crs(4326))
                coast_fixes = pygeos.from_shapely(coast_fixes['geometry'])
                coast_fixes = pygeos.to_shapely(coast_fixes[pygeos.intersects(coast_fixes, area_selection)])
                coast_fixes = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(coast_fixes)}).set_crs(4326)

                #coast_fixes['geometry'] = coast_fixes['geometry'].to_crs(3395).buffer(0.001).buffer(-0.001).to_crs(4326)
                coast_fixes = coast_fixes[coast_fixes.geom_type == "Polygon"]
                coast_fixes = coast_fixes.assign(fix_id = [str(x) for x in list(coast_fixes.index)])
                coast_fixes = gpd.sjoin_nearest(left_df = coast_fixes.to_crs(3395), right_df = gadm_country.to_crs(3395), how = 'left').to_crs(4326)
                coast_fixes = coast_fixes.drop_duplicates(subset=['fix_id'], keep='first')
                
                gadm_country = pd.concat([gadm_country, coast_fixes[['gadm_code','country_code','geometry']]], ignore_index=True)
                gadm_country['geometry'] = gadm_country['geometry'].make_valid()
                gadm_country = gadm_country[round(gadm_country['geometry'].to_crs(3395).area,0) > 0]
                gadm_country = gadm_country.dissolve(by=['gadm_code', 'country_code'], as_index=False)
                if gadm_country['gadm_code'].isnull().values.any(): 
                    raise TypeError(f"{country_code}: GADM column contains null.")
                gadm_country['geometry'] = gadm_country['geometry'].make_valid()

        # Remove non-polygons from GeometryCollection GADMs
        if not all(x in ['Polygon','MultiPolygon'] for x in gadm_country['geometry'].geom_type.unique()):
            gadm_country = gadm_country.explode(index_parts=False)
            gadm_country = gadm_country[gadm_country['geometry'].geom_type == 'Polygon']
            gadm_country = gadm_country[round(gadm_country['geometry'].to_crs(3395).area,0) > 0]
            gadm_country = gadm_country.dissolve(by='gadm_code', as_index=False)
            gadm_country['geometry'] = gadm_country['geometry'].make_valid()
            
        # Write countries
        gadm_country.to_parquet(Path(gadm_output_dir) / f'gadm_{country_code}.parquet', compression='snappy')
        gadm_country.to_file(Path(gpkg_dir) / f'gadm_{country_code}.gpkg', driver="GPKG")

        # Check for overlaps
        gadm_country_overlaps = gpd.overlay(df1 = gadm_country, df2 = gadm_country, how = 'intersection', keep_geom_type = True, make_valid = True)
        gadm_country_overlaps = gadm_country_overlaps[gadm_country_overlaps['gadm_code_1'] != gadm_country_overlaps['gadm_code_2']]
        if gadm_country_overlaps.shape[0] > 0:
            gadm_country_overlaps = gadm_country_overlaps.assign(area = round(gadm_country_overlaps.to_crs(3395).area,0))
            overlap_count = gadm_country_overlaps[gadm_country_overlaps['area'] > 100].shape[0]
            if overlap_count > 0:
                gadm_country_overlaps['gadm_code'] = gadm_country_overlaps[['gadm_code_1', 'gadm_code_2']].bfill(axis=1).iloc[:, 0]
                gadm_country_overlaps['country_code'] = gadm_country_overlaps[['country_code_1', 'country_code_2']].bfill(axis=1).iloc[:, 0]
                logging.info(f"Overlaps over 100m2: {gadm_country_overlaps[gadm_country_overlaps['area'] > 100].shape[0]}")
                logging.info(f"{gadm_country_overlaps['gadm_code'].unique()}")
                gadm_country_overlaps = gadm_country_overlaps.dissolve(by=['gadm_code', 'country_code'], as_index=False)
                gadm_country_overlaps = gadm_country_overlaps[['gadm_code','country_code','geometry']].reset_index(drop = True)
                gadm_country_overlaps['geometry'] = gadm_country_overlaps['geometry'].make_valid()
                gadm_country_overlaps.to_file(Path(overlaps_dir) / f'overlaps_{country_code}.gpkg', driver="GPKG")

        t1 = time.time()
        logging.info(f"Finished {country_code}: {gadm_country.shape} {str(round(t1-t0,3)/60)} minutes")

    # Consolidate GADM data into one file
    gadm_output_list = list(filter(re.compile("gadm_").match, sorted(list(os.listdir(Path(gadm_output_dir))))))
    gadm_output_list = [(re.sub('gadm_', '', re.sub('.parquet', '', i))) for i in gadm_output_list] 
    gadm_combo = gpd.GeoDataFrame({'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326) 
    for country_code in gadm_output_list: 
        gadm_clean = gpd.read_parquet(Path(gadm_output_dir) / f'gadm_{country_code}.parquet')
        gadm_combo = pd.concat([gadm_combo, gadm_clean], ignore_index=True)   
    
    gadm_combo.to_parquet(Path(combined_dir) / f'all_gadm.parquet', compression='snappy')
    gadm_combo.to_file(Path(combined_dir) / f'all_gadm.gpkg', driver="GPKG")

    logging.info(f"Finished")

def setup(args=None):
    parser = argparse.ArgumentParser(description='Clean up and combine GADM delineations.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file.") 
    parser.add_argument('--country_chunk', required=True, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format.")
    parser.add_argument('--gadm_dir', required=True, type=Path, dest="gadm_dir", help="Path to GADM geojson directory.")
    parser.add_argument('--daylight_dir', required=True, type=Path, dest="daylight_dir", help="Path to Daylight coastline directory.")
    parser.add_argument('--osm_dir', required=True, type=Path, dest="osm_dir", help="Path to OSM parquet directory.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to top level output directory. Creates gadm/parquet subdirectory.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))
