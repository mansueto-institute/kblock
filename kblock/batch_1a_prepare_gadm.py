
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

import os
from pathlib import Path
import warnings
import itertools
import re
import time
import argparse
import logging
import math
import dask_geopandas
np.set_printoptions(suppress=True)
warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')
warnings.filterwarnings('ignore', message='.*dropped geometries of different geometry types than.*')

def remove_overlaps(data: gpd.GeoDataFrame, group_column: str, partition_count: int = 10) -> gpd.GeoDataFrame:
    """ 
    Args:
        data: GeoDataFrame, containing delineations, requires CRS WGS 84 EPSG 4326.
        group_column: str, column name with label uniquely identifying each row geometry.
        partition_count: number of partitions to use with dask_geopandas.sjoin(), defaults to 10.
    Returns:
        GeoDataFrame with all overlapping geometries removed. 
        Assigns overlapping sections to the geometry group_column with most area in common or is closest.
    """
    column_list = list(data.columns)

    column_list_no_geo = list(data.columns)
    column_list_no_geo.remove('geometry')

    column_list_id = list(data.columns)
    column_list_id.remove('geometry')
    column_list_id.append('overlap_id')

    data_overlap = dask_geopandas.from_geopandas(data, npartitions = partition_count)
    data_overlap = dask_geopandas.sjoin(left = data_overlap, right = data_overlap, predicate="overlaps")
    data_overlap = data_overlap.compute()

    if data_overlap.shape[0] > 0: 
        print(f'{data_overlap.shape[0]} overlaps.')

        data_overlap = data_overlap.rename(columns={str(group_column + '_left'): group_column})
        data_overlap = data_overlap[[group_column,'geometry']]
        
        overlap_list = data_overlap[group_column].unique()
        data_overlap = data[data[group_column].isin(overlap_list)]
        
        all_intersections = [a.intersection(b) for a, b in list(itertools.combinations(data_overlap['geometry'], 2))]
        data_overlap = pd.concat([data_overlap['geometry'], gpd.GeoSeries(all_intersections).set_crs(4326)])
        data_overlap = list(shapely.ops.polygonize(data_overlap.boundary.unary_union))
        data_overlap = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(data_overlap)}).set_crs(4326).reset_index(drop=True)  
        data_overlap['geometry'] = data_overlap['geometry'].make_valid()
        data_overlap = data_overlap.assign(overlap_id = data_overlap.index.astype(int))
        
        data_overlay = gpd.overlay(df1 = data_overlap, df2 = data, how='intersection', keep_geom_type = True, make_valid = True)
        data_overlay = data_overlay.assign(area = round(data_overlay['geometry'].to_crs(3395).area,0))
        data_overlay['area_rank'] = data_overlay.groupby('overlap_id')['area'].rank(method='first', ascending=False)
        data_overlay = data_overlay[data_overlay[group_column].notnull()]
        data_overlay = data_overlay[data_overlay['area_rank'] == 1]
        data_overlay = data_overlay[column_list_id]
        
        data_sjoin = data_overlap[~data_overlap['overlap_id'].isin(data_overlay['overlap_id'].unique())]
        data_sjoin = gpd.sjoin_nearest(left_df = data_sjoin.to_crs(3395), right_df = data.to_crs(3395), how = 'left').to_crs(4326)
        data_sjoin = data_sjoin.drop_duplicates(subset=['overlap_id'], keep='first')
        data_sjoin = data_sjoin[column_list_id]
        
        data_corrected = pd.concat([data_overlay, data_sjoin], ignore_index=True)
        data_corrected = pd.merge(left = data_overlap, right = data_corrected, how='left', on='overlap_id')
        data_corrected = pd.concat([data_corrected[column_list], data[~data[group_column].isin(overlap_list)]])
        
        data_corrected['geometry'] = data_corrected['geometry'].make_valid()
        data_corrected = data_corrected.dissolve(by=column_list_no_geo, as_index=False)
        data_corrected['geometry'] = data_corrected['geometry'].make_valid()
        
        check = dask_geopandas.from_geopandas(data_corrected, npartitions = partition_count)
        check = dask_geopandas.sjoin(left = check, right = check, predicate="overlaps")
        check = check.compute()
        if check.shape[0] > 0: 
            warnings.warn(f'Unable to remove all overlaps. {check.shape[0]} overlaps remain.')
        else:
            print('All overlaps resolved.')
    else:
        print('No overlaps found.')
        data_corrected = data

    return data_corrected 

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

        # Ensure no overlaps with surrounding countries
        # other_countries = all_gadm_gpd[~all_gadm_gpd['GID_0'].isin([country_code])].to_crs(4326)
        # gadm_country = gpd.overlay(df1 = gadm_country, df2 = other_countries, how = 'difference', keep_geom_type = True, make_valid = True)
        # gadm_country['geometry'] = gadm_country['geometry'].make_valid()

        # Remove non-polygons from GeometryCollection GADMs
        if not all(x in ['Polygon','MultiPolygon'] for x in gadm_country['geometry'].geom_type.unique()):
            gadm_country = gadm_country.explode(index_parts=False)
            gadm_country = gadm_country[gadm_country['geometry'].geom_type == 'Polygon']
            gadm_country = gadm_country[round(gadm_country['geometry'].to_crs(3395).area,0) > 0]
            gadm_country = gadm_country.dissolve(by='gadm_code', as_index=False)
            gadm_country['geometry'] = gadm_country['geometry'].make_valid()

        # Remove overlapping geometries
        gadm_country = remove_overlaps(data = gadm_country, group_column = 'gadm_code')

        # Write countries
        gadm_country.to_parquet(Path(gadm_output_dir) / f'gadm_{country_code}.parquet', compression='snappy')
        gadm_country.to_file(Path(gpkg_dir) / f'gadm_{country_code}.gpkg', driver="GPKG")

        # Check for overlaps and write to file if they exist
        check = dask_geopandas.from_geopandas(gadm_country, npartitions = 10)
        check = dask_geopandas.sjoin(left = check, right = check, predicate="overlaps")
        check = check.compute()
        if check.shape[0] > 0: 
            logging.info(f"Number of overlaps: {check[0].shape}")
            check.to_file(Path(overlaps_dir) / f'overlaps_{country_code}.gpkg', driver="GPKG")

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
