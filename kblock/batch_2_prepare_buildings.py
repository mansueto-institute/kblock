
import geopandas as gpd
import pandas as pd
import logging
from pathlib import Path
import zipfile
import os
import shutil
import argparse
import time
import re 
import pygeohash

import dask_geopandas 
import psutil
import math

gpd.options.use_pygeos = True


def transform_buildings(gpd_buildings: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Convert pandas.DataFrame containing bounding box columns to a geopandas.GeoDataFrame
    Args:
        gpd_buildings: pandas.DataFrame returned from merge_pixels() function
    Returns:
        pandas.DataFrame converted to geopandas.GeoDataFrame with pixels formatted as geometries
    """
    assert gpd_buildings.crs in ['epsg:3395','epsg:4326'], "gpd_buildings is not epsg:4326 or epsg:3395."
    if gpd_buildings.crs == 'epsg:4326': gpd_buildings = gpd_buildings.to_crs(epsg=3395)

    gpd_buildings['building_area'] = gpd_buildings.area
    gpd_buildings['geometry'] = gpd_buildings['geometry'].centroid
    gpd_buildings = gpd_buildings.to_crs(4326)
    return gpd_buildings

def sjoin_buildings(block_data: gpd.GeoDataFrame, bldg_id_col: str, bldg_data: gpd.GeoDataFrame, use_dask_geopandas = False) -> gpd.GeoDataFrame:
    """
    Joins GeoDataFrame of enclosed block geometries to a GeoDataFrame containing building footprints, uses an inner join with 'intersects' predicate
    Args:
        block_data: GeoDataFrame, output returned from build_blocks() function, requires CRS WGS 84 EPSG 4326 or 3395
        bldg_id_col: str, unique ID column in bldg_data
        bldg_data: GeoDataFrame, containing building geometries with building_area column, requires CRS WGS 84 EPSG 4326 or 3395 (accepts 'Polygon' or 'Point' geometries)
        use_dask_geopandas:, boolean, defaults to False, use dask_geopandas to peform join operation
    Returns:
        GeoDataFrame with building geometries mapped to bldg_id_col,'building_area','block_id','gadm_code','country_code'.
        Geometry projected in CRS WGS 84 EPSG 4326. 
    """
    assert block_data.crs in ['epsg:3395','epsg:4326'], "block_data is not epsg:4326 or epsg:3395."
    assert bldg_data.crs in ['epsg:3395','epsg:4326'], "bldg_data is not epsg:4326 or epsg:3395."

    if block_data.crs == 'epsg:4326': block_data = block_data.to_crs(epsg=3395)
    if bldg_data.crs == 'epsg:4326': bldg_data = bldg_data.to_crs(epsg=3395)

    assert all(x in ['Point','Polygon'] for x in bldg_data['geometry'].geom_type.unique()), 'bldg_data must contain only Point or Polygon geometries.'

    block_data = block_data.reset_index(drop=True)
    bldg_data = bldg_data.reset_index(drop=True)

    if bldg_data['geometry'].geom_type.unique().all() != 'Point':
        bldg_points = bldg_data.copy()
        bldg_points['geometry'] = bldg_points.centroid
    else:
        bldg_points = bldg_data.copy()

    if use_dask_geopandas == False:

        bldg_index = bldg_points.sindex
        index_bulk = bldg_index.query_bulk(block_data['geometry'], predicate="intersects")  
        blocks_buildings_map = pd.DataFrame({'index_blocks': index_bulk[0], 'index_buildings': index_bulk[1]})
        blocks_buildings_map = blocks_buildings_map.merge(bldg_data, how = 'left', left_on='index_buildings', right_index=True)
        blocks_buildings_map = blocks_buildings_map.merge(block_data[list(block_data.columns.difference(['geometry']))], left_on='index_blocks', right_index=True)

        column_list = list(bldg_data.columns.difference(['geometry'])) + list(block_data.columns.difference(['geometry'])) + ['geometry']
        data = gpd.GeoDataFrame(blocks_buildings_map[column_list]).set_crs(epsg=3395)

    else:

        file_size = bldg_points.memory_usage(deep=True).sum() / (1024.0 ** 3)
        mem_avail = psutil.virtual_memory().available / (1024.0 ** 3)
        min_partition_setting = math.ceil(file_size/(mem_avail*.1))
        if min_partition_setting < 4: min_partition_setting = 4

        bldg_points_dask = dask_geopandas.from_geopandas(bldg_points, npartitions = min_partition_setting)
        blocks_buildings_join = dask_geopandas.sjoin(left = bldg_points_dask, right = block_data, predicate="intersects")
        data = blocks_buildings_join.compute()

    data = data.to_crs(epsg=4326)

    return data


def main(log_file: Path, country_chunk: list, codes_file: Path, progress_file: Path, blocks_dir: Path, buildings_dir: Path, output_dir: Path):

    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info('Started')

    # Make directories
    output_folder = Path(output_dir) / 'buildings' / 'polygons'
    output_folder.mkdir(parents=True, exist_ok=True)
    output_folder = Path(output_dir) / 'buildings' / 'points'
    output_folder.mkdir(parents=True, exist_ok=True)

    # Read in table with file names and country codes
    if os.path.isfile(codes_file):
        codes = pd.read_csv(codes_file)
    else: 
        # Default DigitizeAfrica Ecopia file names
        codes = pd.DataFrame.from_dict({
            'file_name': ['africa_angola_building_32732.dl.zip', 'africa_angola_building_32733.dl.zip', 'africa_angola_building_32734.dl.zip', 'africa_angola_building_32735.dl.zip', 'africa_benin_building_32631.dl.zip', 'africa_botswana_building_32734.dl.zip', 'africa_botswana_building_32735.dl.zip', 'africa_burkina_faso_building_32630.dl.zip', 'africa_burkina_faso_building_32631.dl.zip', 'africa_burundi_building_32735.dl.zip', 'africa_burundi_building_32736.dl.zip', 'africa_caboverde_building_32626.dl.zip', 'africa_caboverde_building_32627.dl.zip', 'africa_cameroon_building_32632.dl.zip', 'africa_cameroon_building_32633.dl.zip', 'africa_car_building_32633.dl.zip', 'africa_car_building_32634.dl.zip', 'africa_car_building_32635.dl.zip', 'africa_chad_building_32633.dl.zip', 'africa_chad_building_32634.dl.zip', 'africa_comoros_building_32738.dl.zip', 'africa_congo_building_32633.dl.zip', 'africa_congo_building_32634.dl.zip', 'africa_congo_building_32732.dl.zip', 'africa_congo_building_32733.dl.zip', 'africa_cote_divoire_building_32629.dl.zip', 'africa_cote_divoire_building_32630.dl.zip', 'africa_dem_rep_congo_p1_building_32733.dl.zip', 'africa_dem_rep_congo_p1_building_32734.dl.zip', 'africa_dem_rep_congo_p2_building_32633.dl.zip', 'africa_dem_rep_congo_p2_building_32634.dl.zip', 'africa_dem_rep_congo_p2_building_32635.dl.zip', 'africa_dem_rep_congo_p2_building_32636.dl.zip', 'africa_dem_rep_congo_p2_building_32733.dl.zip', 'africa_dem_rep_congo_p2_building_32734.dl.zip', 'africa_dem_rep_congo_p2_building_32735.dl.zip', 'africa_dem_rep_congo_p2_building_32736.dl.zip', 'africa_djibouti_building_32637.dl.zip', 'africa_djibouti_building_32638.dl.zip', 'africa_equatorial_guinea_building_32632.dl.zip', 'africa_equatorial_guinea_building_32731.dl.zip', 'africa_eritrea_building_32637.dl.zip', 'africa_eritrea_building_32638.dl.zip', 'africa_ethiopia_building_1_32637.dl.zip', 'africa_ethiopia_building_21_32637.dl.zip', 'africa_ethiopia_building_22_32637.dl.zip', 'africa_ethiopia_building_32636.dl.zip', 'africa_ethiopia_building_32638.dl.zip', 'africa_ethiopia_building_3_32637.dl.zip', 'africa_gabon_building_32632.dl.zip', 'africa_gabon_building_32633.dl.zip', 'africa_gabon_building_32732.dl.zip', 'africa_gabon_building_32733.dl.zip', 'africa_gambia_building_32628.dl.zip', 'africa_ghana_building_32630.dl.zip', 'africa_ghana_building_32631.dl.zip', 'africa_guineabissau_building_32628.dl.zip', 'africa_guinea_building_32628.dl.zip', 'africa_guinea_building_32629.dl.zip', 'africa_kenya_building_32636.dl.zip', 'africa_kenya_building_32637.dl.zip', 'africa_kenya_building_32736.dl.zip', 'africa_kenya_building_32737.dl.zip', 'africa_lesotho_building_32735.dl.zip', 'africa_liberia_building_32629.dl.zip', 'africa_madagascar_building_32738.dl.zip', 'africa_madagascar_building_32739.dl.zip', 'africa_malawi_building_32736.dl.zip', 'africa_mali_building_32628.dl.zip', 'africa_mali_building_32629.dl.zip', 'africa_mali_building_32630.dl.zip', 'africa_mali_building_32631.dl.zip', 'africa_mauritania_building_32628.dl.zip', 'africa_mauritania_building_32629_v2.dl.zip', 'africa_mauritania_building_32630_v2.dl.zip', 'africa_mauritius_building_32740.dl.zip', 'africa_mauritius_building_32741.dl.zip', 'africa_mozambique_building_32736.dl.zip', 'africa_mozambique_building_32737.dl.zip', 'africa_namibia_building_32732.dl.zip', 'africa_namibia_building_32733.dl.zip', 'africa_namibia_building_32734.dl.zip', 'africa_namibia_building_32735.dl.zip', 'africa_niger_building_32631.dl.zip', 'africa_niger_building_32632.dl.zip', 'africa_niger_building_32633.dl.zip', 'africa_nigeria_p1_building_32631.dl.zip', 'africa_nigeria_p1_building_32632.dl.zip', 'africa_nigeria_p2_building_32631.dl.zip', 'africa_nigeria_p2_building_32632_1.dl.zip', 'africa_nigeria_p2_building_32632_2.dl.zip', 'africa_nigeria_p2_building_32632_3.dl.zip', 'africa_nigeria_p2_building_32633.dl.zip', 'africa_reunion_building_32740.dl.zip', 'africa_rwanda_building_32735.dl.zip', 'africa_rwanda_building_32736.dl.zip', 'africa_saotomeandprincipe_building_32632.dl.zip', 'africa_saotomeandprincipe_building_32732.dl.zip', 'africa_senegal_building_32628.dl.zip', 'africa_senegal_building_32629.dl.zip', 'africa_seychelles_building_32738.dl.zip', 'africa_seychelles_building_32739.dl.zip', 'africa_seychelles_building_32740.dl.zip', 'africa_sierra_leone_building_32628.dl.zip', 'africa_sierra_leone_building_32629.dl.zip', 'africa_somalia_building_32637.dl.zip', 'africa_somalia_building_32638_v2.dl.zip', 'africa_somalia_building_32639.dl.zip', 'africa_somalia_building_32737.dl.zip', 'africa_somalia_building_32738.dl.zip', 'africa_south_africa_building_32733.dl.zip', 'africa_south_africa_building_32734.dl.zip', 'africa_south_africa_building_32735_1.dl.zip', 'africa_south_africa_building_32735_2.dl.zip', 'africa_south_africa_building_32736.dl.zip', 'africa_south_sudan_building_32635.dl.zip', 'africa_south_sudan_building_32636.dl.zip', 'africa_sudan_building_32634.dl.zip', 'africa_sudan_building_32635.dl.zip', 'africa_sudan_building_32636.dl.zip', 'africa_sudan_building_32637.dl.zip', 'africa_swaziland_building_32736.dl.zip', 'africa_tanzania_building_32735.dl.zip', 'africa_tanzania_building_32736.dl.zip', 'africa_tanzania_building_32737.dl.zip', 'africa_togo_building_32630.dl.zip', 'africa_togo_building_32631.dl.zip', 'africa_uganda_building_32635_v2.dl.zip', 'africa_uganda_building_32636_v2.dl.zip', 'africa_uganda_building_32735_v2.dl.zip', 'africa_uganda_building_32736_v2.dl.zip', 'africa_western_sahara_building_32628_v2.dl.zip', 'africa_western_sahara_building_32629_v2.dl.zip', 'africa_zambia_building_32734.dl.zip', 'africa_zambia_building_32735.dl.zip', 'africa_zambia_building_32736.dl.zip', 'africa_zimbabwe_building_32735_v2.dl.zip', 'africa_zimbabwe_building_32736.dl.zip'],
            'country_code': ['AGO', 'AGO', 'AGO', 'AGO', 'BEN', 'BWA', 'BWA', 'BFA', 'BFA', 'BDI', 'BDI', 'CPV', 'CPV', 'CMR', 'CMR', 'CAF', 'CAF', 'CAF', 'TCD', 'TCD', 'COM', 'COG', 'COG', 'COG', 'COG', 'CIV', 'CIV', 'COD', 'COD', 'COD', 'COD', 'COD', 'COD', 'COD', 'COD', 'COD', 'COD', 'DJI', 'DJI', 'GNQ', 'GNQ', 'ERI', 'ERI', 'ETH', 'ETH', 'ETH', 'ETH', 'ETH', 'ETH', 'GAB', 'GAB', 'GAB', 'GAB', 'GMB', 'GHA', 'GHA', 'GNB', 'GIN', 'GIN', 'KEN', 'KEN', 'KEN', 'KEN', 'LSO', 'LBR', 'MDG', 'MDG', 'MWI', 'MLI', 'MLI', 'MLI', 'MLI', 'MRT', 'MRT', 'MRT', 'MUS', 'MUS', 'MOZ', 'MOZ', 'NAM', 'NAM', 'NAM', 'NAM', 'NER', 'NER', 'NER', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'RWA', 'RWA', 'STP', 'STP', 'SEN', 'SEN', 'SYC', 'SYC', 'SYC', 'SLE', 'SLE', 'SOM', 'SOM', 'SOM', 'SOM', 'SOM', 'ZAF', 'ZAF', 'ZAF', 'ZAF', 'ZAF', 'SSD', 'SSD', 'SDN', 'SDN', 'SDN', 'SDN', 'SWZ', 'TZA', 'TZA', 'TZA', 'TGO', 'TGO', 'UGA', 'UGA', 'UGA', 'UGA', 'ESH', 'ESH', 'ZMB', 'ZMB', 'ZMB', 'ZWE', 'ZWE']
            })

    # Automatically generated file tracking progress in completing job (skips if file is finished)
    # (necessary because there are many inputs files to one combined output file so if an output file exists it doesn't mean its necessarily completed)
    if os.path.isfile(progress_file):
        progress = pd.read_csv(progress_file)
        skip_list = list(progress[progress['outcome'] == 'completed']['file_name'])
        finished_codes = list(progress[progress['outcome'] == 'completed']['country_code'].unique())
        codes = codes[~codes['file_name'].isin(skip_list)]
        if len(skip_list) > 0 and len(codes['file_name'].unique()) == 0:
            raise ValueError('All countries have finished.')

    # Remove countries from country_chunk argument that are marked as completed in progress_file 
    if country_chunk: 
        codes = codes[codes['country_code'].isin(country_chunk)]
        if codes.shape[0] == 0: 
            raise ValueError('All countries in country_chunk argument already finished.')
    else: 
        raise ValueError('Empty country_chunk arg.')
    
    # Check input files in buildings_dir against codes_file
    input_file_dir_list = list(filter(re.compile("africa_").match, sorted(list(os.listdir(Path(buildings_dir))))))
    input_file_config_list = list(codes['file_name'].unique())
    file_in_config_not_in_dir = [x for x in input_file_config_list if x not in set(input_file_dir_list)]
    if len(file_in_config_not_in_dir) > 0:
        raise ValueError(f'Files in codes_file not in input buildings_dir: {file_in_config_not_in_dir}')

    # Check blocks_dir
    input_blocks_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    input_blocks_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_blocks_list] 
    input_blocks_list = [x for x in country_chunk if x not in set(input_blocks_list)]  
    if input_blocks_list: 
        raise ValueError(f'File for {input_blocks_list} missing from blocks_dir.')

    # Invalid country codes in country_chunk argument
    if country_chunk: 
        invalid_country_arg = [x for x in country_chunk if x not in set(codes['country_code'].unique()) and x not in set(finished_codes)]
        if invalid_country_arg: 
            raise ValueError(f'Invalid country code in country_chunk: {invalid_country_arg}')

    logging.info(f"Countries to process: {codes['country_code'].unique()}")
    logging.info(f"Files to process: {codes['file_name'].unique()}")

    # Iterate through country codes
    for code, files in codes.groupby('country_code'): 
        logging.info(f'Country: {code}')
        # gadm_file = Path(gadm_dir) / f'gadm_{code}.parquet'
        blocks_file = Path(blocks_dir) / f'blocks_{code}.parquet'
        codes_sub = codes[codes['country_code'] == code]

        # Iterate through building files for each country code
        for index, row in codes_sub.iterrows():
            logging.info(f'Ecopia file: {row[0]}')
            print(f'Ecopia file: {row[0]}')
            t0_country = time.time()

            # Extract shapefile
            file_path = Path(buildings_dir) / row[0]
            try:
                z = zipfile.ZipFile(file_path)
            except: 
                logging.info("Invalid Zip file : %s", str(row[0]))
                entry = pd.DataFrame({'file_name': [row[0]], 'country_code': [code], 'outcome':['invalid']})   #, 'polygons':[0], 'points':[0]
                if not os.path.isfile(progress_file):
                    entry.to_csv(progress_file, header='column_names', index=False)
                else: 
                    file_log = pd.read_csv(progress_file)
                    file_log = pd.concat([file_log, entry]) 
                    file_log.to_csv(progress_file, mode='w', index=False)
                continue
            z.extractall(path=buildings_dir) 

            file_name = [y for y in sorted(z.namelist()) if y.endswith('shp')] 
            bldgs_polygons = gpd.read_file(Path(buildings_dir) / file_name[0].replace("//", "/"))
            bldgs_polygons = bldgs_polygons.to_crs(epsg=4326)
    
            folder_name = os.path.dirname(Path(buildings_dir) / file_name[0].replace("//", "/"))
            shutil.rmtree(folder_name) 

            # Read GADM file for country
            #gadms = gpd.read_parquet(gadm_file) 

            # Read blocks file for country
            blocks = gpd.read_parquet(blocks_file)
            blocks = blocks[['block_id','gadm_code','country_code','geometry']].to_crs(3395)

            # Create unique id code
            bldgs_polygons['building_id'] = range(len(bldgs_polygons))
            bldgs_polygons = bldgs_polygons[['building_id','geometry']]
            bldgs_polygons['source'] = row[0]

            # Compute area 
            bldgs_points = transform_buildings(gpd_buildings = bldgs_polygons)

            # Confirm formats
            if 'building_area' not in bldgs_points.columns and bldgs_points['geometry'].geom_type.unique().all() == 'Polygon':
                assert bldgs_points.crs in ['epsg:3395','epsg:4326'], "bldgs_points is not epsg:4326 or epsg:3395."
                if bldgs_points.crs == 'epsg:4326': bldgs_points = bldgs_points.to_crs(epsg=3395)
                bldgs_points['building_area'] = bldgs_points.area
            else: 
                assert 'building_area' in bldgs_points.columns and bldgs_points['geometry'].geom_type.unique().all() == 'Point'

            # Create hash code
            precision_level = 18
            bldgs_points['building_geohash'] = list(map(lambda x: pygeohash.encode(x.x,x.y, precision=precision_level), bldgs_points['geometry'].to_list()))
            bldgs_points['hash_group'] = bldgs_points.groupby(['building_geohash']).ngroup()
            bldgs_points['rank_id'] = bldgs_points.groupby(['building_geohash'])['hash_group'].rank(method="first")
            bldgs_points['building_geohash'] = list(map(lambda x, z: str(int(z)) + '_' + pygeohash.encode(x.x,x.y, precision=precision_level), bldgs_points['geometry'].to_list(), bldgs_points['rank_id'].to_list()))

            bldgs_points = bldgs_points[['building_id','building_geohash','building_area','source','geometry']]
            
            # Join blocks file to building points
            #bldgs_points = gpd.sjoin(left_df = bldgs_points, right_df = gadms, how='left', predicate='intersects')
            bldgs_points = sjoin_buildings(block_data = blocks, bldg_id_col = 'building_geohash', bldg_data = bldgs_points, use_dask_geopandas = True) 

            # Remove duplicates
            if bldgs_points['building_geohash'].duplicated().any():
                warnings.warn(f"Dropping {bldgs_points[bldgs_points['building_geohash'].duplicated()].shape[0]} duplicate rows for column building_geohash. Keeping first occurrence.", UserWarning)
                bldgs_points = bldgs_points.drop_duplicates(subset=['building_geohash'], keep='first')
            
            bldgs_points = bldgs_points[['building_id', 'building_geohash', 'block_id', 'gadm_code', 'country_code', 'building_area', 'source', 'geometry']]

            # Join blocks in points file to polygon file
            bldgs_polygons = bldgs_polygons.merge(bldgs_points[['building_id','building_geohash','gadm_code', 'country_code']], how='inner', left_on='building_id', right_on='building_id')
            bldgs_polygons = bldgs_polygons[['building_id', 'building_geohash', 'gadm_code', 'country_code', 'source', 'geometry']]
            
            # Write polygon building file (and append if country file exists)
            out_poly = Path(output_dir) / 'buildings' / 'polygons' / f"buildings_polygons_{code}.parquet"
            if out_poly.exists():
                bldgs_append = gpd.read_parquet(out_poly)
                if row[0] in bldgs_append['source'].unique(): 
                    raise ValueError(f"{row[0]} is in polygon output file.")
                bldgs_polygons = pd.concat([bldgs_append, bldgs_polygons]) 
                del bldgs_append
            bldgs_polygons.to_parquet(out_poly)
            del bldgs_polygons
            logging.info(f'{str(row[0])} polygons success.')
            print(str(row[0]) + ' poly success.')

            # Write point building file (and append if country file exists)
            out_point = Path(output_dir) / 'buildings' / 'points' / f"buildings_points_{code}.parquet"
            if out_point.exists():
                bldgs_append = gpd.read_parquet(out_point)
                if row[0] in bldgs_append['source'].unique(): 
                    raise ValueError(f"{row[0]} is in point output file.")
                bldgs_points = pd.concat([bldgs_append, bldgs_points]) 
                del bldgs_append
            bldgs_points.to_parquet(out_point)
            del bldgs_points
            logging.info(f'{str(row[0])} points success.')
            print(str(row[0]) + ' point success.')

            # Update progress file (in case job needs to pick up from checkpoint)
            entry = pd.DataFrame({'file_name': [row[0]], 'country_code': [code], 'outcome':['completed']})
            if not os.path.isfile(progress_file):
                entry.to_csv(progress_file, header='column_names', index=False)
            else: 
                file_log = pd.read_csv(progress_file)
                file_log = pd.concat([file_log, entry]) 
                file_log.to_csv(progress_file, mode='w', index=False)
            logging.info(str(row[0]) + ' finished.')
            print(str(row[0]) + ' finished.')
            t1_country = time.time()
            logging.info(f"Finished: {str(row[0])}, {str(round((t1_country-t0_country)/60,3))} minutes")
        print(str(code) + ' finished.')
        logging.info(f'Finished: {str(code)}')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Combine building files.')
    parser.add_argument('--log_file', required=True, type=Path, dest="log_file", help="Path to log file.")
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="Array of country codes in ISO 3166-1 alpha-3 format.")
    parser.add_argument('--codes_file', required=True, type=Path, dest="codes_file", help="Path file linking file names and ISO country codes." )
    parser.add_argument('--progress_file', required=True, type=Path, dest="progress_file", help="Path file tracking progress.")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks file directory.")
    parser.add_argument('--buildings_dir', required=True, type=Path, dest="buildings_dir", help="Path to Ecopia building file directory.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to top level output directory. Creates /buildings/points/ and /buildings/polygons/ subdirectories.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))



