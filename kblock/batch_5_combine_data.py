
import numpy as np
import pandas as pd
import geopandas as gpd
import pygeos
import momepy
gpd.options.use_pygeos = True
np.seterr(divide = 'ignore')

import dask
import dask.dataframe

import psutil
import logging
import re
import os
from pathlib import Path
import argparse
import warnings; warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

def mem_profile() -> str:
    """
    Return memory usage, str
    """
    mem_use = str(round(100 - psutil.virtual_memory().percent,4))+'% of '+str(round(psutil.virtual_memory().total/1e+9,3))+' GB RAM'
    return mem_use

def main(log_file: Path, country_chunk: list, blocks_dir: Path, population_dir: Path, buildings_dir: Path, complexity_dir: Path, streets_dir: Path, crosswalks_dir: Path, output_dir: Path):

    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    logging.info(f'----------------------')
    logging.info(f'----------------------')

    output_dir_africa = str(output_dir) + '/africa_data'
    output_dir_region = str(output_dir) + '/region_data'
    output_dir_country = str(output_dir) + '/country_data'
    dask_dir = str(output_dir_region) + '/dask'

    Path(output_dir_africa).mkdir(parents=True, exist_ok=True)
    Path(output_dir_region).mkdir(parents=True, exist_ok=True)
    Path(output_dir_country).mkdir(parents=True, exist_ok=True)
    Path(dask_dir).mkdir(parents=True, exist_ok=True)

    input_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_file_list]
    logging.info(f'Country list: {country_list}')

    in_chunk_not_in_inputs = [x for x in country_chunk if x not in set(country_list)]
    if len(in_chunk_not_in_inputs) > 0:
        raise ValueError(f'Input data does not exist for {in_chunk_not_in_inputs} in country_chunk argument.')

    if country_chunk:
        country_list = [x for x in country_chunk if x in set(country_list)]
        if not country_list:
            raise ValueError('Missing input data in blocks directory.')
    else:
        raise ValueError('Empty country_chunk arg.')

    logging.info(f'----------------------')

    africa_files_exist = os.path.isfile(Path(output_dir_africa) / f'africa_geodata.parquet')
    if africa_files_exist:
        completed_country_list = pd.read_parquet(Path(output_dir_africa) / f'africa_geodata.parquet', columns=['country_code'])
        completed_country_list = completed_country_list['country_code'].unique()
        logging.info(f'Completed countries: {completed_country_list}')
        completed_country_list = [x for x in country_list if x not in set(completed_country_list)]
        if len(completed_country_list) > 0:
            africa_files_exist = False

    regional_files_exist = os.path.isfile(Path(output_dir_region) / f'aggregate_regional_geodata.parquet')

    if (africa_files_exist is not True):

        logging.info(f'Processing Africa and regional data.')
        logging.info(f'----------------------')

        # Buildings
        if not os.path.isfile(Path(output_dir_africa) / f'africa_buildings_data.parquet'):

            logging.info(f'Processing buildings.')
            all_buildings = pd.DataFrame({
                'block_id': pd.Series(dtype='str'),
                'building_count': pd.Series(dtype= 'int'),
                'building_area': pd.Series(dtype= 'float64'),

                'bldg_area_count_bin_01_0.50__log10_3.2': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_02_0.75__log10_5.6': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_03_1.00__log10_10': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_04_1.25__log10_17.8': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_05_1.50__log10_31.6': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_06_1.75__log10_56.2': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_07_2.00__log10_100': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_08_2.25__log10_177.8': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_09_2.50__log10_316.2': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_10_2.75__log10_562.3': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_11_3.00__log10_1000': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_12_3.25__log10_1778.3': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_13_3.50__log10_3162.3': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_14_3.75__log10_5623.4': pd.Series(dtype= 'int'),
                'bldg_area_count_bin_15_4.00__log10_10000': pd.Series(dtype= 'int'),

                'bldg_area_m2_bin_01_0.50__log10_3.2': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_02_0.75__log10_5.6': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_03_1.00__log10_10': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_04_1.25__log10_17.8': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_05_1.50__log10_31.6': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_06_1.75__log10_56.2': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_07_2.00__log10_100': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_08_2.25__log10_177.8': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_09_2.50__log10_316.2': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_10_2.75__log10_562.3': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_11_3.00__log10_1000': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_12_3.25__log10_1778.3': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_13_3.50__log10_3162.3': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_14_3.75__log10_5623.4': pd.Series(dtype= 'int'),
                'bldg_area_m2_bin_15_4.00__log10_10000': pd.Series(dtype= 'int')
                })

            for country_code in country_list:
                print(country_code)
                logging.info(f'{country_code}')
                buildings = pd.read_parquet(path = Path(buildings_dir) / 'points' / f'buildings_points_{country_code}.parquet', columns=['block_id', 'building_area'])
                #buildings = buildings.drop(columns='geometry')
                buildings = buildings[['block_id', 'building_area']]
                buildings['building_count'] = int(1)
                buildings['building_area_log'] = np.log10(buildings['building_area']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)

                conditions = [(buildings['building_area_log'] < 0.75),
                                (buildings['building_area_log'] >= 0.75) & (buildings['building_area_log'] < 1),
                                (buildings['building_area_log'] >= 1) & (buildings['building_area_log'] < 1.25),
                                (buildings['building_area_log'] >= 1.25) & (buildings['building_area_log'] < 1.5),
                                (buildings['building_area_log'] >= 1.5) & (buildings['building_area_log'] < 1.75),
                                (buildings['building_area_log'] >= 1.75) & (buildings['building_area_log'] < 2),
                                (buildings['building_area_log'] >= 2) & (buildings['building_area_log'] < 2.25),
                                (buildings['building_area_log'] >= 2.25) & (buildings['building_area_log'] < 2.5),
                                (buildings['building_area_log'] >= 2.5) & (buildings['building_area_log'] < 2.75),
                                (buildings['building_area_log'] >= 2.75) & (buildings['building_area_log'] < 3),
                                (buildings['building_area_log'] >= 3) & (buildings['building_area_log'] < 3.25),
                                (buildings['building_area_log'] >= 3.25) & (buildings['building_area_log'] < 3.5),
                                (buildings['building_area_log'] >= 3.5) & (buildings['building_area_log'] < 3.75),
                                (buildings['building_area_log'] >= 3.75) & (buildings['building_area_log'] < 4),
                                (buildings['building_area_log'] >= 4)]
                bin_area_labels = ['01_0.50__log10_3.2', '02_0.75__log10_5.6', '03_1.00__log10_10', '04_1.25__log10_17.8', '05_1.50__log10_31.6', '06_1.75__log10_56.2', '07_2.00__log10_100', '08_2.25__log10_177.8', '09_2.50__log10_316.2', '10_2.75__log10_562.3', '11_3.00__log10_1000', '12_3.25__log10_1778.3', '13_3.50__log10_3162.3', '14_3.75__log10_5623.4', '15_4.00__log10_10000']
                buildings["bin_area"] = np.select(conditions, bin_area_labels, default=None)

                buildings = buildings[['block_id','bin_area','building_count','building_area']]
                # buildings = buildings.groupby(['block_id','bin_area']).sum(['building_count','building_area']).reset_index()
                buildings = buildings.groupby(['block_id','bin_area'])[['building_count','building_area']].agg('sum').reset_index()

                buildings = buildings.join(pd.get_dummies(buildings["bin_area"], prefix='bldg_area_count_bin'))
                buildings = buildings.join(pd.get_dummies(buildings["bin_area"], prefix='bldg_area_m2_bin'))
                logging.info(f"Memory usage {buildings.shape}, {mem_profile()}")
                bin_area_count_list = ['bldg_area_count_bin_' + s for s in bin_area_labels]
                bin_area_m2_list = ['bldg_area_m2_bin_' + s for s in bin_area_labels]

                buildings[buildings.columns[buildings.columns.isin(bin_area_m2_list)]] = buildings[buildings.columns[buildings.columns.isin(bin_area_m2_list)]].multiply(buildings['building_area'], axis="index")
                building_col_list = ['building_area', 'building_count'] + bin_area_count_list + bin_area_m2_list
                building_col_list = list(buildings.columns[buildings.columns.isin(building_col_list)])
                buildings = buildings[['block_id'] + building_col_list]
                buildings = buildings.groupby(['block_id'])[building_col_list].agg('sum').reset_index()
                all_buildings = pd.concat([all_buildings, buildings], ignore_index=True)

            all_buildings.to_parquet(path = Path(output_dir_africa) / f'africa_buildings_data.parquet')
            del all_buildings
            logging.info(f"Memory usage {mem_profile()}")
            logging.info(f'----------------------')

        # Population
        if not os.path.isfile(Path(output_dir_africa) / f'africa_population_data.parquet'):
            logging.info(f'Processing population.')
            all_population = pd.DataFrame({
                'block_id': pd.Series(dtype='str'),
                'landscan_population': pd.Series(dtype= 'float64'),
                'landscan_population_un': pd.Series(dtype= 'float64'),
                'worldpop_population': pd.Series(dtype= 'float64'),
                'worldpop_population_un': pd.Series(dtype= 'float64')
                })

            for country_code in country_list:
                print(country_code)
                logging.info(f'{country_code}')
                population = pd.read_parquet(path = Path(population_dir) / f'population_{country_code}.parquet')
                population = population[['block_id', 'landscan_population', 'landscan_population_un', 'worldpop_population', 'worldpop_population_un']]
                all_population = pd.concat([all_population, population], ignore_index=True)

            all_population.to_parquet(path = Path(output_dir_africa) / f'africa_population_data.parquet')
            del all_population
            logging.info(f"Memory usage {mem_profile()}")
            logging.info(f'----------------------')

        # Complexity
        if not os.path.isfile(Path(output_dir_africa) / f'africa_complexity_data.parquet'):
            logging.info(f'Processing complexity.')
            all_complexity = pd.DataFrame({
                'block_id': pd.Series(dtype='str'),
                'on_network_street_length': pd.Series(dtype='float64'),
                'off_network_street_length': pd.Series(dtype='float64'),
                'nearest_external_street': pd.Series(dtype='float64'),
                'parcel_count': pd.Series(dtype='int64'),
                'parcel_layers': pd.Series(dtype='str'),
                'k_complexity': pd.Series(dtype='int64')
                })

            for country_code in country_list:
                print(country_code)
                logging.info(f'{country_code}')
                complexity = pd.read_parquet(path = Path(complexity_dir) / f'complexity_{country_code}.parquet')
                complexity = complexity.rename(columns={'building_count':'parcel_count', 'building_layers':'parcel_layers'})
                complexity['parcel_layers'] = complexity['parcel_layers'].astype('str')
                complexity = complexity[['block_id', 'on_network_street_length', 'off_network_street_length', 'nearest_external_street', 'parcel_count', 'parcel_layers', 'k_complexity']]
                all_complexity = pd.concat([all_complexity, complexity], ignore_index=True)

            all_complexity.to_parquet(path = Path(output_dir_africa) / f'africa_complexity_data.parquet')
            del all_complexity
            logging.info(f"Memory usage {mem_profile()}")
            logging.info(f'----------------------')

        # Blocks
        if not os.path.isfile(Path(output_dir_africa) / f'africa_blocks_geodata.parquet'):
            logging.info(f'Processing blocks.')
            all_blocks = gpd.GeoDataFrame({
                'block_id': pd.Series(dtype='str'),
                'block_geohash': pd.Series(dtype='str'),
                'gadm_code': pd.Series(dtype='str'),
                'country_code': pd.Series(dtype='str'),
                'block_area': pd.Series(dtype='float64'),
                'block_perimeter': pd.Series(dtype='float64'),
                'geometry': gpd.GeoSeries(dtype='geometry')}).set_crs(epsg=4326)

            for country_code in country_list:
                print(country_code)
                logging.info(f'{country_code}')
                blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
                blocks = blocks[['block_id', 'block_geohash', 'gadm_code', 'country_code', 'block_area', 'block_perimeter', 'geometry']]
                all_blocks = pd.concat([all_blocks, blocks], ignore_index=True)

            logging.info(f"Memory usage {mem_profile()}")
            all_blocks.to_parquet(path = Path(output_dir_africa) / f'africa_blocks_geodata.parquet')
            del all_blocks
            logging.info(f"Memory usage {mem_profile()}")
            logging.info(f'----------------------')

        # Merge data
        logging.info(f'Merge together.')
        all_buildings = pd.read_parquet(path = Path(output_dir_africa) / f'africa_buildings_data.parquet')
        bin_area_col_list = list(all_buildings.columns[all_buildings.columns.str.contains('_bin_')])
        all_population = pd.read_parquet(path = Path(output_dir_africa) / f'africa_population_data.parquet')
        all_complexity = pd.read_parquet(path = Path(output_dir_africa) / f'africa_complexity_data.parquet')
        all_data = gpd.read_parquet(path = Path(output_dir_africa) / f'africa_blocks_geodata.parquet')

        logging.info(f"Memory usage {mem_profile()}")
        all_data = all_data.set_index('block_id').join(all_population.set_index('block_id'), how='left')
        logging.info(f"Memory usage {mem_profile()}")
        all_data = all_data.join(all_buildings.set_index('block_id'), how='left')
        logging.info(f"Memory usage {mem_profile()}")
        all_data = all_data.join(all_complexity.set_index('block_id'), how='left')
        logging.info(f"Memory usage {mem_profile()}")
        all_data = all_data.reset_index()
        # all_data = pd.merge(left = all_blocks, right = all_population, how = 'left', on = 'block_id')
        # all_data = pd.merge(left = all_data, right = all_buildings, how = 'left', on = 'block_id')
        # all_data = pd.merge(left = all_data, right = all_complexity, how = 'left', on = 'block_id')
        del all_population, all_buildings, all_complexity # all_blocks
        logging.info(f"Memory usage {mem_profile()}")
        all_data['on_network_street_length_na'] = all_data['on_network_street_length'].isnull().astype(int)
        all_data['off_network_street_length_na'] = all_data['off_network_street_length'].isnull().astype(int)
        all_data[bin_area_col_list] = all_data[bin_area_col_list].fillna(value=0)
        all_data['k_complexity'] = all_data['k_complexity'].fillna(value=1)
        logging.info(f'Check for missing.')
        for col in all_data.columns:
            if all_data[col].isnull().sum() > 0:
                all_data[col] = all_data[col].fillna(0)
                print(f'{col}: {all_data[col].isnull().sum()}')
                logging.info(f'{col}: {all_data[col].isnull().sum()}')
        logging.info(f'----------------------')

        # Generate metrics
        logging.info(f'Generate additional metrics.')
        logging.info(f"Memory usage {mem_profile()}")
        all_data = all_data.rename(columns={'block_area':'block_area_m2',
                                            'building_area':'building_area_m2',
                                            'block_perimeter': 'block_perimeter_meters',
                                            'on_network_street_length': 'on_network_street_length_meters',
                                            'off_network_street_length': 'off_network_street_length_meters',
                                            'nearest_external_street': 'nearest_external_street_meters'
                                            })

        logging.info(f"Memory usage {mem_profile()}")
        all_data['block_hectares'] = all_data['block_area_m2'] * 0.0001
        all_data['block_area_km2'] = all_data['block_area_m2'] * 1e-6
        all_data['average_building_area_m2'] = (all_data['building_area_m2'] / all_data['building_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['average_parcel_area_m2'] = (all_data['block_area_m2'] / all_data['parcel_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['average_buildings_per_hectare'] = (all_data['building_count'] / all_data['block_hectares']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['building_to_block_area_ratio'] = (all_data['building_area_m2'] / all_data['block_area_m2']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['k_complexity_weighted_landscan_un'] = all_data['k_complexity'] * all_data['landscan_population_un']
        all_data['k_complexity_weighted_worldpop_un'] = all_data['k_complexity'] * all_data['worldpop_population_un']
        all_data['landscan_population_un_log'] = np.log10(all_data['landscan_population_un']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_data['landscan_population_un_density_hectare'] = (all_data['landscan_population_un'] / all_data['block_hectares']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['landscan_population_un_density_hectare_log'] = np.log10(all_data['landscan_population_un_density_hectare']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_data['landscan_population_un_per_building_area_m2'] = (all_data['landscan_population_un'] / all_data['building_area_m2']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['landscan_population_un_per_building'] = (all_data['landscan_population_un'] / all_data['building_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['worldpop_population_un_log'] = np.log10(all_data['worldpop_population_un']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_data['worldpop_population_un_density_hectare'] = (all_data['worldpop_population_un'] / all_data['block_hectares']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['worldpop_population_un_density_hectare_log'] = np.log10(all_data['worldpop_population_un_density_hectare']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_data['worldpop_population_un_per_building_area_m2'] = (all_data['worldpop_population_un'] / all_data['building_area_m2']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['worldpop_population_un_per_building'] = (all_data['worldpop_population_un'] / all_data['building_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_data['parcel_layers'] = all_data['parcel_layers'].astype('str')

        logging.info(f"Memory usage {mem_profile()}")
        mokoko = ['NGA.25.14_1_802', 'NGA.25.14_1_696', 'NGA.25.14_1_697', 'NGA.25.14_1_698', 'NGA.25.14_1_699', 'NGA.25.14_1_700', 'NGA.25.14_1_714', 'NGA.25.14_1_715', 'NGA.25.14_1_716', 'NGA.25.14_1_721', 'NGA.25.14_1_722', 'NGA.25.14_1_723', 'NGA.25.14_1_724', 'NGA.25.14_1_727', 'NGA.25.14_1_728', 'NGA.25.14_1_729', 'NGA.25.14_1_730', 'NGA.25.14_1_731', 'NGA.25.14_1_732', 'NGA.25.14_1_733', 'NGA.25.14_1_734', 'NGA.25.14_1_735', 'NGA.25.14_1_736', 'NGA.25.14_1_737', 'NGA.25.14_1_738', 'NGA.25.14_1_739', 'NGA.25.14_1_740', 'NGA.25.14_1_741', 'NGA.25.14_1_742', 'NGA.25.14_1_743', 'NGA.25.14_1_744', 'NGA.25.14_1_745', 'NGA.25.14_1_746', 'NGA.25.14_1_747', 'NGA.25.14_1_748', 'NGA.25.14_1_749', 'NGA.25.14_1_750', 'NGA.25.14_1_751', 'NGA.25.14_1_752', 'NGA.25.14_1_753', 'NGA.25.14_1_754', 'NGA.25.14_1_755', 'NGA.25.14_1_756', 'NGA.25.14_1_757', 'NGA.25.14_1_758', 'NGA.25.14_1_759', 'NGA.25.14_1_760', 'NGA.25.14_1_761', 'NGA.25.14_1_762', 'NGA.25.14_1_763', 'NGA.25.14_1_764', 'NGA.25.14_1_765', 'NGA.25.14_1_766', 'NGA.25.14_1_767', 'NGA.25.14_1_768', 'NGA.25.14_1_769', 'NGA.25.14_1_770', 'NGA.25.14_1_771', 'NGA.25.14_1_772', 'NGA.25.14_1_773', 'NGA.25.14_1_774', 'NGA.25.14_1_775', 'NGA.25.14_1_776', 'NGA.25.14_1_777', 'NGA.25.14_1_778', 'NGA.25.14_1_779', 'NGA.25.14_1_780', 'NGA.25.14_1_781', 'NGA.25.14_1_782', 'NGA.25.14_1_783', 'NGA.25.14_1_784', 'NGA.25.14_1_785', 'NGA.25.14_1_786', 'NGA.25.14_1_787', 'NGA.25.14_1_788', 'NGA.25.14_1_789', 'NGA.25.14_1_790', 'NGA.25.14_1_791', 'NGA.25.14_1_792', 'NGA.25.14_1_793', 'NGA.25.14_1_794', 'NGA.25.14_1_795', 'NGA.25.14_1_796', 'NGA.25.14_1_797', 'NGA.25.14_1_798', 'NGA.25.14_1_799', 'NGA.25.14_1_800', 'NGA.25.14_1_801', 'NGA.25.14_1_803', 'NGA.25.14_1_804', 'NGA.25.14_1_805', 'NGA.25.14_1_807', 'NGA.25.14_1_808', 'NGA.25.14_1_809', 'NGA.25.14_1_810', 'NGA.25.14_1_811', 'NGA.25.14_1_812', 'NGA.25.14_1_813', 'NGA.25.14_1_814', 'NGA.25.14_1_815', 'NGA.25.14_1_816', 'NGA.25.14_1_817', 'NGA.25.14_1_818', 'NGA.25.14_1_819', 'NGA.25.14_1_820', 'NGA.25.14_1_821', 'NGA.25.14_1_822', 'NGA.25.14_1_823', 'NGA.25.14_1_824', 'NGA.25.14_1_825', 'NGA.25.14_1_826', 'NGA.25.14_1_827', 'NGA.25.14_1_828', 'NGA.25.14_1_829', 'NGA.25.14_1_830', 'NGA.25.14_1_831', 'NGA.25.14_1_832', 'NGA.25.14_1_852', 'NGA.25.14_1_853', 'NGA.25.14_1_854', 'NGA.25.14_1_855', 'NGA.25.14_1_856']
        kibera = ['KEN.30.11_1_1031', 'KEN.30.11_1_1166', 'KEN.30.11_1_1047', 'KEN.30.11_1_1033', 'KEN.30.10_1_751', 'KEN.30.11_1_1172', 'KEN.30.11_1_1167', 'KEN.30.11_1_1030', 'KEN.30.11_1_1065', 'KEN.30.10_1_754', 'KEN.30.10_1_1005', 'KEN.30.10_1_652', 'KEN.30.10_1_651', 'KEN.30.11_1_1214', 'KEN.30.11_1_1103', 'KEN.30.10_1_929', 'KEN.30.11_1_1029', 'KEN.30.10_1_721', 'KEN.30.11_1_1032', 'KEN.30.11_1_1041', 'KEN.30.11_1_1522', 'KEN.30.11_1_1073', 'KEN.30.10_1_979', 'KEN.30.10_1_745', 'KEN.30.10_1_649', 'KEN.30.11_1_1238', 'KEN.30.11_1_1147', 'KEN.30.11_1_1021', 'KEN.30.11_1_1243', 'KEN.30.10_1_916', 'KEN.30.10_1_1200', 'KEN.30.10_1_1210', 'KEN.30.10_1_1155', 'KEN.30.10_1_650', 'KEN.30.11_1_1257', 'KEN.30.10_1_635', 'KEN.30.10_1_984', 'KEN.30.10_1_1153', 'KEN.30.10_1_960', 'KEN.30.11_1_1244', 'KEN.30.10_1_691', 'KEN.30.10_1_739', 'KEN.30.11_1_1141', 'KEN.30.10_1_1202', 'KEN.30.10_1_711']
        # fix spots (will be resolved in future versions by excluding barrier from block definition)
        offnet_list = mokoko + kibera

        conditions = [(all_data['nearest_external_street_meters'] >= 200),
                        (all_data['k_complexity'] == 1),
                        (all_data['k_complexity'] == 2),
                        (all_data['k_complexity'] == 3),
                        (all_data['k_complexity'] == 4),
                        (all_data['k_complexity'] == 5),
                        (all_data['k_complexity'] == 6),
                        (all_data['k_complexity'] == 7),
                        (all_data['k_complexity'] == 8),
                        (all_data['k_complexity'] == 9),
                        (all_data['k_complexity'] >= 10)]
        labels = ['Off-network','1', '2', '3', '4', '5', '6', '7', '8', '9', '10+']
        all_data['k_labels'] = np.select(conditions, labels, default='Off-network')
        all_data.loc[all_data['block_id'].isin(offnet_list), 'k_labels'] = 'Off-network'

        conditions = [(all_data['nearest_external_street_meters'] >= 200),
                        (all_data['k_complexity'] >= 30)]
        labels = ['Off-network','30+']
        all_data['k_labels_detailed'] = np.select(conditions, labels, default= all_data['k_complexity'].astype(int).astype(str))
        all_data.loc[all_data['block_id'].isin(offnet_list), 'k_labels_detailed'] = 'Off-network'

        logging.info(f"Memory usage {mem_profile()}")
        logging.info(f"Merge crosswalks")
        all_xwalk = pd.read_parquet(path = Path(crosswalks_dir) / f'ghsl_crosswalk.parquet')
        all_xwalk = all_xwalk[['block_id', 'country_name', 'continent', 'area_type', 'class_urban_hierarchy','class_urban_periurban_nonurban', 'class_urban_nonurban', 'urban_id','urban_center_name', 'urban_country_code', 'urban_country_name','conurbation_id', 'conurbation_area_name_short', 'conurbation_country_code','conurbation_country_name', 'agglosid', 'agglosname', 'metropole', 'urban_layer_code']]
        #all_data = pd.merge(left = all_data, right = all_xwalk, how = 'left', on = 'block_id')
        all_data = all_data.set_index('block_id').join(all_xwalk.set_index('block_id'), how='left')
        logging.info(f"Memory usage {mem_profile()}")
        all_data = all_data.reset_index()

        all_data_col_list = ['block_id', 'block_geohash', 'block_area_m2', 'block_hectares', 'block_area_km2', 'block_perimeter_meters', 'building_area_m2', 'building_count', 'average_building_area_m2', 'building_to_block_area_ratio', 'parcel_count', 'average_parcel_area_m2', 'parcel_layers', 'k_complexity', 'k_labels', 'k_labels_detailed', 'k_complexity_weighted_landscan_un', 'k_complexity_weighted_worldpop_un', 'landscan_population', 'landscan_population_un', 'landscan_population_un_log', 'landscan_population_un_density_hectare', 'landscan_population_un_density_hectare_log', 'landscan_population_un_per_building_area_m2', 'landscan_population_un_per_building', 'worldpop_population', 'worldpop_population_un', 'worldpop_population_un_log', 'worldpop_population_un_density_hectare', 'worldpop_population_un_density_hectare_log', 'worldpop_population_un_per_building_area_m2', 'worldpop_population_un_per_building', 'on_network_street_length_meters', 'off_network_street_length_meters', 'nearest_external_street_meters', 'on_network_street_length_na', 'off_network_street_length_na', 'gadm_code', 'country_code', 'country_name', 'continent', 'area_type', 'class_urban_hierarchy', 'class_urban_periurban_nonurban', 'class_urban_nonurban', 'urban_id', 'urban_center_name', 'urban_country_code', 'urban_country_name', 'conurbation_id', 'conurbation_area_name_short', 'conurbation_country_code', 'conurbation_country_name', 'agglosid', 'agglosname', 'metropole', 'urban_layer_code'] + bin_area_col_list + ['geometry']
        all_data = all_data[all_data_col_list]

        logging.info(f'Check for missing.')
        for col in all_data.columns:
            if all_data[col].isnull().sum() > 0:
                print(f'{col}: {all_data[col].isnull().sum()}')
                logging.info(f'{col}: {all_data[col].isnull().sum()}')

        dup_count = all_data.duplicated(subset=['block_id']).sum()
        logging.info(f'Duplicates {dup_count}')
        if dup_count > 0:
            logging.info(f'Dropping duplicates.')
            all_data = all_data.drop_duplicates(subset=['block_id'], keep='first')
            assert all_data.duplicated(subset=['block_id']).sum() == 0

        # Write files
        logging.info(f'Writing Africa files.')
        logging.info(f"Memory usage {mem_profile()}")
        all_data.to_parquet(path = Path(output_dir_africa) / f'africa_geodata.parquet')
        
        map_col_list = ['block_id', 'area_type', 'country_name', 'agglosname', 'k_complexity', 'k_labels', 'k_labels_detailed', 'landscan_population_un', 'landscan_population_un_density_hectare', 'worldpop_population_un', 'worldpop_population_un_density_hectare', 'block_hectares', 'building_count', 'average_building_area_m2', 'geometry']
        all_data[map_col_list].to_parquet(path = Path(output_dir_africa) / f'blocktile.parquet')
        all_data.drop(columns='geometry').to_parquet(path = Path(output_dir_africa) / f'africa_data.parquet')
        del all_data
        logging.info(f"Memory usage {mem_profile()}")

        logging.info(f'----------------------')

    if (regional_files_exist is not True):

        all_data = pd.read_parquet(path = Path(output_dir_africa) / f'africa_data.parquet')

        bin_area_col_list = list(all_data.columns[all_data.columns.str.contains('_bin_')])

        # Urban and regional aggregates
        logging.info(f'Processing regional aggregates.')
        all_regions = all_data
        del all_data
        
        mokoko = ['NGA.25.14_1_802', 'NGA.25.14_1_696', 'NGA.25.14_1_697', 'NGA.25.14_1_698', 'NGA.25.14_1_699', 'NGA.25.14_1_700', 'NGA.25.14_1_714', 'NGA.25.14_1_715', 'NGA.25.14_1_716', 'NGA.25.14_1_721', 'NGA.25.14_1_722', 'NGA.25.14_1_723', 'NGA.25.14_1_724', 'NGA.25.14_1_727', 'NGA.25.14_1_728', 'NGA.25.14_1_729', 'NGA.25.14_1_730', 'NGA.25.14_1_731', 'NGA.25.14_1_732', 'NGA.25.14_1_733', 'NGA.25.14_1_734', 'NGA.25.14_1_735', 'NGA.25.14_1_736', 'NGA.25.14_1_737', 'NGA.25.14_1_738', 'NGA.25.14_1_739', 'NGA.25.14_1_740', 'NGA.25.14_1_741', 'NGA.25.14_1_742', 'NGA.25.14_1_743', 'NGA.25.14_1_744', 'NGA.25.14_1_745', 'NGA.25.14_1_746', 'NGA.25.14_1_747', 'NGA.25.14_1_748', 'NGA.25.14_1_749', 'NGA.25.14_1_750', 'NGA.25.14_1_751', 'NGA.25.14_1_752', 'NGA.25.14_1_753', 'NGA.25.14_1_754', 'NGA.25.14_1_755', 'NGA.25.14_1_756', 'NGA.25.14_1_757', 'NGA.25.14_1_758', 'NGA.25.14_1_759', 'NGA.25.14_1_760', 'NGA.25.14_1_761', 'NGA.25.14_1_762', 'NGA.25.14_1_763', 'NGA.25.14_1_764', 'NGA.25.14_1_765', 'NGA.25.14_1_766', 'NGA.25.14_1_767', 'NGA.25.14_1_768', 'NGA.25.14_1_769', 'NGA.25.14_1_770', 'NGA.25.14_1_771', 'NGA.25.14_1_772', 'NGA.25.14_1_773', 'NGA.25.14_1_774', 'NGA.25.14_1_775', 'NGA.25.14_1_776', 'NGA.25.14_1_777', 'NGA.25.14_1_778', 'NGA.25.14_1_779', 'NGA.25.14_1_780', 'NGA.25.14_1_781', 'NGA.25.14_1_782', 'NGA.25.14_1_783', 'NGA.25.14_1_784', 'NGA.25.14_1_785', 'NGA.25.14_1_786', 'NGA.25.14_1_787', 'NGA.25.14_1_788', 'NGA.25.14_1_789', 'NGA.25.14_1_790', 'NGA.25.14_1_791', 'NGA.25.14_1_792', 'NGA.25.14_1_793', 'NGA.25.14_1_794', 'NGA.25.14_1_795', 'NGA.25.14_1_796', 'NGA.25.14_1_797', 'NGA.25.14_1_798', 'NGA.25.14_1_799', 'NGA.25.14_1_800', 'NGA.25.14_1_801', 'NGA.25.14_1_803', 'NGA.25.14_1_804', 'NGA.25.14_1_805', 'NGA.25.14_1_807', 'NGA.25.14_1_808', 'NGA.25.14_1_809', 'NGA.25.14_1_810', 'NGA.25.14_1_811', 'NGA.25.14_1_812', 'NGA.25.14_1_813', 'NGA.25.14_1_814', 'NGA.25.14_1_815', 'NGA.25.14_1_816', 'NGA.25.14_1_817', 'NGA.25.14_1_818', 'NGA.25.14_1_819', 'NGA.25.14_1_820', 'NGA.25.14_1_821', 'NGA.25.14_1_822', 'NGA.25.14_1_823', 'NGA.25.14_1_824', 'NGA.25.14_1_825', 'NGA.25.14_1_826', 'NGA.25.14_1_827', 'NGA.25.14_1_828', 'NGA.25.14_1_829', 'NGA.25.14_1_830', 'NGA.25.14_1_831', 'NGA.25.14_1_832', 'NGA.25.14_1_852', 'NGA.25.14_1_853', 'NGA.25.14_1_854', 'NGA.25.14_1_855', 'NGA.25.14_1_856']
        kibera = ['KEN.30.11_1_1031', 'KEN.30.11_1_1166', 'KEN.30.11_1_1047', 'KEN.30.11_1_1033', 'KEN.30.10_1_751', 'KEN.30.11_1_1172', 'KEN.30.11_1_1167', 'KEN.30.11_1_1030', 'KEN.30.11_1_1065', 'KEN.30.10_1_754', 'KEN.30.10_1_1005', 'KEN.30.10_1_652', 'KEN.30.10_1_651', 'KEN.30.11_1_1214', 'KEN.30.11_1_1103', 'KEN.30.10_1_929', 'KEN.30.11_1_1029', 'KEN.30.10_1_721', 'KEN.30.11_1_1032', 'KEN.30.11_1_1041', 'KEN.30.11_1_1522', 'KEN.30.11_1_1073', 'KEN.30.10_1_979', 'KEN.30.10_1_745', 'KEN.30.10_1_649', 'KEN.30.11_1_1238', 'KEN.30.11_1_1147', 'KEN.30.11_1_1021', 'KEN.30.11_1_1243', 'KEN.30.10_1_916', 'KEN.30.10_1_1200', 'KEN.30.10_1_1210', 'KEN.30.10_1_1155', 'KEN.30.10_1_650', 'KEN.30.11_1_1257', 'KEN.30.10_1_635', 'KEN.30.10_1_984', 'KEN.30.10_1_1153', 'KEN.30.10_1_960', 'KEN.30.11_1_1244', 'KEN.30.10_1_691', 'KEN.30.10_1_739', 'KEN.30.11_1_1141', 'KEN.30.10_1_1202', 'KEN.30.10_1_711']
        offnet_list = mokoko + kibera

        conditions = [(all_regions['nearest_external_street_meters'] >= 200),
                        (all_regions['k_complexity'] == 1), (all_regions['k_complexity'] == 2), (all_regions['k_complexity'] == 3),
                        (all_regions['k_complexity'] == 4), (all_regions['k_complexity'] == 5), (all_regions['k_complexity'] == 6),
                        (all_regions['k_complexity'] == 7), (all_regions['k_complexity'] == 8), (all_regions['k_complexity'] == 9),
                        (all_regions['k_complexity'] >= 10)]
        k_bucket = ['off_network','01', '02', '03', '04', '05', '06', '07', '08', '09', '10_plus']
        all_regions['k_bucket'] = np.select(conditions, k_bucket, default='Off-network')
        all_regions.loc[all_regions['block_id'].isin(offnet_list), 'k_bucket'] = 'Off-network'
        
        all_regions['block_count'] = int(1)
        agg_col_list = ['block_count', 'block_area_m2', 'block_hectares', 'block_area_km2', 'block_perimeter_meters', 'building_area_m2', 'building_count', 'parcel_count', 'k_complexity_weighted_landscan_un', 'k_complexity_weighted_worldpop_un', 'landscan_population', 'landscan_population_un', 'worldpop_population', 'worldpop_population_un'] + bin_area_col_list
        all_regions = all_regions[list(all_regions.columns[all_regions.columns.isin(['urban_layer_code', 'k_bucket'] + agg_col_list)])]

        #all_regions = all_regions.groupby(['urban_layer_code','k_bucket']).sum(agg_col_list).reset_index()
        all_regions = all_regions.groupby(['urban_layer_code','k_bucket'])[agg_col_list].agg('sum').reset_index()

        logging.info(f"Memory usage {mem_profile()}")

        all_regions = all_regions.join(pd.get_dummies(all_regions['k_bucket'], prefix='k_ls'))
        all_regions = all_regions.join(pd.get_dummies(all_regions['k_bucket'], prefix='k_wp'))
        k_ls_list = ['k_ls_01','k_ls_02','k_ls_03','k_ls_04','k_ls_05','k_ls_06','k_ls_07','k_ls_08','k_ls_09','k_ls_10_plus','k_ls_off_network']
        k_wp_list = ['k_wp_01','k_wp_02','k_wp_03','k_wp_04','k_wp_05','k_wp_06','k_wp_07','k_wp_08','k_wp_09','k_wp_10_plus','k_wp_off_network']

        k_ls_list = list(all_regions.columns[all_regions.columns.isin(k_ls_list)])
        k_wp_list = list(all_regions.columns[all_regions.columns.isin(k_wp_list)])

        all_regions[k_ls_list] = all_regions[k_ls_list].multiply(all_regions['landscan_population_un'], axis="index")
        all_regions[k_wp_list] = all_regions[k_wp_list].multiply(all_regions['worldpop_population_un'], axis="index")

        agg_col_list = agg_col_list + k_ls_list + k_wp_list

        #all_regions = all_regions.groupby(['urban_layer_code']).sum(agg_col_list).reset_index()
        all_regions = all_regions.groupby(['urban_layer_code'])[agg_col_list].agg('sum').reset_index()

        all_regions['average_building_area_m2'] = (all_regions['building_area_m2'] / all_regions['building_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['average_parcel_area_m2'] = (all_regions['block_area_m2'] / all_regions['parcel_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['average_buildings_per_hectare'] = (all_regions['building_count'] / all_regions['block_hectares']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['average_buildings_per_block'] = (all_regions['building_count'] / all_regions['block_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['buildings_to_block_area_ratio'] = (all_regions['building_area_m2'] / all_regions['block_area_m2']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['k_complexity_weighted_landscan_un'] = (all_regions['k_complexity_weighted_landscan_un'] / all_regions['landscan_population_un']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['k_complexity_weighted_worldpop_un'] = (all_regions['k_complexity_weighted_worldpop_un'] / all_regions['worldpop_population_un']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['landscan_population_un_log'] = np.log10(all_regions['landscan_population_un']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_regions['landscan_population_un_density_hectare'] = (all_regions['landscan_population_un'] / all_regions['block_hectares']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['landscan_population_un_density_hectare_log'] = np.log10(all_regions['landscan_population_un_density_hectare']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_regions['landscan_population_un_per_building_area_m2'] = (all_regions['landscan_population_un'] / all_regions['building_area_m2']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['landscan_population_un_per_building'] = (all_regions['landscan_population_un'] / all_regions['building_count']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['worldpop_population_un_log'] = np.log10(all_regions['worldpop_population_un']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_regions['worldpop_population_un_density_hectare'] = (all_regions['worldpop_population_un'] / all_regions['block_hectares']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['worldpop_population_un_density_hectare_log'] = np.log10(all_regions['worldpop_population_un_density_hectare']).replace([np.inf, -np.inf, np.nan], 0).clip(lower=0)
        all_regions['worldpop_population_un_per_building_area_m2'] = (all_regions['worldpop_population_un'] / all_regions['building_area_m2']).replace([np.inf, -np.inf, np.nan], 0)
        all_regions['worldpop_population_un_per_building'] = (all_regions['worldpop_population_un'] / all_regions['building_count']).replace([np.inf, -np.inf, np.nan], 0)

        conditions = [(all_regions['k_complexity_weighted_landscan_un'] >= 0) & (all_regions['k_complexity_weighted_landscan_un'] < 2),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 2) & (all_regions['k_complexity_weighted_landscan_un'] < 3),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 3) & (all_regions['k_complexity_weighted_landscan_un'] < 4),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 4) & (all_regions['k_complexity_weighted_landscan_un'] < 5),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 5) & (all_regions['k_complexity_weighted_landscan_un'] < 6),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 6) & (all_regions['k_complexity_weighted_landscan_un'] < 7),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 7) & (all_regions['k_complexity_weighted_landscan_un'] < 8),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 8) & (all_regions['k_complexity_weighted_landscan_un'] < 9),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 9) & (all_regions['k_complexity_weighted_landscan_un'] < 10),
                        (all_regions['k_complexity_weighted_landscan_un'] >= 10)]
        labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10+']
        all_regions['k_ls_labels'] = np.select(conditions, labels, default='Off-network')

        conditions = [(all_regions['k_complexity_weighted_worldpop_un'] >= 0) & (all_regions['k_complexity_weighted_worldpop_un'] < 2),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 2) & (all_regions['k_complexity_weighted_worldpop_un'] < 3),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 3) & (all_regions['k_complexity_weighted_worldpop_un'] < 4),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 4) & (all_regions['k_complexity_weighted_worldpop_un'] < 5),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 5) & (all_regions['k_complexity_weighted_worldpop_un'] < 6),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 6) & (all_regions['k_complexity_weighted_worldpop_un'] < 7),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 7) & (all_regions['k_complexity_weighted_worldpop_un'] < 8),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 8) & (all_regions['k_complexity_weighted_worldpop_un'] < 9),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 9) & (all_regions['k_complexity_weighted_worldpop_un'] < 10),
                        (all_regions['k_complexity_weighted_worldpop_un'] >= 10)]
        labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10+']
        all_regions['k_wp_labels'] = np.select(conditions, labels, default='Off-network')

        all_boundaries = gpd.read_parquet(Path(crosswalks_dir) / f'urban_boundaries.parquet')
        all_regions = pd.merge(left = all_boundaries, right = all_regions, how = 'right', on = 'urban_layer_code')

        dup_count = all_regions.duplicated(subset=['urban_layer_code']).sum()

        logging.info(f'Duplicates {dup_count}')
        if dup_count > 0:
            logging.info(f'Dropping duplicates.')
            all_regions = all_regions.drop_duplicates(subset=['urban_layer_code'], keep='first')
            assert all_regions.duplicated(subset=['urban_layer_code']).sum() == 0

        for col in all_regions.columns:
            if all_regions[col].isnull().sum() > 0:
                print(f'{col}: {all_regions[col].isnull().sum()}')

        # Write files
        logging.info(f'Writing regional aggregates.')
        all_regions.to_parquet(path = Path(output_dir_region) / 'aggregate_regional_geodata.parquet')
        all_regions.drop(columns='geometry').to_parquet(path = Path(output_dir_region) /  'aggregate_regional_data.parquet')
        all_regions.to_file(filename = Path(output_dir_region) / 'aggregate_regional_geodata.gpkg', driver='GPKG')
        
        reg_col_list = ['urban_layer_code', 'country_name', 'area_type', 'block_count', 'block_area_km2', 'building_count', 'average_building_area_m2', 'k_ls_labels', 'k_wp_labels', 'k_complexity_weighted_landscan_un', 'k_complexity_weighted_worldpop_un', 'landscan_population_un', 'worldpop_population_un', 'landscan_population_un_density_hectare', 'worldpop_population_un_density_hectare', 'geometry']
        all_regions[reg_col_list].to_parquet(path = Path(output_dir_region) / f'regiontile.parquet')
        # all_regions.drop(columns='geometry').to_csv(path_or_buf = Path(output_dir_region) / 'aggregate_regional_data.csv', index=False)
        logging.info(f"Memory usage {mem_profile()}")

    logging.info(f'----------------------')

    logging.info(f'Writing country files.')
    africa_files_exist = os.path.isfile(Path(output_dir_africa) / f'africa_geodata.parquet')

    if africa_files_exist:

        completed_country_list = os.listdir(Path(output_dir_country))
        if len(completed_country_list) > 0:
            remaining_country_list = [x for x in country_list if x not in set(completed_country_list)]
        else:
            remaining_country_list = country_list

        if len(remaining_country_list) > 0:
            for country_code in remaining_country_list:
                logging.info(f'{country_code}')
                print(country_code)
                Path(str(output_dir_country) + f'/{country_code}').mkdir(parents=True, exist_ok=True)
                country_data = gpd.read_parquet(Path(output_dir_africa) / f'africa_geodata.parquet', filters = [('country_code', 'in', [country_code])])
                country_data_sub = country_data.loc[country_data['country_code'] == country_code]
                country_data_sub = country_data_sub.drop(columns=['k_labels', 'k_labels_detailed'])
                country_data_sub.to_parquet(Path(output_dir_country) / f'{country_code}' / f'{country_code}_geodata.parquet')
                country_data_sub.to_file(Path(output_dir_country) / f'{country_code}' / f'{country_code}_geodata.gpkg', driver="GPKG")
                # country_data.drop(columns='geometry').to_csv(path_or_buf = Path(output_dir_country) / f'{country_code}' / f'{country_code}_data.csv', index=False)

    logging.info(f'----------------------')

    # Street network metrics
    logging.info(f'Processing street networks.')
    all_boundaries = gpd.read_parquet(Path(crosswalks_dir) / f'urban_boundaries.parquet')
    all_boundaries = all_boundaries[['urban_layer_code','country_code','geometry']].reset_index(drop= True)

    streets_metrics = pd.DataFrame({
        'urban_layer_code': pd.Series(dtype='str'),
        'country_code': pd.Series(dtype='str'),
        'osm_highway_tag': pd.Series(dtype='str'),
        'street_length_meters': pd.Series(dtype='float64'),
        'mean_linearity': pd.Series(dtype='float64')
        })

    if not os.path.isfile(Path(output_dir_region) / 'streets_metrics.parquet'):
        street_list = country_list
    else:
        if len(os.listdir(Path(dask_dir))) > 0:
            streets_metrics = dask.dataframe.read_parquet(path = Path(dask_dir) / f'streets.parquet').compute()
            completed_streets_list = list(streets_metrics['country_code'].unique())
            street_list = [x for x in country_list if x not in set(completed_streets_list)]
        else:
            street_list = country_list

    if len(street_list) > 0:
        for country_code in street_list:
            print(country_code)
            logging.info(f'{country_code}')
            streets = gpd.read_parquet(path = Path(streets_dir) / f'streets_{country_code}.parquet', columns=['id','highway','geometry'])
            streets = streets[['id','highway','geometry']].reset_index(drop=True)
            streets = streets.rename(columns={'highway':'osm_highway_tag'})
            streets = streets.explode(ignore_index=True)
            streets.sindex
            logging.info(f"Memory usage {mem_profile()}")
            subset_boundaries = all_boundaries[all_boundaries['country_code'] == country_code].copy()
            subset_boundaries['area'] = subset_boundaries.to_crs(3395).area
            subset_boundaries = subset_boundaries.sort_values(by='area', ascending=True).reset_index(drop=True)
            subset_boundaries.sindex
            boundary_list = subset_boundaries['urban_layer_code'].unique()

            for boundary in boundary_list:
                print(boundary)
                boundary_polys = subset_boundaries[subset_boundaries['urban_layer_code'] == boundary][['urban_layer_code','country_code','geometry']]
                if boundary_polys['geometry'].to_crs(3395).area.sum() > 1000000000:
                    boundary_polys = pygeos.multipolygons(pygeos.get_parts(pygeos.from_shapely(boundary_polys['geometry'])))
                    boundary_polys = pygeos.make_valid(boundary_polys)
                    for tag in streets['osm_highway_tag'].unique():
                        streets_lines = pygeos.multilinestrings(pygeos.get_parts(pygeos.from_shapely(streets[streets['osm_highway_tag'] == tag]['geometry'])))
                        streets_lines = pygeos.make_valid(streets_lines)
                        streets_lines = pygeos.intersection_all([streets_lines, boundary_polys])
                        if not pygeos.is_empty(streets_lines).all():
                            streets_lines = gpd.GeoDataFrame.from_dict({'geometry': gpd.GeoSeries(pygeos.to_shapely(streets_lines))}).set_crs(4326)
                            streets_lines = streets_lines.explode(ignore_index=True)
                            streets_lines['urban_layer_code'] = boundary
                            streets_lines['country_code'] = country_code
                            streets_lines['osm_highway_tag'] = tag
                            streets_lines['street_length_meters'] = (streets_lines.to_crs(3395).length).replace([np.inf, -np.inf, np.nan], 0)
                            streets_lines['mean_linearity'] = (momepy.Linearity(streets_lines['geometry'].to_crs(3395)).series * streets_lines['street_length_meters']).replace([np.inf, -np.inf, np.nan], 0)
                            streets_lines = streets_lines[['urban_layer_code','country_code','osm_highway_tag','street_length_meters','mean_linearity']]
                            streets_lines[['street_length_meters','mean_linearity']] = streets_lines[['street_length_meters','mean_linearity']].fillna(0)
                            # streets_lines = streets_lines.groupby(['urban_layer_code','country_code','osm_highway_tag']).sum(['street_length_meters','mean_linearity']).reset_index()
                            streets_lines = streets_lines.groupby(['urban_layer_code','country_code','osm_highway_tag'])[['street_length_meters','mean_linearity']].agg('sum').reset_index()
                            try:
                                streets_lines['mean_linearity'] = (streets_lines['mean_linearity'] / streets_lines['street_length_meters']).replace([np.inf, -np.inf, np.nan], 0)
                            except:
                                streets_lines['mean_linearity'] = 0
                            streets_metrics = pd.concat([streets_metrics, streets_lines])
                else:
                    streets_lines = gpd.overlay(df1 = streets, df2 = boundary_polys[['urban_layer_code','country_code','geometry']], how='intersection', keep_geom_type = True, make_valid = True)
                    streets_lines['street_length_meters'] = (streets_lines.to_crs(3395).length).replace([np.inf, -np.inf, np.nan], 0)
                    streets_lines['mean_linearity'] = (momepy.Linearity(streets_lines['geometry'].to_crs(3395)).series * streets_lines['street_length_meters']).replace([np.inf, -np.inf, np.nan], 0)
                    streets_lines = streets_lines[['urban_layer_code','country_code','osm_highway_tag','street_length_meters','mean_linearity']]
                    streets_lines[['street_length_meters','mean_linearity']] = streets_lines[['street_length_meters','mean_linearity']].fillna(0)
                    # streets_lines = streets_lines.groupby(['urban_layer_code','country_code','osm_highway_tag']).sum(['street_length_meters','mean_linearity']).reset_index()
                    streets_lines = streets_lines.groupby(['urban_layer_code','country_code','osm_highway_tag'])[['street_length_meters','mean_linearity']].agg('sum').reset_index()
                    try:
                        streets_lines['mean_linearity'] = (streets_lines['mean_linearity'] / streets_lines['street_length_meters']).replace([np.inf, -np.inf, np.nan], 0)
                    except:
                        streets_lines['mean_linearity'] = 0
                    streets_metrics = pd.concat([streets_metrics, streets_lines])
                # boundary_residual = subset_boundaries[subset_boundaries['urban_layer_code'].isin([streets_metrics['urban_layer_code'].unique()])][['geometry']]
                # streets = gpd.overlay(df1 = streets, df2 = boundary_residual, how='difference', keep_geom_type=True, make_valid=True)
            streets_metrics_dask = dask.dataframe.from_pandas(data = streets_metrics, npartitions = 1)
            dask.dataframe.to_parquet(df = streets_metrics_dask, path = Path(dask_dir) / f'streets.parquet', append=True, ignore_divisions=True)
            # Write to dask dataframe
            del streets, subset_boundaries, boundary_list, boundary_polys, streets_lines
        principal_roads = ['motorway','trunk','primary','secondary','tertiary','unclassified','residential']
        link_and_special_roads = ['motorway_link','trunk_link','road','primary_link','secondary_link','tertiary_link','busway','bus_guideway','track','living_street','escape','service','pedestrian']
        other_roads = ['raceway','footway','bridleway','steps','corridor','path','via_ferrata']
        conditions = [(streets_metrics['osm_highway_tag'].isin(principal_roads)),
                        (streets_metrics['osm_highway_tag'].isin(link_and_special_roads)),
                        (streets_metrics['osm_highway_tag'].isin(other_roads))]
        labels = ['Principal roads', "Link and special roads", "Other roads"]
        streets_metrics["highway_classification"] = np.select(conditions, labels, default='Other roads')
        vehicular_roads = ['motorway', 'trunk', 'primary', 'secondary', 'tertiary',  'unclassified', 'residential', 'motorway_link',  'trunk_link',  'road', 'primary_link',  'secondary_link',  'living_street']
        streets_metrics["vehicular_classification"] = np.select([(streets_metrics['osm_highway_tag'].isin(vehicular_roads))], ['Vehicular roads'], default='Non-vehicular roads')
        streets_metrics = streets_metrics[['urban_layer_code', 'country_code', 'osm_highway_tag', 'highway_classification', 'vehicular_classification', 'street_length_meters','mean_linearity']]
        logging.info(f'Writing files.')
        streets_metrics.to_parquet(Path(output_dir_region) / 'streets_metrics.parquet')
        del streets_metrics

    logging.info(f'----------------------')
    logging.info(f'----------------------')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Combine files.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file.")
    parser.add_argument('--country_chunk', required=True, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format.")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory.")
    parser.add_argument('--population_dir', required=True, type=Path, dest="population_dir", help="Path to population directory.")
    parser.add_argument('--buildings_dir', required=True, type=Path, dest="buildings_dir", help="Path to buildings directory.")
    parser.add_argument('--complexity_dir', required=True, type=Path, dest="complexity_dir", help="Path to complexity directory.")
    parser.add_argument('--streets_dir', required=True, type=Path, dest="streets_dir", help="Path to streets directory.")
    parser.add_argument('--crosswalks_dir', required=True, type=Path, dest="crosswalks_dir", help="Path to crosswalks directory.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to output directory.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))
