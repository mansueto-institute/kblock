
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

import re
import os
import shutil
from typing import Union
from pathlib import Path
import argparse
import dask
import dask.dataframe 
import dask_geopandas
import warnings; warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

# ghsl_file = '/Users/nm/Downloads/update/inputs/ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg'
# blocks_dir = '/Users/nm/Downloads/update/outputs/blocks'
# output_dir = '/Users/nm/Downloads/update/outputs'

def main(ghsl_file: Path, blocks_dir: Path, output_dir: Path):

    # Create directories
    dask_dir =  str(output_dir) + '/regions' +  '/dask'
    partition_dir =  str(output_dir) + '/partitions'
    if os.path.isdir(dask_dir): shutil.rmtree(dask_dir)
    if os.path.isdir(partition_dir): shutil.rmtree(partition_dir)
    Path(dask_dir).mkdir(parents=True, exist_ok=True)
    Path(partition_dir).mkdir(parents=True, exist_ok=True)

    # Counrty labels
    country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}
    country_names = pd.DataFrame.from_dict(country_dict).T.reset_index().rename(columns= {'index':'country_code', 0:'country',1:'continent',2:'country_name'}).drop(['country'], axis=1)

    # Prepare GHSL data
    ghsl_data = gpd.read_file(Path(ghsl_file))
    ghsl_data = ghsl_data[ghsl_data['GRGN_L2'].isin(['Western Africa', 'Northern Africa', 'Middle Africa', 'Southern Africa', 'Eastern Africa'])][['ID_HDC_G0','UC_NM_MN','geometry']].reset_index()
    col_recode = {'ID_HDC_G0':'urban_id', 'UC_NM_MN':'urban_center_name'}
    ghsl_data = ghsl_data.rename(columns=col_recode)
    ghsl_data = ghsl_data[['urban_id','urban_center_name','geometry']] 
    ghsl_data = ghsl_data[ghsl_data['urban_center_name'] != 'N/A']
    ghsl_data = ghsl_data.explode(index_parts=False)
    ghsl_data['urban_id'] = ghsl_data['urban_id'].astype('int64').astype('string')
    ghsl_data['urban_id'] = 'ghsl_' + ghsl_data['urban_id']
    ghsl_data['urban_id'][ghsl_data['urban_id'] == 'ghsl_2544'] = 'ghsl_2565' # fix Abuja
    ghsl_data = ghsl_data[['urban_id','urban_center_name','geometry']].dissolve(by=['urban_id','urban_center_name'], as_index=False)
    ghsl_data = ghsl_data.to_crs(3395)
    ghsl_data['urban_area'] = ghsl_data['geometry'].area
    assert ghsl_data[ghsl_data['urban_id'].duplicated()].shape[0] == 0

    # Define urban centers
    urban_centers = ghsl_data[['urban_id','geometry']].to_crs(4326)
    urban_centers['area_type'] = 'Urban'

    # Create conurbations 10km buffer
    conurbation_buffer = pygeos.from_shapely(ghsl_data['geometry'].to_crs(3395).buffer(10000))
    conurbation_buffer = pygeos.get_parts(pygeos.union_all(conurbation_buffer))
    conurbation_buffer = gpd.GeoDataFrame.from_dict({'geometry': pygeos.to_shapely(conurbation_buffer)}).set_crs(epsg=3395)
    conurbation_buffer = conurbation_buffer.assign(conurbation_id = [str(x) for x in list(conurbation_buffer.index)])
    conurbation_buffer = conurbation_buffer[['conurbation_id','geometry']].to_crs(4326)

    # Define peri-urban areas (difference between conurban and urban)
    periurban_buffer = gpd.overlay(df1 = conurbation_buffer, df2 = urban_centers, how = 'difference')
    periurban_buffer['urban_id'] = 'periurban_' + periurban_buffer['conurbation_id']
    periurban_buffer['area_type'] = 'Peri-urban'
    periurban_buffer = periurban_buffer[['urban_id','area_type','geometry']]
    urban_periurban_buffer = pd.concat([urban_centers, periurban_buffer])

    # Create list of countries based on block directory
    input_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_file_list]
    full_xwalk = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'urban_id': pd.Series(dtype='float'), 'conurbation_id': pd.Series(dtype='float')})

    region_delineations = gpd.GeoDataFrame({'urban_regional_layer_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326)
    urban_delineations = gpd.GeoDataFrame({'urban_layer_code': pd.Series(dtype='str'), 'geometry': pd.Series(dtype='geometry')}).set_crs(epsg=4326)

    # Loop through countries
    for country_code in country_list:
        print(country_code)
    
        # Overlay urban, periurban layers on blocks
        blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
        urban_periurban_blocks = gpd.overlay(df1 = blocks, df2 = urban_periurban_buffer, how='intersection', keep_geom_type=None, make_valid=True)
        non_urban_blocks = gpd.overlay(df1 = blocks, df2 = urban_periurban_buffer, how='difference', keep_geom_type=None, make_valid=True)
        non_urban_blocks['urban_id'] = 'nonurban_' + non_urban_blocks['country_code']
        non_urban_blocks['area_type'] = 'Non-urban'
        blocks_class = pd.concat([urban_periurban_blocks, non_urban_blocks])
        blocks_class['area'] = blocks_class['geometry'].to_crs(3395).area
        blocks_class['rank'] = blocks_class.groupby('block_id')['area'].rank(method='first', ascending=False)

        # Ensure urban centers are contiguous
        urban_centers_refined = blocks_class[blocks_class['rank'] == 1]
        urban_centers_refined = urban_centers_refined[['block_id','urban_id','area_type']]
        urban_centers_refined = urban_centers_refined[urban_centers_refined['area_type'] == 'Urban']
        urban_centers_refined = blocks.merge(urban_centers_refined, how = 'inner', left_on=['block_id'], right_on = ['block_id'])
        urban_centers_refined = urban_centers_refined[['urban_id','geometry']].dissolve(by=['urban_id'], as_index=False)
        urban_centers_refined['geometry'] = urban_centers_refined['geometry'].make_valid()
        urban_centers_refined = urban_centers_refined.explode(index_parts=False)
        urban_centers_refined = urban_centers_refined[urban_centers_refined['geometry'].geom_type == 'Polygon']
        urban_centers_refined['area'] = urban_centers_refined['geometry'].to_crs(3395).area
        urban_centers_refined['rank'] = urban_centers_refined.groupby('urban_id')['area'].rank(method='first', ascending=False)
        urban_centers_refined = urban_centers_refined[urban_centers_refined['rank'] == 1][['urban_id','geometry']]
        urban_centers_refined['area_type'] = 'Urban'

        # Fit periurban areas to contiguous urban areas
        periurban_buffer_refined = gpd.overlay(df1 = conurbation_buffer, df2 = urban_centers_refined, how = 'difference')
        periurban_buffer_refined['urban_id'] = 'periurban_' + periurban_buffer_refined['conurbation_id']
        periurban_buffer_refined['area_type'] = 'Peri-urban'
        periurban_buffer_refined = periurban_buffer_refined[['urban_id','area_type','geometry']]

        # Combine urban and periurban, join in conurbations
        urban_periurban_refined = pd.concat([urban_centers_refined, periurban_buffer_refined])
        urban_periurban_refined = gpd.sjoin(left_df = urban_periurban_refined, right_df = conurbation_buffer, how = 'left')
        urban_periurban_refined['conurbation_id'] = 'conurban_' + urban_periurban_refined['conurbation_id']
        
        assert urban_periurban_refined[urban_periurban_refined['urban_id'].duplicated()].shape[0] == 0
        assert urban_periurban_refined[urban_periurban_refined['conurbation_id'].isnull()].shape[0] == 0

        # Ensure blocks are one to one with boundaries
        urban_periurban_blocks = gpd.overlay(df1 = blocks, df2 = urban_periurban_refined, how='intersection', keep_geom_type=True, make_valid=True)
        urban_periurban_blocks['area'] = urban_periurban_blocks['geometry'].to_crs(3395).area
        urban_periurban_blocks['rank'] = urban_periurban_blocks.groupby('block_id')['area'].rank(method='first', ascending=False)
        urban_periurban_blocks = urban_periurban_blocks[urban_periurban_blocks['rank'] == 1][['block_id','urban_id','area_type','conurbation_id']]

        # Join boundary labels to blocks
        blocks_xwalk = blocks.merge(urban_periurban_blocks, how = 'left', left_on=['block_id'], right_on = ['block_id'])
        blocks_xwalk['urban_id'] = blocks_xwalk['urban_id'].fillna('nonurban_' + blocks_xwalk['country_code'])
        blocks_xwalk['conurbation_id'] = blocks_xwalk['conurbation_id'].fillna('nonurban_' + blocks_xwalk['country_code'])
        blocks_xwalk['area_type'] = blocks_xwalk['area_type'].fillna('Non-urban')
        assert blocks_xwalk[blocks_xwalk['block_id'].duplicated()].shape[0] == 0
        assert blocks_xwalk[blocks_xwalk['block_id'].isnull()].shape[0] == 0
            
        # Urban and regional delineations
        blocks_xwalk['urban_regional_layer_code'] = blocks_xwalk['country_code'] + '_' + blocks_xwalk['conurbation_id'] + '_' + blocks_xwalk['urban_id']
        blocks_xwalk['urban_regional_layer_code'] = blocks_xwalk['urban_regional_layer_code'].fillna(blocks_xwalk['gadm_code'] + '_nonurban')
        blocks_xwalk['urban_layer_code'] = blocks_xwalk['country_code'] + '_' + blocks_xwalk['conurbation_id'] + '_' + blocks_xwalk['urban_id']
        blocks_xwalk['urban_layer_code'] = blocks_xwalk['urban_layer_code'].fillna(blocks_xwalk['country_code'] + '_nonurban')
    
        region_boundary = dask_geopandas.from_geopandas(blocks_xwalk, npartitions = 50)
        region_boundary = region_boundary.dissolve(by = 'urban_regional_layer_code').compute()
        region_boundary.reset_index(inplace = True)
        region_boundary = region_boundary[['urban_regional_layer_code','geometry']]
        region_delineations = pd.concat([region_delineations, region_boundary], ignore_index=True)

        urban_boundary = dask_geopandas.from_geopandas(blocks_xwalk, npartitions = 50)
        urban_boundary = urban_boundary.dissolve(by = 'urban_layer_code').compute()
        urban_boundary.reset_index(inplace = True)
        urban_boundary = urban_boundary[['urban_layer_code','geometry']]
        urban_delineations = pd.concat([urban_delineations, urban_boundary], ignore_index=True)

        blocks_xwalk = blocks_xwalk[['block_id','block_geohash','gadm_code','country_code','urban_id','conurbation_id','area_type','urban_regional_layer_code','urban_layer_code']]
        blocks_xwalk = dask.dataframe.from_pandas(data = blocks_xwalk, npartitions = 1) 
        dask.dataframe.to_parquet(df = blocks_xwalk, path = Path(dask_dir) / f'{country_code}.parquet', engine='pyarrow', compression='snappy', append=True, ignore_divisions=True)

    # Read parquet data
    urban_delineations.to_parquet(path = Path(output_dir) / f'urban_boundaries.parquet')
    region_delineations.to_parquet(path = Path(output_dir) / f'urban_regional_boundaries.parquet')

    print('Reading countries')
    full_xwalk = dask.dataframe.read_parquet(path = Path(dask_dir)).compute()

    # Join GHSL data to crosswalk
    full_xwalk = full_xwalk.merge(ghsl_data[['urban_id','urban_center_name','urban_area']], how = 'left', left_on=['urban_id'], right_on = ['urban_id'])

    # Coerce all conurbations to have at least one urban_id (applies to blocks on fringes of )
    full_xwalk['urban_conurbation_count'] = full_xwalk.groupby(['conurbation_id'])['urban_id'].transform('count')
    full_xwalk.loc[full_xwalk['urban_conurbation_count'] == 0, 'conurbation_id'] = np.nan
    
    # Create conurbation names
    conurbation_dict = full_xwalk[(full_xwalk['conurbation_id'].notnull()) & (full_xwalk['urban_center_name'].notnull())][['conurbation_id','urban_center_name','urban_area']].drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['conurbation_id'])['urban_center_name'].apply(list).to_dict()
    conurbation_labels = pd.DataFrame({'conurbation_id': conurbation_dict.keys(), "conurbation_area_name": conurbation_dict.values()})
    conurbation_labels['conurbation_area_name'] = ['-'.join(map(str, l)) for l in conurbation_labels ['conurbation_area_name']]
    
    # Join in country and conurbation names to crosswalk
    full_xwalk = full_xwalk.merge(conurbation_labels, how = 'left', left_on=['conurbation_id'], right_on = ['conurbation_id'])
    full_xwalk = full_xwalk.merge(country_names, how = 'left', left_on=['country_code'], right_on = ['country_code'])
    
    # Create country conurbation names
    conurbation_country_dict = full_xwalk.groupby(['conurbation_id','country_name']).agg({'urban_area': 'sum'}).reset_index().drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['conurbation_id'])['country_name'].apply(list).to_dict()
    conurbation_country_name_labels = pd.DataFrame({'conurbation_id': conurbation_country_dict.keys(), "conurbation_country_name": conurbation_country_dict.values()})
    conurbation_country_name_labels["conurbation_country_name"] = ['–'.join(map(str, l)) for l in conurbation_country_name_labels["conurbation_country_name"]]
    
    conurbation_country_dict = full_xwalk.groupby(['conurbation_id','country_code']).agg({'urban_area': 'sum'}).reset_index().drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['conurbation_id'])['country_code'].apply(list).to_dict()
    conurbation_country_code_labels = pd.DataFrame({'conurbation_id': conurbation_country_dict.keys(), "conurbation_country_code": conurbation_country_dict.values()})
    conurbation_country_code_labels["conurbation_country_code"] = ['–'.join(map(str, l)) for l in conurbation_country_code_labels["conurbation_country_code"]]
    
    # Create country urban names
    urban_country_dict = full_xwalk.groupby(['urban_id','country_name']).agg({'urban_area': 'sum'}).reset_index().drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['urban_id'])['country_name'].apply(list).to_dict()
    urban_country_name_labels = pd.DataFrame({'urban_id': urban_country_dict.keys(), 'urban_country_name': urban_country_dict.values()})
    urban_country_name_labels['urban_country_name'] = ['–'.join(map(str, l)) for l in urban_country_name_labels['urban_country_name']]

    urban_country_dict = full_xwalk.groupby(['urban_id','country_code']).agg({'urban_area': 'sum'}).reset_index().drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['urban_id'])['country_code'].apply(list).to_dict()
    urban_country_code_labels = pd.DataFrame({'urban_id': urban_country_dict.keys(), 'urban_country_code': urban_country_dict.values()})
    urban_country_code_labels['urban_country_code'] = ['–'.join(map(str, l)) for l in urban_country_code_labels['urban_country_code']]
    
    # Join country conurbation names
    full_xwalk = full_xwalk.merge(conurbation_country_name_labels, how = 'left', left_on=['conurbation_id'], right_on = ['conurbation_id'])
    full_xwalk = full_xwalk.merge(conurbation_country_code_labels, how = 'left', left_on=['conurbation_id'], right_on = ['conurbation_id'])
    full_xwalk = full_xwalk.merge(urban_country_name_labels, how = 'left', left_on=['urban_id'], right_on = ['urban_id'])
    full_xwalk = full_xwalk.merge(urban_country_code_labels, how = 'left', left_on=['urban_id'], right_on = ['urban_id'])

    full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name'].str.split('-', n = 2, expand = True)[[0,1]].dropna().astype(str).apply('-'.join, 1)
    full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name_short'].fillna(full_xwalk['conurbation_area_name'])
    
    # Fill in for primary country fields
    full_xwalk['conurbation_country_name'] = full_xwalk['conurbation_country_name'].fillna(full_xwalk['country_name'])
    full_xwalk['urban_country_name'] = full_xwalk['urban_country_name'].fillna(full_xwalk['country_name'])
 
    # Join in conurbation ID ranking for labels
    urban_area_rank = full_xwalk[full_xwalk['urban_id'].notnull()].groupby(['urban_id','conurbation_id']).agg({'urban_area': 'sum'}).reset_index()
    urban_area_rank['rank'] = urban_area_rank.groupby('conurbation_id')['urban_area'].rank(method='first', ascending=False)
    full_xwalk = full_xwalk.merge(urban_area_rank[['urban_id','conurbation_id','rank']], how = 'left', left_on = ['urban_id','conurbation_id'], right_on = ['urban_id','conurbation_id'])

    print('Fully merged')
    # 'Core urban', 'Peripheral urban', 'Peri-urban', 'Non-urban'
    conditions = [(full_xwalk['area_type'] == 'Urban') & (full_xwalk['rank'] == 1),
                  (full_xwalk['area_type'] == 'Urban') & (full_xwalk['rank'] > 1),
                  (full_xwalk['area_type'] == 'Peri-urban'),
                  (full_xwalk['area_type'] == 'Non-urban')]
    labels = ['1 - Core urban', '2 - Peripheral urban', '3 - Peri-urban', '4 - Non-urban']
    full_xwalk["class_urban_hierarchy"] = np.select(conditions, labels, default='4 - Non-urban')
    
    # 'Urban', 'Peri-urban', 'Non-urban'
    conditions = [full_xwalk["class_urban_hierarchy"].isin(['1 - Core urban', '2 - Peripheral urban']),
                  full_xwalk["class_urban_hierarchy"].isin(['3 - Peri-urban']),
                  full_xwalk["class_urban_hierarchy"].isin(['4 - Non-urban'])]
    labels = ['1 - Core & peripheral urban', '2 - Peri-urban', '3 - Non-urban']
    full_xwalk["class_urban_periurban_nonurban"] = np.select(conditions, labels, default='4 - Non-urban')
    
    # 'Urban', 'Non-urban'
    conditions = [full_xwalk["class_urban_hierarchy"].isin(['1 - Core urban', '2 - Peripheral urban', '3 - Peri-urban']),
                  full_xwalk["class_urban_hierarchy"].isin(['4 - Non-urban'])]
    labels = ['1 - Core, peripheral, & peri-urban', '2 - Non-urban']
    full_xwalk["class_urban_nonurban"] = np.select(conditions, labels, default='4 - Non-urban')

    # Fill in non-urban labels
    full_xwalk['urban_center_name'] = full_xwalk['urban_center_name'].fillna('Rest of ' + full_xwalk['country_name'])
    full_xwalk['urban_country_name'] = full_xwalk['urban_country_name'].fillna(full_xwalk['country_name'])
    full_xwalk['urban_country_code'] = full_xwalk['urban_country_code'].fillna(full_xwalk['country_code'])
    full_xwalk['conurbation_area_name'] = full_xwalk['conurbation_area_name'].fillna('Rest of ' + full_xwalk['country_name'])
    full_xwalk['conurbation_country_name'] = full_xwalk['conurbation_country_name'].fillna(full_xwalk['country_name'])
    full_xwalk['conurbation_country_code'] = full_xwalk['conurbation_country_code'].fillna(full_xwalk['country_code'])
    full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name_short'].fillna('Rest of ' + full_xwalk['country_name'])
    
    # Create special aggregation labels
    full_xwalk['c4_code'] = full_xwalk['class_urban_hierarchy'].str.lower().str.split(' - ', n = 2, expand = True)[[0]] 
    full_xwalk['c3_code'] = full_xwalk['class_urban_periurban_nonurban'].str.lower().str.split(' - ', n = 2, expand = True)[[0]] 
    full_xwalk['c2_code'] = full_xwalk['class_urban_nonurban'].str.lower().str.split(' - ', n = 2, expand = True)[[0]] 
    
    full_xwalk['c4_label'] = ' (' + full_xwalk['class_urban_hierarchy'].str.lower().str.split(' - ', n = 2, expand = True)[[1]] + ')'
    full_xwalk['c3_label'] = ' (' + full_xwalk['class_urban_periurban_nonurban'].str.lower().str.split(' - ', n = 2, expand = True)[[1]] + ')'
    full_xwalk['c2_label'] = ' (' + full_xwalk['class_urban_nonurban'].str.lower().str.split(' - ', n = 2, expand = True)[[1]] + ')'
    
    full_xwalk['country_code_c4'] = full_xwalk['country_code'] + '_' + full_xwalk['c4_code']
    full_xwalk['country_code_c3'] = full_xwalk['country_code'] + '_' + full_xwalk['c3_code']
    full_xwalk['country_code_c2'] = full_xwalk['country_code'] + '_' + full_xwalk['c2_code']
    
    full_xwalk['urban_id_c4'] = full_xwalk['urban_id'] + '_' + full_xwalk['c4_code']
    full_xwalk['urban_id_c3'] = full_xwalk['urban_id'] + '_' + full_xwalk['c3_code']
    full_xwalk['urban_id_c2'] = full_xwalk['urban_id'] + '_' + full_xwalk['c2_code']
    
    full_xwalk['conurbation_id_c4'] = full_xwalk['conurbation_id']  + '_' + full_xwalk['c4_code']
    full_xwalk['conurbation_id_c3'] = full_xwalk['conurbation_id']  + '_' + full_xwalk['c3_code']
    full_xwalk['conurbation_id_c2'] = full_xwalk['conurbation_id']  + '_' + full_xwalk['c2_code']
    
    full_xwalk['country_area_class_urban_hierarchy'] = full_xwalk['country_name'] + full_xwalk['c4_label']
    full_xwalk['country_area_class_urban_periurban_nonurban'] = full_xwalk['country_name'] + full_xwalk['c3_label']
    full_xwalk['country_area_class_urban_nonurban'] = full_xwalk['country_name'] + full_xwalk['c2_label']
    
    full_xwalk['urban_center_class_urban_hierarchy'] = full_xwalk['urban_center_name'] + full_xwalk['c4_label']
    full_xwalk['urban_center_class_urban_periurban_nonurban'] = full_xwalk['urban_center_name'] + full_xwalk['c3_label']
    full_xwalk['urban_center_class_urban_nonurban'] = full_xwalk['urban_center_name'] + full_xwalk['c2_label']
    
    full_xwalk['conurbation_class_urban_hierarchy'] = full_xwalk['conurbation_area_name_short'] + full_xwalk['c4_label']
    full_xwalk['conurbation_class_urban_periurban_nonurban'] = full_xwalk['conurbation_area_name_short'] + full_xwalk['c3_label']
    full_xwalk['conurbation_class_urban_nonurban'] = full_xwalk['conurbation_area_name_short'] + full_xwalk['c2_label']

    full_xwalk = full_xwalk[['block_id', 'block_geohash', 'gadm_code', 'country_code', 'country_name', 'continent', 'area_type', 'class_urban_hierarchy', 'class_urban_periurban_nonurban', 'class_urban_nonurban', 'urban_id', 'urban_center_name', 'urban_country_code', 'urban_country_name', 'conurbation_id', 'conurbation_area_name', 'conurbation_area_name_short', 'conurbation_country_code', 'conurbation_country_name', 'c4_code', 'country_code_c4', 'country_area_class_urban_hierarchy', 'urban_id_c4', 'urban_center_class_urban_hierarchy', 'conurbation_id_c4', 'conurbation_class_urban_hierarchy', 'c3_code', 'country_code_c3', 'country_area_class_urban_periurban_nonurban', 'urban_id_c3', 'urban_center_class_urban_periurban_nonurban', 'conurbation_id_c3', 'conurbation_class_urban_periurban_nonurban', 'c2_code', 'country_code_c2', 'country_area_class_urban_nonurban', 'urban_id_c2', 'urban_center_class_urban_nonurban', 'conurbation_id_c2', 'conurbation_class_urban_nonurban', 'urban_regional_layer_code', 'urban_layer_code']]

    assert full_xwalk[full_xwalk['block_id'].duplicated()].shape[0] == 0
    assert full_xwalk[full_xwalk['block_id'].isnull()].shape[0] == 0

    # Write file to parquet CSV
    print('Writing files')
    full_xwalk.to_parquet(path = Path(output_dir) / f'ghsl_crosswalk.parquet')
    full_xwalk.to_parquet(Path(partition_dir), partition_cols=['country_code'])
    print('Finished')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute crosswalk.')
    parser.add_argument('--ghsl_file', required=True, type=Path, dest="ghsl_file", help="Path to GHSL file.")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to outputs directory.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))
