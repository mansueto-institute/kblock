
import pygeos
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
gpd.options.use_pygeos = True

import re
import os
import shutil
from typing import Union
from pathlib import Path
import argparse
import pyarrow
import dask
import dask.dataframe 


def main(ghsl_dir: Path, blocks_dir: Path, output_dir: Path):

    dask_dir =  str(output_dir) + '/dask'
    partition_dir =  str(output_dir) + '/partitions'
    if os.path.isdir(dask_dir): shutil.rmtree(dask_dir)
    if os.path.isdir(partition_dir): shutil.rmtree(partition_dir)
    Path(dask_dir).mkdir(parents=True, exist_ok=True)
    Path(partition_dir).mkdir(parents=True, exist_ok=True)

    country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}
    country_names = pd.DataFrame.from_dict(country_dict).T.reset_index().rename(columns= {'index':'country_code', 0:'country',1:'continent',2:'country_name'}).drop(['country'], axis=1)

    # Setup data
    ghsl_data = gpd.read_file(Path(ghsl_dir))
    ghsl_data = ghsl_data[ghsl_data['GRGN_L2'].isin(['Western Africa', 'Northern Africa', 'Middle Africa', 'Southern Africa', 'Eastern Africa'])][['ID_HDC_G0','GCPNT_LAT', 'GCPNT_LON', 'XBRDR', 'XCTR_NBR', 'CTR_MN_ISO', 'CTR_MN_NM', 'UC_NM_MN','UC_NM_LST','XC_NM_LST','XC_ISO_LST','GRGN_L1','GRGN_L2','geometry']].reset_index()
    ghsl_data = ghsl_data.to_crs(3395)
    ghsl_data['urban_area'] = ghsl_data['geometry'].area
    ghsl_data['urban_radius'] = pygeos.minimum_bounding_radius(pygeos.from_shapely(ghsl_data['geometry']))
    col_recode = {'ID_HDC_G0':'urban_id', 'GCPNT_LAT': 'latitude' , 'GCPNT_LON': 'longitude', 'XBRDR':'crossborder', 'XCTR_NBR':'num_country_crossborder', 'CTR_MN_ISO':'primary_country_code', 'CTR_MN_NM':'primary_country_name', 'UC_NM_MN':'urban_center_name', 'UC_NM_LST':'list_of_urban_center_names', 'XC_NM_LST':'intersected_countries', 'XC_ISO_LST':'intersected_codes', 'GRGN_L1':'geographical_region', 'GRGN_L2':'geographical_subregion'}
    ghsl_data = ghsl_data.rename(columns=col_recode)
    ghsl_data = ghsl_data[['urban_id','urban_center_name','list_of_urban_center_names','primary_country_code','primary_country_name','geographical_subregion','geographical_region','urban_area','urban_radius','latitude', 'longitude', 'geometry']]
    ghsl_data = ghsl_data[ghsl_data['urban_center_name'] != 'N/A']
    ghsl_data = ghsl_data.explode(index_parts=False)
    assert ghsl_data[ghsl_data['urban_id'].duplicated()].shape[0] == 0
    
    # Urban centers
    urban_centers = ghsl_data[['urban_id','geometry']].to_crs(4326)
    
    # Periurban zones based on minimum bounding circle
    periurban_inner = pygeos.from_shapely(ghsl_data['geometry'].buffer(ghsl_data['urban_radius']*1))
    periurban_inner = pygeos.get_parts(pygeos.union_all(periurban_inner))
    periurban_inner = gpd.GeoDataFrame.from_dict({'geometry': pygeos.to_shapely(periurban_inner)}).set_crs(epsg=3395)
    periurban_inner = periurban_inner.assign(conurbation_id = [str(x) for x in list(periurban_inner.index)])
    periurban_inner = periurban_inner[['conurbation_id','geometry']].to_crs(4326)
    
    periurban_residual = periurban_inner.overlay(urban_centers, how='difference')
    urban_centers = pd.concat([urban_centers,periurban_residual])
    
    # Create list of countries based on block directory
    input_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    country_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_file_list]
    full_xwalk = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'urban_id': pd.Series(dtype='float'), 'conurbation_id': pd.Series(dtype='float')})

    for country_code in country_list:
        print(country_code)
    
        blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
    
        # Urban centers 
        urban_overlay = gpd.overlay(df1 = blocks, df2 = urban_centers, how='intersection')
        urban_overlay['area'] = urban_overlay['geometry'].to_crs(3395).area
        urban_overlay['rank'] = urban_overlay.groupby('block_id')['area'].rank(method='first', ascending=False)
        urban_overlay = urban_overlay[(urban_overlay['urban_id'].notnull()) & (urban_overlay['rank'] == 1)]
        urban_overlay = urban_overlay[['block_id','urban_id','rank']]
    
        blocks_xwalk = blocks.merge(urban_overlay, how = 'left', left_on=['block_id'], right_on = ['block_id'])
        blocks_xwalk = blocks_xwalk[['block_id','gadm_code','country_code','geometry','urban_id']]
    
        # Conurbations
        conurbation_overlay = gpd.overlay(df1 = blocks, df2 = periurban_inner, how='intersection')
        conurbation_overlay['area'] = conurbation_overlay['geometry'].to_crs(3395).area
        conurbation_overlay['rank'] = conurbation_overlay.groupby('block_id')['area'].rank(method='first', ascending=False)
        conurbation_overlay = conurbation_overlay[(conurbation_overlay['conurbation_id'].notnull()) & (conurbation_overlay['rank'] == 1)]
        conurbation_overlay = conurbation_overlay[['block_id','conurbation_id','rank']]
    
        blocks_xwalk = blocks_xwalk.merge(conurbation_overlay, how = 'left', left_on=['block_id'], right_on = ['block_id'])
        blocks_xwalk = blocks_xwalk[['block_id','gadm_code','country_code','geometry','urban_id','conurbation_id']]
    
        if blocks_xwalk[blocks_xwalk['block_id'].duplicated()].shape[0] > 0:
            print('Duplicated rows')
    
        # Blocks crosswalk
        blocks_xwalk = blocks_xwalk[['block_id','gadm_code','country_code','urban_id','conurbation_id']]
    
        blocks_xwalk = dask.dataframe.from_pandas(data = blocks_xwalk, npartitions = 1) 
        dask.dataframe.to_parquet(df = blocks_xwalk, path = Path(dask_dir) / f'{country_code}.parquet', engine='pyarrow', compression='snappy', append=True, ignore_divisions=True)

    # Read parquet data
    print('Reading countries')
    full_xwalk = dask.dataframe.read_parquet(path = Path(dask_dir)).compute()
    
    # Join GHSL data to crosswalk
    full_xwalk = full_xwalk.merge(ghsl_data[['urban_id','urban_center_name','urban_radius','urban_area']], how = 'left', left_on=['urban_id'], right_on = ['urban_id'])
    
    # Coerce all conurbations to have at least one urban_id (applies to blocks on fringes of )
    full_xwalk['urban_conurbation_count'] = full_xwalk.groupby(['conurbation_id'])['urban_id'].transform('count')
    full_xwalk.loc[full_xwalk['urban_conurbation_count'] == 0, 'conurbation_id'] = np.nan
    
    # Create conurbation names
    conurbation_dict = full_xwalk[(full_xwalk['conurbation_id'].notnull()) & (full_xwalk['urban_center_name'].notnull())][['conurbation_id','urban_center_name','urban_area']].drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['conurbation_id'])['urban_center_name'].apply(list).to_dict()
    conurbation_labels = pd.DataFrame({'conurbation_id': conurbation_dict.keys(), "conurbation_area_name": conurbation_dict.values()})
    conurbation_labels['conurbation_area_name'] = [' — '.join(map(str, l)) for l in conurbation_labels ['conurbation_area_name']]
    
    # Join in country and conurbation names to crosswalk
    full_xwalk = full_xwalk.merge(conurbation_labels, how = 'left', left_on=['conurbation_id'], right_on = ['conurbation_id'])
    full_xwalk = full_xwalk.merge(country_names, how = 'left', left_on=['country_code'], right_on = ['country_code'])
    
    # Create country conurbation names
    conurbation_country_dict = full_xwalk.groupby(['conurbation_id','country_name']).agg({'urban_area': 'sum'}).reset_index().drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['conurbation_id'])['country_name'].apply(list).to_dict()
    conurbation_country_labels = pd.DataFrame({'conurbation_id': conurbation_country_dict.keys(), "conurbation_country_name": conurbation_country_dict.values()})
    conurbation_country_labels["conurbation_country_name"] = [' - '.join(map(str, l)) for l in conurbation_country_labels["conurbation_country_name"]]
    
    # Create country urban names
    urban_country_dict = full_xwalk.groupby(['urban_id','country_name']).agg({'urban_area': 'sum'}).reset_index().drop_duplicates().reset_index(drop = True).sort_values(by=['urban_area'], ascending = False).groupby(['urban_id'])['country_name'].apply(list).to_dict()
    urban_country_labels = pd.DataFrame({'urban_id': urban_country_dict.keys(), 'urban_country_name': urban_country_dict.values()})
    urban_country_labels['urban_country_name'] = [' - '.join(map(str, l)) for l in urban_country_labels['urban_country_name']]
    
    # Join country conurbation names
    full_xwalk = full_xwalk.merge(conurbation_country_labels, how = 'left', left_on=['conurbation_id'], right_on = ['conurbation_id'])
    full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name'].str.split(' — ', n = 2, expand = True)[[0,1]].dropna().astype(str).apply(' — '.join, 1)
    full_xwalk = full_xwalk.merge(urban_country_labels, how = 'left', left_on=['urban_id'], right_on = ['urban_id'])
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
    conditions = [(full_xwalk['urban_id'].notnull()) & (full_xwalk['rank'] == 1),
                  (full_xwalk['urban_id'].notnull()) & (full_xwalk['rank'] > 1),
                  (full_xwalk['urban_id'].isnull()) & (full_xwalk['conurbation_id'].notnull()),
                  (full_xwalk['urban_id'].isnull()) & (full_xwalk['conurbation_id'].isnull())]
    labels = ['Core urban', 'Peripheral urban', 'Peri-urban', 'Non-urban']
    full_xwalk["area_type_4"] = np.select(conditions, labels, default='Non-urban')
    
    # 'Urban', 'Peri-urban', 'Non-urban'
    conditions = [full_xwalk["area_type_4"].isin(['Core urban', 'Peripheral urban']),
                  full_xwalk["area_type_4"].isin(['Peri-urban']),
                  full_xwalk["area_type_4"].isin(['Non-urban'])]
    labels = ['Core & peripheral urban', 'Peri-urban', 'Non-urban']
    full_xwalk["area_type_3"] = np.select(conditions, labels, default='Non-urban')
    
    # 'Urban', 'Non-urban'
    conditions = [full_xwalk["area_type_4"].isin(['Core urban', 'Peripheral urban', 'Peri-urban']),
                  full_xwalk["area_type_4"].isin(['Non-urban'])]
    labels = ['Core, peripheral, & peri-urban', 'Non-urban']
    full_xwalk["area_type_2"] = np.select(conditions, labels, default='Non-urban')
    
    full_xwalk = full_xwalk[['block_id', 'area_type_2', 'area_type_3', 'area_type_4', 'gadm_code', 'country_code', 'country_name', 'continent', 'urban_id', 'urban_center_name', 'urban_country_name', 'conurbation_id', 'conurbation_area_name', 'conurbation_area_name_short', 'conurbation_country_name']]
    
    # Fill urban using rest of country
    full_xwalk['urban_id'] = full_xwalk['urban_id'].fillna('r_' + full_xwalk['country_code'])
    full_xwalk['urban_center_name'] = full_xwalk['urban_center_name'].fillna('Rest of ' + full_xwalk['country_name'])
    full_xwalk['urban_country_name'] = full_xwalk['urban_country_name'].fillna(full_xwalk['country_name'])
    
    # Fill conurbation with rest of country
    full_xwalk['conurbation_id'] = full_xwalk['conurbation_id'].fillna('r_' + full_xwalk['country_code'])
    full_xwalk['conurbation_area_name'] = full_xwalk['conurbation_area_name'].fillna('Rest of ' + full_xwalk['country_name'])
    full_xwalk['conurbation_country_name'] = full_xwalk['conurbation_country_name'].fillna(full_xwalk['country_name'])
    full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name_short'].fillna('Rest of ' + full_xwalk['country_name'])
    
    full_xwalk['urban_id'] = full_xwalk['urban_id'].astype("string")
    full_xwalk['conurbation_id'] = full_xwalk['conurbation_id'].astype("string")
    
    # Write file to parquet CSV
    print('Writing files')
    full_xwalk.to_parquet(path = Path(output_dir) / f'ghsl_crosswalk.parquet')
    full_xwalk.to_parquet(Path(partition_dir), partition_cols=['country_code'])
    print('Finished')

def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute crosswalk.')
    parser.add_argument('--ghsl_dir', required=True, type=Path, dest="ghsl_dir", help="Path to GHSL directory")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to outputs directory")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))

