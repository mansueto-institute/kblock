
import numpy as np
import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = True

import re
import os
import logging
from pathlib import Path
import argparse
import dask
import dask.dataframe
import dask_geopandas
import warnings; warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

def main(log_file: Path, country_chunk: list, africapolis_file: Path, ghsl_file: Path, gadm_dir: Path, blocks_dir: Path, output_dir: Path):

    logging.basicConfig(filename=Path(log_file), format='%(asctime)s:%(message)s: ', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    # Create directories
    crosswalk_dir =  str(output_dir) + '/crosswalks'
    Path(crosswalk_dir).mkdir(parents=True, exist_ok=True)
    dask_dir =  str(output_dir) +  '/crosswalks' +  '/dask'
    if os.path.isdir(dask_dir) & os.path.isdir(Path(dask_dir) / f'crosswalks.parquet'):
        finished_countries_list = dask.dataframe.read_parquet(path = Path(dask_dir)).compute()
        finished_countries_list = finished_countries_list['country_code'].unique()
        if os.path.isfile(Path(crosswalk_dir) / f'ghsl_crosswalk.parquet'):
            finished_countries_list_2 = pd.read_parquet(Path(crosswalk_dir) / f'ghsl_crosswalk.parquet')
            finished_countries_list_2 = finished_countries_list_2['country_code'].unique()
            finished_countries_list = [x for x in finished_countries_list if x in set(finished_countries_list_2)]
    else:
        Path(dask_dir).mkdir(parents=True, exist_ok=True)
        finished_countries_list = []

    # Create list of countries based on block directory
    input_file_list = list(filter(re.compile("blocks_").match, sorted(list(os.listdir(Path(blocks_dir))))))
    input_file_list = [(re.sub('blocks_', '', re.sub('.parquet', '', i))) for i in input_file_list]

    output_file_list = list(filter(re.compile("block_xwalk_").match, sorted(list(os.listdir(Path(dask_dir))))))
    output_file_list = [(re.sub('block_xwalk_', '', re.sub('.gpkg', '', i))) for i in output_file_list]

    # Country chunk list
    if country_chunk:
        country_list = [x for x in country_chunk if x in set(input_file_list)]
        country_list = [x for x in country_list if x not in set(finished_countries_list)]
        country_list_regions = [x for x in country_chunk if x in set(input_file_list)]
        country_list_regions = [x for x in country_list_regions if x not in set(output_file_list)]
        if len(country_list) == 0 and len(country_list_regions) == 0 and os.path.isfile(Path(crosswalk_dir) / f'urban_boundaries.gpkg'):
            raise ValueError('No more files to process.')
    else:
        raise ValueError('Empty country_chunk arg.')

    logging.info(f'Country list: {country_list}')

    # Country labels
    country_dict = {'DZA' : ['algeria', 'Africa', 'Algeria'], 'AGO' : ['angola', 'Africa', 'Angola'], 'BEN' : ['benin', 'Africa', 'Benin'], 'BWA' : ['botswana', 'Africa', 'Botswana'], 'BFA' : ['burkina-faso', 'Africa', 'Burkina Faso'], 'BDI' : ['burundi', 'Africa', 'Burundi'], 'CPV' : ['cape-verde', 'Africa', 'Cabo Verde'], 'CMR' : ['cameroon', 'Africa', 'Cameroon'], 'CAF' : ['central-african-republic', 'Africa', 'Central African Republic'], 'TCD' : ['chad', 'Africa', 'Chad'], 'COM' : ['comores', 'Africa', 'Comoros'], 'COG' : ['congo-brazzaville', 'Africa', 'Congo'], 'CIV' : ['ivory-coast', 'Africa', "Côte d'Ivoire"], 'COD' : ['congo-democratic-republic', 'Africa', 'Democratic Republic of the Congo '], 'DJI' : ['djibouti', 'Africa', 'Djibouti'], 'EGY' : ['egypt', 'Africa', 'Egypt'], 'GNQ' : ['equatorial-guinea', 'Africa', 'Equatorial Guinea'], 'ERI' : ['eritrea', 'Africa', 'Eritrea'], 'SWZ' : ['swaziland', 'Africa', 'Eswatini'], 'ETH' : ['ethiopia', 'Africa', 'Ethiopia'], 'GAB' : ['gabon', 'Africa', 'Gabon'], 'GMB' : ['senegal-and-gambia', 'Africa', 'Gambia'], 'GHA' : ['ghana', 'Africa', 'Ghana'], 'GIN' : ['guinea', 'Africa', 'Guinea'], 'GNB' : ['guinea-bissau', 'Africa', 'Guinea-Bissau'], 'KEN' : ['kenya', 'Africa', 'Kenya'], 'LSO' : ['lesotho', 'Africa', 'Lesotho'], 'LBR' : ['liberia', 'Africa', 'Liberia'], 'LBY' : ['libya', 'Africa', 'Libya'], 'MDG' : ['madagascar', 'Africa', 'Madagascar'], 'MWI' : ['malawi', 'Africa', 'Malawi'], 'MLI' : ['mali', 'Africa', 'Mali'], 'MRT' : ['mauritania', 'Africa', 'Mauritania'], 'MUS' : ['mauritius', 'Africa', 'Mauritius'], 'MAR' : ['morocco', 'Africa', 'Morocco'], 'ESH' : ['morocco', 'Africa', 'Morocco'], 'MOZ' : ['mozambique', 'Africa', 'Mozambique'], 'NAM' : ['namibia', 'Africa', 'Namibia'], 'NER' : ['niger', 'Africa', 'Niger'], 'NGA' : ['nigeria', 'Africa', 'Nigeria'], 'RWA' : ['rwanda', 'Africa', 'Rwanda'], 'SHN' : ['saint-helena-ascension-and-tristan-da-cunha', 'Africa', 'Saint Helena'], 'STP' : ['sao-tome-and-principe', 'Africa', 'Sao Tome and Principe'], 'SEN' : ['senegal-and-gambia', 'Africa', 'Senegal'], 'SYC' : ['seychelles', 'Africa', 'Seychelles'], 'SLE' : ['sierra-leone', 'Africa', 'Sierra Leone'], 'SOM' : ['somalia', 'Africa', 'Somalia'], 'ZAF' : ['south-africa', 'Africa', 'South Africa'], 'SSD' : ['south-sudan', 'Africa', 'South Sudan'], 'SDN' : ['sudan', 'Africa', 'Sudan'], 'TZA' : ['tanzania', 'Africa', 'Tanzania'], 'TGO' : ['togo', 'Africa', 'Togo'], 'TUN' : ['tunisia', 'Africa', 'Tunisia'], 'UGA' : ['uganda', 'Africa', 'Uganda'], 'ZMB' : ['zambia', 'Africa', 'Zambia'], 'ZWE' : ['zimbabwe', 'Africa', 'Zimbabwe'], 'AFG' : ['afghanistan', 'Asia', 'Afghanistan'], 'ARM' : ['armenia', 'Asia', 'Armenia'], 'AZE' : ['azerbaijan', 'Asia', 'Azerbaijan'], 'BHR' : ['gcc-states', 'Asia', 'Bahrain'], 'BGD' : ['bangladesh', 'Asia', 'Bangladesh'], 'BTN' : ['bhutan', 'Asia', 'Bhutan'], 'BRN' : ['malaysia-singapore-brunei', 'Asia', 'Brunei Darussalam'], 'MYS' : ['malaysia-singapore-brunei', 'Asia', 'Malaysia'], 'SGP' : ['malaysia-singapore-brunei', 'Asia', 'Singapore'], 'KHM' : ['cambodia', 'Asia', 'Cambodia'], 'CHN' : ['china', 'Asia', 'China'], 'IND' : ['india', 'Asia', 'India'], 'IDN' : ['indonesia', 'Asia', 'Indonesia'], 'IRQ' : ['iraq', 'Asia', 'Iraq'], 'IRN' : ['iran', 'Asia', 'Islamic Republic of Iran'], 'PSE' : ['israel-and-palestine', 'Asia', 'Palestine'], 'ISR' : ['israel-and-palestine', 'Asia', 'Israel'], 'JPN' : ['japan', 'Asia', 'Japan'], 'JOR' : ['jordan', 'Asia', 'Jordan'], 'KAZ' : ['kazakhstan', 'Asia', 'Kazakhstan'], 'KGZ' : ['kyrgyzstan', 'Asia', 'Kyrgyzstan'], 'LAO' : ['laos', 'Asia', 'Laos'], 'LBN' : ['lebanon', 'Asia', 'Lebanon'], 'MDV' : ['maldives', 'Asia', 'Maldives'], 'MNG' : ['mongolia', 'Asia', 'Mongolia'], 'MMR' : ['myanmar', 'Asia', 'Myanmar'], 'NPL' : ['nepal', 'Asia', 'Nepal'], 'PRK' : ['north-korea', 'Asia', 'North Korea'], 'PAK' : ['pakistan', 'Asia', 'Pakistan'], 'PHL' : ['philippines', 'Asia', 'Philippines'], 'KOR' : ['south-korea', 'Asia', 'South Korea'], 'LKA' : ['sri-lanka', 'Asia', 'Sri Lanka'], 'SYR' : ['syria', 'Asia', 'Syrian Arab Republic'], 'TWN' : ['taiwan', 'Asia', 'Taiwan'], 'TJK' : ['tajikistan', 'Asia', 'Tajikistan'], 'THA' : ['thailand', 'Asia', 'Thailand'], 'TKM' : ['turkmenistan', 'Asia', 'Turkmenistan'], 'UZB' : ['uzbekistan', 'Asia', 'Uzbekistan'], 'VNM' : ['vietnam', 'Asia', 'Vietnam'], 'YEM' : ['yemen', 'Asia', 'Yemen'], 'ALB' : ['albania', 'Europe', 'Albania'], 'AND' : ['andorra', 'Europe', 'Andorra'], 'AUT' : ['austria', 'Europe', 'Austria'], 'BLR' : ['belarus', 'Europe', 'Belarus'], 'BEL' : ['belgium', 'Europe', 'Belgium'], 'BIH' : ['bosnia-herzegovina', 'Europe', 'Bosnia and Herzegovina'], 'BGR' : ['bulgaria', 'Europe', 'Bulgaria'], 'HRV' : ['croatia', 'Europe', 'Croatia'], 'CYP' : ['cyprus', 'Europe', 'Cyprus'], 'CZE' : ['czech-republic', 'Europe', 'Czechia'], 'DNK' : ['denmark', 'Europe', 'Denmark'], 'EST' : ['estonia', 'Europe', 'Estonia'], 'FRO' : ['faroe-islands', 'Europe', 'Faroe Islands'], 'FIN' : ['finland', 'Europe', 'Finland'], 'FRA' : ['france', 'Europe', 'France'], 'GEO' : ['georgia', 'Europe', 'Georgia'], 'DEU' : ['germany', 'Europe', 'Germany'], 'GRC' : ['greece', 'Europe', 'Greece'], 'HUN' : ['hungary', 'Europe', 'Hungary'], 'ISL' : ['iceland', 'Europe', 'Iceland'], 'IRL' : ['ireland-and-northern-ireland', 'Europe', 'Ireland'], 'IMN' : ['isle-of-man', 'Europe', 'Isle of Man'], 'ITA' : ['italy', 'Europe', 'Italy'], 'KOS' : ['kosovo', 'Europe', 'Kosovo'], 'LVA' : ['latvia', 'Europe', 'Latvia'], 'LIE' : ['liechtenstein', 'Europe', 'Liechtenstein'], 'LTU' : ['lithuania', 'Europe', 'Lithuania'], 'LUX' : ['luxembourg', 'Europe', 'Luxembourg'], 'MLT' : ['malta', 'Europe', 'Malta'], 'MDA' : ['moldova', 'Europe', 'Moldova'], 'MCO' : ['monaco', 'Europe', 'Monaco'], 'MNE' : ['montenegro', 'Europe', 'Montenegro'], 'NLD' : ['netherlands', 'Europe', 'Netherlands'], 'MKD' : ['macedonia', 'Europe', 'North Macedonia'], 'NOR' : ['norway', 'Europe', 'Norway'], 'POL' : ['poland', 'Europe', 'Poland'], 'PRT' : ['portugal', 'Europe', 'Portugal'], 'ROU' : ['romania', 'Europe', 'Romania'], 'RUS' : ['russia', 'Europe', 'Russian Federation'], 'SRB' : ['serbia', 'Europe', 'Serbia'], 'SVK' : ['slovakia', 'Europe', 'Slovakia'], 'SVN' : ['slovenia', 'Europe', 'Slovenia'], 'ESP' : ['spain', 'Europe', 'Spain'], 'SWE' : ['sweden', 'Europe', 'Sweden'], 'CHE' : ['switzerland', 'Europe', 'Switzerland'], 'TUR' : ['turkey', 'Europe', 'Turkey'], 'UKR' : ['ukraine', 'Europe', 'Ukraine'], 'GBR' : ['great-britain', 'Europe', 'United Kingdom'], 'CAN' : ['canada', 'North America', 'Canada'], 'GRL' : ['greenland', 'North America', 'Greenland'], 'MEX' : ['mexico', 'North America', 'Mexico'], 'USA' : ['usa', 'North America', 'United States of America'], 'AUS' : ['australia', 'Oceania', 'Australia'], 'COK' : ['cook-islands', 'Oceania', 'Cook Islands'], 'FJI' : ['fiji', 'Oceania', 'Fiji'], 'KIR' : ['kiribati', 'Oceania', 'Kiribati'], 'MHL' : ['marshall-islands', 'Oceania', 'Marshall Islands'], 'FSM' : ['micronesia', 'Oceania', 'Micronesia'], 'NRU' : ['nauru', 'Oceania', 'Nauru'], 'NCL' : ['new-caledonia', 'Oceania', 'New Caledonia'], 'NZL' : ['new-zealand', 'Oceania', 'New Zealand'], 'NIU' : ['niue', 'Oceania', 'Niue'], 'PLW' : ['palau', 'Oceania', 'Palau'], 'PNG' : ['papua-new-guinea', 'Oceania', 'Papua New Guinea'], 'WSM' : ['samoa', 'Oceania', 'Samoa'], 'SLB' : ['solomon-islands', 'Oceania', 'Solomon Islands'], 'TON' : ['tonga', 'Oceania', 'Tonga'], 'TUV' : ['tuvalu', 'Oceania', 'Tuvalu'], 'VUT' : ['vanuatu', 'Oceania', 'Vanuatu'], 'ARG' : ['argentina', 'South America', 'Argentina'], 'BOL' : ['bolivia', 'South America', 'Bolivia'], 'BRA' : ['brazil', 'South America', 'Brazil'], 'CHL' : ['chile', 'South America', 'Chile'], 'COL' : ['colombia', 'South America', 'Colombia'], 'ECU' : ['ecuador', 'South America', 'Ecuador'], 'PRY' : ['paraguay', 'South America', 'Paraguay'], 'PER' : ['peru', 'South America', 'Peru'], 'SUR' : ['suriname', 'South America', 'Suriname'], 'URY' : ['uruguay', 'South America', 'Uruguay'], 'VEN' : ['venezuela', 'South America', 'Venezuela']}
    country_names = pd.DataFrame.from_dict(country_dict).T.reset_index().rename(columns= {'index':'country_code', 0:'country',1:'continent',2:'country_name'}).drop(['country'], axis=1)

    # Build crosswalks
    if len(country_list) > 0:

        # Prepare GHSL data
        ghsl_data = gpd.read_file(Path(ghsl_file)).to_crs(4326)
        ghsl_data = ghsl_data[ghsl_data['GRGN_L2'].isin(['Western Africa', 'Northern Africa', 'Middle Africa', 'Southern Africa', 'Eastern Africa'])][['ID_HDC_G0','UC_NM_MN','P15','geometry']].reset_index()
        col_recode = {'ID_HDC_G0':'urban_id', 'UC_NM_MN':'urban_center_name','P15':'population_15'}
        ghsl_data = ghsl_data.rename(columns=col_recode)
        ghsl_data = ghsl_data[['urban_id','urban_center_name','population_15','geometry']]
        ghsl_data = ghsl_data[ghsl_data['urban_center_name'] != 'N/A']
        ghsl_data = ghsl_data.explode(index_parts=False)
        ghsl_data['urban_id'] = ghsl_data['urban_id'].astype('int64').astype('string')
        ghsl_data['urban_id'] = 'ghsl_' + ghsl_data['urban_id']
        ghsl_data['urban_id'][ghsl_data['urban_id'] == 'ghsl_2544'] = 'ghsl_2565' # fix Abuja
        ghsl_data = ghsl_data[['urban_id','urban_center_name','population_15','geometry']].dissolve(by=['urban_id','urban_center_name'], as_index=False)
        assert ghsl_data[ghsl_data['urban_id'].duplicated()].shape[0] == 0

        # Prepare Africapolis data
        africapolis_data = gpd.read_file(Path(africapolis_file)).to_crs(4326)
        africapolis_data = africapolis_data[['agglosID', 'agglosName','Metropole','Pop2015', 'geometry']].rename(columns={'agglosID': 'agglosid','agglosName': 'agglosname','Metropole': 'metropole','Pop2015': 'pop2015'})
        africapolis_data['agglosid'] = africapolis_data['agglosid'].astype('int64').astype('str')
        africapolis_data['agglosid'] = 'africapolis_' + africapolis_data['agglosid']
        africapolis_data['geometry'] = africapolis_data['geometry'].to_crs(3395).simplify(100).to_crs(4326).make_valid()

        # Combine GHSL and Africapolis data
        ghsl_africapolis_data = pd.concat([ghsl_data[['geometry']],africapolis_data[['geometry']]]).reset_index(drop = True)
        ghsl_africapolis_data = dask_geopandas.from_geopandas(ghsl_africapolis_data, npartitions = 50)
        ghsl_africapolis_data = ghsl_africapolis_data.dissolve().compute()
        ghsl_africapolis_data = ghsl_africapolis_data.reset_index(drop=True)
        ghsl_africapolis_data = ghsl_africapolis_data.explode(index_parts=False).reset_index(drop=True)
        ghsl_africapolis_data = ghsl_africapolis_data.assign(contiguous_urban_id = [str(x) for x in list(ghsl_africapolis_data.index)])

        # Define urban centers
        urban_centers = ghsl_data[['urban_id','urban_center_name','geometry']].to_crs(4326)
        urban_centers['area_type'] = 'Urban'
        urban_centers['urban_area'] = urban_centers['geometry'].to_crs(3395).area

        # Create conurbation buffer
        logging.info(f'Create conurbation buffer')
        conurbation_buffer = urban_centers.copy()
        # Set the conurbation buffer to 10 KM or the equivalent circular diameter, whichever is smaller. The periurban buffer also incorporates blocks with most of their area intersecting Africapolis urban delineations.
        # Equivalent Circular Diameter The Diameter of the circle that would have the equivalent area as this object.
        conurbation_buffer['urban_diameter'] = (np.sqrt(conurbation_buffer['urban_area']/np.pi)*2).clip(lower=0, upper=10000)
        conurbation_buffer['geometry'] = conurbation_buffer['geometry'].to_crs(3395).buffer(conurbation_buffer['urban_diameter']*1).to_crs(4326)
        # Set the conurbation to 10km (simple method but connects areas that are too sparse and spread out)
        #conurbation_buffer['geometry'] = conurbation_buffer['geometry'].to_crs(3395).buffer(10000)
        conurbation_buffer = conurbation_buffer[['geometry']].to_crs(4326)
        conurbation_buffer = conurbation_buffer.dissolve(as_index=False)
        conurbation_buffer = conurbation_buffer.explode(index_parts=False).reset_index(drop=True)
        conurbation_buffer = conurbation_buffer.assign(conurbation_id = [str(x) for x in list(conurbation_buffer.index)])
        conurbation_buffer = conurbation_buffer[['conurbation_id','geometry']].to_crs(4326)

        # Create buffer around Africapolis areas that intersect with conurbation buffer
        conurbation_buffer_africapolis = gpd.sjoin(left_df = ghsl_africapolis_data[['geometry']], right_df = conurbation_buffer, predicate = 'intersects', how = 'inner').drop(['index_right'], axis=1)
        conurbation_buffer_africapolis['geometry'] = conurbation_buffer_africapolis['geometry'].to_crs(3395).buffer(1000).to_crs(4326)
        conurbation_buffer_africapolis = dask_geopandas.from_geopandas(conurbation_buffer_africapolis, npartitions = 50)
        conurbation_buffer_africapolis = conurbation_buffer_africapolis.dissolve('conurbation_id').compute()
        conurbation_buffer_africapolis = conurbation_buffer_africapolis.reset_index()

        # Combine conurbation buffer and Africapolis buffer
        conurbation_buffer = pd.concat([conurbation_buffer, conurbation_buffer_africapolis]).reset_index(drop=True)
        conurbation_buffer = dask_geopandas.from_geopandas(conurbation_buffer, npartitions = 50)
        #conurbation_buffer = conurbation_buffer.dissolve('conurbation_id').compute()
        #conurbation_buffer = conurbation_buffer.reset_index()
        conurbation_buffer = conurbation_buffer.dissolve().compute()
        conurbation_buffer = conurbation_buffer.explode(index_parts=False).reset_index(drop=True)
        conurbation_buffer = conurbation_buffer.assign(conurbation_id = [str(x) for x in list(conurbation_buffer.index)])
        conurbation_buffer = conurbation_buffer[['conurbation_id','geometry']].to_crs(4326)

        # Define peri-urban areas (difference between conurban and urban)
        logging.info(f'Create periurban buffer')
        periurban_buffer = gpd.overlay(df1 = conurbation_buffer, df2 = urban_centers, how = 'difference')
        periurban_buffer['urban_id'] = 'periurban_' + periurban_buffer['conurbation_id']
        periurban_buffer['area_type'] = 'Peri-urban'
        urban_periurban_buffer = pd.concat([urban_centers[['urban_id','area_type','geometry']], periurban_buffer[['urban_id','area_type','geometry']]])
        # urban_periurban_buffer.explore('area_type', cmap='viridis')

        logging.info(f'Input data prepared')
        full_xwalk = pd.DataFrame({'block_id': pd.Series(dtype='str'), 'gadm_code': pd.Series(dtype='str'), 'country_code': pd.Series(dtype='str'), 'urban_id': pd.Series(dtype='float'), 'conurbation_id': pd.Series(dtype='float')})

        # Loop through countries
        for country_code in country_list:
            logging.info(f'Processing country: {country_code}')
            print(country_code)

            # Overlay GHSL geometry and 10km buffered geometry on blocks
            logging.info(f'Overlay GHSL geometry and 10km buffered geometry on blocks')
            blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
            urban_periurban_blocks = gpd.overlay(df1 = blocks, df2 = urban_periurban_buffer, how='intersection', keep_geom_type=None, make_valid=True)
            non_urban_blocks = gpd.overlay(df1 = blocks, df2 = urban_periurban_buffer, how='difference', keep_geom_type=None, make_valid=True)
            non_urban_blocks['urban_id'] = 'nonurban_' + non_urban_blocks['country_code']
            non_urban_blocks['area_type'] = 'Non-urban'
            blocks_class = pd.concat([urban_periurban_blocks, non_urban_blocks])

            # Assign blocks to GHSL or periurban based on whichever has greater share of block area
            logging.info(f'Assign blocks to GHSL or periurban based on whichever has greater share of block area')
            blocks_class['area'] = blocks_class['geometry'].to_crs(3395).area
            blocks_class['rank'] = blocks_class.groupby('block_id')['area'].rank(method='first', ascending=False)
            urban_centers_refined = blocks_class[blocks_class['rank'] == 1]
            urban_centers_refined = urban_centers_refined[['block_id','urban_id','area_type']]
            urban_centers_refined = urban_centers_refined[urban_centers_refined['area_type'] == 'Urban']

            # Merge assigned area label and urban_id to blocks (1 to 1 match)
            urban_centers_refined = blocks.merge(urban_centers_refined, how = 'inner', left_on=['block_id'], right_on = ['block_id'])

            # Dissolve blocks into urban areas so GHSL contour matches block outlines
            logging.info(f'Dissolve blocks into urban areas so GHSL contour matches block outlines')
            urban_centers_refined = urban_centers_refined[['urban_id','geometry']]
            urban_centers_refined = dask_geopandas.from_geopandas(urban_centers_refined, npartitions = 50)
            urban_centers_refined = urban_centers_refined.dissolve(by = 'urban_id').compute()
            urban_centers_refined.reset_index(inplace = True)
            urban_centers_refined['geometry'] = urban_centers_refined['geometry'].make_valid()

            # Explode dissolved blocks in urban area and rank order chunks by size (deals with discontguous areas)
            logging.info(f'Explode dissolved blocks in urban area and rank order chunks by size')
            urban_centers_refined = urban_centers_refined.explode(index_parts=False).reset_index(drop=True)
            urban_centers_refined = urban_centers_refined.reset_index(drop = False)
            urban_centers_refined = urban_centers_refined[urban_centers_refined['geometry'].geom_type == 'Polygon']
            urban_centers_refined['area'] = urban_centers_refined['geometry'].to_crs(3395).area
            urban_centers_refined['area_total'] = urban_centers_refined.groupby(['urban_id','index'])['area'].transform('sum')
            urban_centers_refined['share'] = urban_centers_refined['area']/urban_centers_refined['area_total']
            urban_centers_refined['rank'] = urban_centers_refined.groupby('urban_id')['area'].rank(method='first', ascending=False)

            # Join in original GHSL label on to contiguous exploded block groups
            logging.info(f'Join in original GHSL label on to contiguous exploded block groups')
            print(1)
            urban_flag = urban_centers.copy()
            urban_flag = urban_flag.assign(urban_center=1)
            urban_flag = urban_flag[['urban_center','geometry']]
            urban_centers_refined = gpd.sjoin(left_df = urban_centers_refined, right_df = urban_flag, predicate = 'intersects', how = 'left').drop(['index_right'], axis=1)

            # Process only if there is a GHSL urban area intersecting country
            if urban_centers_refined.shape[0] > 0:

                # Keep only largest contiguous blocks for each GHSL ID (keeps largest contiguous chunk or in urban area and over 1km square)
                logging.info(f'Keep only largest contiguous blocks for each GHSL ID')
                print(2)
                urban_centers_refined = urban_centers_refined[((urban_centers_refined['urban_center'] == 1) & (urban_centers_refined['area'] >= 1000000) & (urban_centers_refined['share'] == 1)) | (urban_centers_refined['rank'] == 1)]
                urban_centers_refined = urban_centers_refined[['urban_id','geometry']]
                urban_centers_refined = dask_geopandas.from_geopandas(urban_centers_refined, npartitions = 50)
                urban_centers_refined = urban_centers_refined.dissolve(by = 'urban_id').compute()
                urban_centers_refined.reset_index(inplace = True)
                urban_centers_refined['area_type'] = 'Urban'

                # Difference the urban area from the conurbation buffered area
                logging.info(f'Difference the urban area from the conurbation buffered area')
                print(3)
                target_conurbations = conurbation_buffer.sjoin(urban_centers_refined, how = 'inner', predicate = 'intersects')[['conurbation_id']]
                target_conurbations = conurbation_buffer[conurbation_buffer['conurbation_id'].isin(target_conurbations['conurbation_id'])]
                periurban_buffer_refined = gpd.overlay(df1 = target_conurbations, df2 = urban_centers_refined, how = 'difference')

                periurban_buffer_refined = periurban_buffer_refined.explode(index_parts=False)
                periurban_buffer_refined['area'] = periurban_buffer_refined['geometry'].to_crs(3395).area
                periurban_buffer_refined['rank'] = periurban_buffer_refined.groupby('conurbation_id')['area'].rank(method='first', ascending=False)
                periurban_buffer_refined = periurban_buffer_refined[periurban_buffer_refined['rank'] == 1][['conurbation_id', 'geometry']]

                periurban_buffer_refined['urban_id'] = 'periurban_' + periurban_buffer_refined['conurbation_id']
                periurban_buffer_refined['area_type'] = 'Peri-urban'
                periurban_buffer_refined = periurban_buffer_refined[['urban_id','area_type','geometry']]

                # Difference periurban areas from target conurbations to get a clean polygon of urban areas without holes
                logging.info(f'Difference periurban areas from target conurbations to get a clean polygon')
                print(4)
                urban_buffer_refined = gpd.overlay(df1 = target_conurbations, df2 = periurban_buffer_refined, how = 'difference')
                urban_buffer_refined = urban_buffer_refined.explode(index_parts=False)
                urban_buffer_refined = urban_buffer_refined.sjoin(urban_centers_refined, how = 'left', predicate = 'intersects')
                urban_buffer_refined = urban_buffer_refined[['urban_id', 'area_type', 'geometry']].reset_index(drop = True)
                urban_buffer_refined = urban_buffer_refined.dissolve(by = ['urban_id','area_type'], as_index=False)

                # Combine urban and periurban polygons, join in conurbation labels
                logging.info(f'Combine urban and periurban polygons, join in conurbation labels')
                print(5)
                urban_periurban_refined = pd.concat([urban_buffer_refined, periurban_buffer_refined])
                urban_periurban_refined = gpd.sjoin(left_df = urban_periurban_refined, right_df = target_conurbations, how = 'left', predicate = 'intersects')
                urban_periurban_refined['conurbation_id'] = 'conurban_' + urban_periurban_refined['conurbation_id']
                urban_periurban_refined = urban_periurban_refined[['urban_id','area_type','conurbation_id','geometry']]

                assert urban_periurban_refined[urban_periurban_refined['urban_id'].duplicated()].shape[0] == 0
                assert urban_periurban_refined[urban_periurban_refined['conurbation_id'].isnull()].shape[0] == 0

                # Assign labels to blocks with relationship of block to urban or peri-urban boundaries
                logging.info(f'Assign labels to blocks with relationship of block to urban or peri-urban boundaries')
                logging.info(f'Assign conurbation labels')
                print(6)
                target_conurbations = target_conurbations.assign(within_conurbation='within_conurbation')
                target_conurbations = target_conurbations.assign(intersects_conurbation='intersects_conurbation')
                block_reference = blocks[['block_id','geometry']]
                block_reference = gpd.sjoin(left_df = block_reference, right_df = target_conurbations[['within_conurbation','geometry']], predicate = 'within', how = 'left').drop(['index_right'], axis=1)
                block_reference = gpd.sjoin(left_df = block_reference, right_df = target_conurbations[['intersects_conurbation','geometry']], predicate = 'intersects', how = 'left').drop(['index_right'], axis=1)

                target_conurbation_maxbuffer = target_conurbations[['geometry']].copy()
                target_conurbation_maxbuffer['geometry'] = target_conurbation_maxbuffer['geometry'].to_crs(3395).buffer(1).to_crs(4326)
                target_conurbation_maxbuffer = gpd.overlay(df1 = target_conurbation_maxbuffer, df2 = target_conurbations[['geometry']], how = 'difference')
                target_conurbation_maxbuffer = target_conurbation_maxbuffer.assign(intersects_conurbation_max='intersects_conurbation_max')
                block_reference = gpd.sjoin(left_df = block_reference, right_df = target_conurbation_maxbuffer[['intersects_conurbation_max','geometry']], predicate = 'intersects', how = 'left').drop(['index_right'], axis=1)

                logging.info(f'Assign urban labels')
                print(7)
                target_urban = urban_periurban_refined[urban_periurban_refined['area_type'] == 'Urban'][['geometry']]
                target_urban = target_urban.explode(index_parts=False)
                target_urban = target_urban.assign(within_urban='within_urban')
                block_reference = gpd.sjoin(left_df = block_reference, right_df = target_urban[['within_urban','geometry']], predicate = 'within', how = 'left').drop(['index_right'], axis=1)

                target_urban_maxbuffer = target_urban[['geometry']].copy()
                target_urban_maxbuffer['geometry'] = target_urban_maxbuffer['geometry'].to_crs(3395).buffer(10).to_crs(4326)
                target_urban['geometry'] = target_urban['geometry'].to_crs(3395).buffer(-10).to_crs(4326)
                target_urban_maxbuffer = gpd.overlay(df1 = target_urban_maxbuffer, df2 = target_urban[['geometry']], how = 'difference')
                target_urban_maxbuffer = target_urban_maxbuffer.assign(intersects_urban_max='intersects_urban_max')
                block_reference = gpd.sjoin(left_df = block_reference, right_df = target_urban_maxbuffer[['intersects_urban_max','geometry']], predicate = 'intersects', how = 'left').drop(['index_right'], axis=1)

                block_reference = block_reference[['block_id', 'within_conurbation', 'intersects_conurbation', 'intersects_conurbation_max', 'within_urban', 'intersects_urban_max']] #'touches_urban',
                block_reference = block_reference.drop_duplicates()

                block_reference['spatial_predicate'] = block_reference['within_urban'].combine_first(
                    block_reference['intersects_urban_max']).combine_first(
                    block_reference['within_conurbation']).combine_first(
                    block_reference['intersects_conurbation']).fillna('non_urban')
                logging.info(block_reference.groupby(['spatial_predicate'], as_index=False).count())
                assert block_reference[block_reference['block_id'].duplicated()].shape[0] == 0

                # Overlay urban-periurban polygons and assign blocks to either category based on share of area
                logging.info(f'Overlay urban-periurban polygons and assign blocks to either category based on share of area')
                print(8)
                overlay_list = block_reference[~block_reference['spatial_predicate'].isin(['non_urban'])]['block_id'].unique().tolist()
                blocks_conurban = blocks[blocks['block_id'].isin(overlay_list)]
                urban_periurban_blocks = gpd.overlay(df1 = blocks_conurban, df2 = urban_periurban_refined, how='intersection', keep_geom_type=True, make_valid=True)
                urban_periurban_blocks['area'] = urban_periurban_blocks['geometry'].to_crs(3395).area
                urban_periurban_blocks['rank'] = urban_periurban_blocks.groupby('block_id')['area'].rank(method='first', ascending=False)
                urban_periurban_blocks['share'] = urban_periurban_blocks['area']/urban_periurban_blocks['block_area']
                urban_periurban_blocks = urban_periurban_blocks.merge(block_reference, how = 'left', left_on=['block_id'], right_on = ['block_id'])

                # Filter to only blocks within courbation and limit to area type that it has most area in
                urban_periurban_blocks = urban_periurban_blocks[(urban_periurban_blocks['spatial_predicate'].isin(['within_conurbation','intersects_conurbation','within_urban','intersects_urban_max'])) & (urban_periurban_blocks['rank'] == 1)]
                # Filter out blocks that are more than 10km away from urban area and not touching urban area
                urban_periurban_blocks = urban_periurban_blocks[~((urban_periurban_blocks['intersects_urban_max'] != 'intersects_urban_max') & (urban_periurban_blocks['intersects_conurbation_max'] == 'intersects_conurbation_max'))]
                # Filter out if conurbation_id is null
                urban_periurban_blocks = urban_periurban_blocks[~urban_periurban_blocks['conurbation_id'].isna()]
                urban_periurban_blocks = urban_periurban_blocks[['block_id','urban_id','area_type','conurbation_id']]

                # Merge urban and periurban labels to blocks
                logging.info(f'Merge urban and periurban labels to blocks')
                print(9)
                urban_periurban_blocks = blocks_conurban.merge(urban_periurban_blocks, how = 'inner', left_on=['block_id'], right_on = ['block_id'])
                assert urban_periurban_blocks[urban_periurban_blocks['block_id'].duplicated()].shape[0] == 0

                # Dissolve blocks into contiguous conurbation area and explode to polygons
                logging.info(f'Dissolve blocks to contiguous conurbation area')
                print(10)
                urban_periurban_contiguous = urban_periurban_blocks[['conurbation_id','geometry']]
                urban_periurban_contiguous = dask_geopandas.from_geopandas(urban_periurban_contiguous, npartitions = 50)
                urban_periurban_contiguous = urban_periurban_contiguous.dissolve(by = 'conurbation_id').compute()
                urban_periurban_contiguous.reset_index(inplace = True)
                urban_periurban_contiguous['geometry'] = urban_periurban_contiguous['geometry'].make_valid()
                urban_periurban_contiguous = urban_periurban_contiguous.explode(index_parts=False)
                urban_periurban_contiguous = urban_periurban_contiguous[urban_periurban_contiguous['geometry'].geom_type == 'Polygon']

                # Create an urban flag and 500 meters buffer around urban area flag, join to contiguous exploded polygons
                logging.info(f'Label dissolved contiguous polygons as urban or peri-urban')
                print(11)
                urban_flag = urban_centers_refined.copy()
                urban_flag = urban_flag.assign(intersects_urban=1)
                urban_periurban_contiguous = gpd.sjoin(left_df = urban_periurban_contiguous, right_df = urban_flag[['intersects_urban','geometry']], predicate = 'intersects', how = 'left').drop(['index_right'], axis=1)

                periurban_minimum_flag = urban_centers_refined.copy()
                periurban_minimum_flag['geometry'] = periurban_minimum_flag['geometry'].to_crs(3395).buffer(500).to_crs(4326)
                periurban_minimum_flag = periurban_minimum_flag.assign(intersects_periurban=1)
                urban_periurban_contiguous = gpd.sjoin(left_df = urban_periurban_contiguous, right_df = periurban_minimum_flag[['intersects_periurban','geometry']], predicate = 'intersects', how = 'left').drop(['index_right'], axis=1)

                urban_periurban_contiguous = urban_periurban_contiguous[(urban_periurban_contiguous['intersects_urban'] == 1) | (urban_periurban_contiguous['intersects_periurban'] == 1)][['conurbation_id','geometry']]
                urban_periurban_contiguous = urban_periurban_contiguous.assign(contiguous_conurbation=1)
                urban_periurban_contiguous = urban_periurban_contiguous[['contiguous_conurbation','geometry']]

                urban_periurban_blocks = gpd.overlay(df1 = urban_periurban_blocks, df2 = urban_periurban_contiguous, how='intersection', keep_geom_type=True, make_valid=True)
                urban_periurban_blocks['area'] = urban_periurban_blocks['geometry'].to_crs(3395).area
                urban_periurban_blocks['rank'] = urban_periurban_blocks.groupby('block_id')['area'].rank(method='first', ascending=False)
                urban_periurban_blocks = urban_periurban_blocks[(urban_periurban_blocks['contiguous_conurbation'] == 1) & (urban_periurban_blocks['rank'] == 1)]

                # Merge urban and periurban polygon labels with blocks
                logging.info(f'Label blocks with contiguous urban or peri-urban areas')
                print(12)
                urban_periurban_blocks = urban_periurban_blocks.merge(block_reference, how = 'left', left_on=['block_id'], right_on = ['block_id'])
                urban_periurban_blocks = urban_periurban_blocks[['block_id','urban_id','area_type','conurbation_id']]

            else:
                urban_periurban_blocks = pd.DataFrame(columns = ['block_id','urban_id','area_type','conurbation_id'])

            blocks_xwalk = blocks.merge(urban_periurban_blocks, how = 'left', left_on=['block_id'], right_on = ['block_id'])

            blocks_xwalk['urban_id'] = blocks_xwalk['urban_id'].fillna('nonurban_' + blocks_xwalk['country_code'])
            blocks_xwalk['conurbation_id'] = blocks_xwalk['conurbation_id'].fillna('nonurban_' + blocks_xwalk['country_code'])
            blocks_xwalk['area_type'] = blocks_xwalk['area_type'].fillna('Non-urban')
            assert blocks_xwalk[blocks_xwalk['block_id'].duplicated()].shape[0] == 0
            assert blocks_xwalk[blocks_xwalk['block_id'].isnull()].shape[0] == 0
            # blocks_xwalk[blocks_xwalk['area_type'] != 'Non-urban'].explore('urban_id')

            # Label with Africapolis
            logging.info(f'Blocks and Africapolis overlay intersect')
            print(13)
            blocks_africapolis_within = gpd.sjoin(left_df = blocks[['block_id','geometry']], right_df = africapolis_data, predicate = 'within', how = 'inner').drop(['index_right'], axis=1)
            blocks_africapolis_within = blocks_africapolis_within[['block_id', 'agglosid', 'agglosname', 'metropole']]
            blocks_africapolis_intersects = gpd.sjoin(left_df = blocks[['block_id','geometry']], right_df = africapolis_data, predicate = 'intersects', how = 'inner').drop(['index_right'], axis=1)
            blocks_africapolis_intersects = blocks_africapolis_intersects[~blocks_africapolis_intersects['block_id'].isin(blocks_africapolis_within['block_id'].unique())]
            blocks_africapolis_intersects = blocks[blocks['block_id'].isin(blocks_africapolis_intersects['block_id'].unique())]
            if blocks_africapolis_intersects.shape[0] > 0:
                blocks_africapolis = gpd.overlay(df1 = blocks_africapolis_intersects, df2 = africapolis_data, how='intersection', keep_geom_type=True, make_valid=True)
                blocks_africapolis['area'] = blocks_africapolis['geometry'].to_crs(3395).area
                blocks_africapolis['rank'] = blocks_africapolis.groupby('block_id')['area'].rank(method='first', ascending=False)
                blocks_africapolis = blocks_africapolis[blocks_africapolis['rank'] == 1]
                blocks_africapolis = blocks_africapolis[['block_id', 'agglosid', 'agglosname', 'metropole']]
                blocks_africapolis = pd.concat([blocks_africapolis, blocks_africapolis_within]).reset_index(drop = True)
            else:
                blocks_africapolis = blocks_africapolis_within

            blocks_xwalk = blocks_xwalk.merge(blocks_africapolis, how = 'left', left_on=['block_id'], right_on = ['block_id'])
            blocks_xwalk['agglosid'] = blocks_xwalk['agglosid'].fillna('no_agglosid')
            blocks_xwalk['agglosname'] = blocks_xwalk['agglosname'].fillna('Non-urban')
            blocks_xwalk['metropole'] = blocks_xwalk['metropole'].fillna('No')
            # blocks_xwalk[blocks_xwalk['agglosname'] != 'Non-urban'].explore('agglosname')
            # blocks_xwalk[blocks_xwalk['area_type'] != 'Non-urban'].explore('urban_id')

            blocks_xwalk.to_file(filename = Path(dask_dir) / f'block_xwalk_{country_code}.gpkg')
            blocks_xwalk = blocks_xwalk[['block_id','block_geohash','gadm_code','country_code','urban_id','conurbation_id','area_type','agglosid', 'agglosname', 'metropole']]

            logging.info(f'Writing to dask partition')
            blocks_xwalk = dask.dataframe.from_pandas(data = blocks_xwalk, npartitions = 1) 
            dask.dataframe.to_parquet(df = blocks_xwalk, path = Path(dask_dir) / f'crosswalks.parquet', engine='pyarrow', compression='snappy', append=True, ignore_divisions=True)

        # Read parquet data
        print('Preparing crosswalk')
        full_xwalk = dask.dataframe.read_parquet(path = Path(dask_dir)).compute()

        # Join GHSL data to crosswalk
        full_xwalk = full_xwalk.merge(urban_centers[['urban_id','urban_center_name','urban_area']], how = 'left', left_on=['urban_id'], right_on = ['urban_id'])

        # Coerce all conurbations to have at least one urban_id (applies to blocks on fringes of )
        full_xwalk['urban_area'] = full_xwalk.groupby(['urban_id'])['urban_area'].transform('sum')
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

        try: full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name'].str.split('-', n = 2, expand = True)[[0,1]].dropna().astype(str).apply('-'.join, 1)
        except: full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name'].fillna(full_xwalk['conurbation_area_name'])
        full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name_short'].fillna(full_xwalk['conurbation_area_name'])

        # Fill in for primary country fields
        full_xwalk['conurbation_country_name'] = full_xwalk['conurbation_country_name'].fillna(full_xwalk['country_name'])
        full_xwalk['urban_country_name'] = full_xwalk['urban_country_name'].fillna(full_xwalk['country_name'])

        # Join in conurbation ID ranking for labels
        urban_area_rank = full_xwalk[full_xwalk['urban_id'].notnull()].groupby(['urban_id','conurbation_id']).agg({'urban_area': 'sum'}).reset_index()
        urban_area_rank['rank'] = urban_area_rank.groupby('conurbation_id')['urban_area'].rank(method='first', ascending=False)
        full_xwalk = full_xwalk.merge(urban_area_rank[['urban_id','conurbation_id','rank']], how = 'left', left_on = ['urban_id','conurbation_id'], right_on = ['urban_id','conurbation_id'])

        print('Fully merged')
        # 'Core urban', 'Peripheral urban', 'Peri-urban', 'Non-urban place', 'Non-urban'
        conditions = [(full_xwalk['area_type'] == 'Urban') & (full_xwalk['rank'] == 1),
                        (full_xwalk['area_type'] == 'Urban') & (full_xwalk['rank'] > 1),
                        (full_xwalk['area_type'] == 'Peri-urban'),
                        (full_xwalk['area_type'] == 'Non-urban') & (full_xwalk['agglosid'] != 'no_agglosid'),
                        (full_xwalk['area_type'] == 'Non-urban')]
        labels = ['1 - Core urban', '2 - Peripheral urban', '3 - Peri-urban', '4 - Non-urban place', '5 - Non-urban']
        full_xwalk["class_urban_hierarchy_detail"] = np.select(conditions, labels, default='5 - Non-urban')
        # full_xwalk.groupby(['class_urban_hierarchy_detail']).size()
        # full_xwalk[full_xwalk['class_urban_hierarchy_detail'] != '5 - Non-urban'].explore('class_urban_hierarchy_detail')

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

        # Fill in non-urban center name labels
        conditions = [(full_xwalk['urban_center_name'].isnull()) & (full_xwalk["class_urban_hierarchy"] ==  '3 - Peri-urban'),
                        (full_xwalk['urban_center_name'].isnull()) & (full_xwalk['area_type'] == 'Non-urban')]
        labels = ['Peri-urban ' + full_xwalk['conurbation_area_name_short'], 'Rest of ' + full_xwalk['country_name']]
        full_xwalk['urban_center_name'] = np.select(conditions, labels, default = full_xwalk['urban_center_name'])
        full_xwalk['urban_center_name'] = full_xwalk['urban_center_name'].fillna('Rest of ' + full_xwalk['country_name']) # not necessary really

        full_xwalk['urban_country_name'] = full_xwalk['urban_country_name'].fillna(full_xwalk['country_name'])
        full_xwalk['urban_country_code'] = full_xwalk['urban_country_code'].fillna(full_xwalk['country_code'])
        full_xwalk['conurbation_area_name'] = full_xwalk['conurbation_area_name'].fillna('Rest of ' + full_xwalk['country_name'])
        full_xwalk['conurbation_country_name'] = full_xwalk['conurbation_country_name'].fillna(full_xwalk['country_name'])
        full_xwalk['conurbation_country_code'] = full_xwalk['conurbation_country_code'].fillna(full_xwalk['country_code'])
        full_xwalk['conurbation_area_name_short'] = full_xwalk['conurbation_area_name_short'].fillna('Rest of ' + full_xwalk['country_name'])

        full_xwalk = full_xwalk[['block_id', 'block_geohash', 'gadm_code', 'country_code', 'country_name', 'continent', 'area_type', 'class_urban_hierarchy_detail', 'class_urban_hierarchy', 'class_urban_periurban_nonurban', 'class_urban_nonurban', 'urban_id', 'urban_center_name', 'urban_country_code', 'urban_country_name', 'conurbation_id', 'conurbation_area_name', 'conurbation_area_name_short', 'conurbation_country_code', 'conurbation_country_name', 'agglosid', 'agglosname', 'metropole']]

        # Urban layer code
        full_xwalk.loc[full_xwalk['area_type'].isin(['Non-urban']), "urban_layer_code"] = 'nonurban' + '_' + full_xwalk['country_code']
        full_xwalk.loc[~full_xwalk['area_type'].isin(['Non-urban']), "urban_layer_code"] = full_xwalk['country_code'] + '_' + full_xwalk['conurbation_id'] + '_' +  full_xwalk['urban_id']

        assert full_xwalk[full_xwalk['block_id'].duplicated()].shape[0] == 0
        assert full_xwalk[full_xwalk['block_id'].isnull()].shape[0] == 0

        # Write file to parquet CSV
        print('Writing files')
        logging.info(f'Writing crosswalk.')
        full_xwalk.to_parquet(path = Path(crosswalk_dir) / f'ghsl_crosswalk.parquet')
        full_xwalk.to_csv(Path(crosswalk_dir) / f'ghsl_crosswalk.csv')

    # Aggregate boundaries
    if len(country_list_regions) > 0:

        full_xwalk = pd.read_parquet(path = Path(crosswalk_dir) / f'ghsl_crosswalk.parquet')
        logging.info(f'Create aggregate boundaries')
        urban_delineations_labels = full_xwalk[['urban_layer_code', 'country_code', 'country_name', 'area_type', 'class_urban_hierarchy', 'class_urban_periurban_nonurban', 'class_urban_nonurban', 'urban_id', 'urban_center_name', 'urban_country_code', 'urban_country_name', 'conurbation_id', 'conurbation_area_name', 'conurbation_area_name_short', 'conurbation_country_code', 'conurbation_country_name']]
        urban_delineations_labels = urban_delineations_labels.drop_duplicates()
        assert urban_delineations_labels[urban_delineations_labels['urban_layer_code'].duplicated()].shape[0] == 0

        urban_delineations = gpd.GeoDataFrame({'urban_layer_code': pd.Series(dtype='str'), 'geometry': gpd.GeoSeries(dtype='geometry')}).set_crs(epsg=4326)
        for country_code in country_list_regions:
            logging.info(f'{country_code}')
            blocks = gpd.read_parquet(path = Path(blocks_dir) / f'blocks_{country_code}.parquet', memory_map = True)
            regions = gpd.read_parquet(path = Path(gadm_dir) / 'parquet' / f'gadm_{country_code}.parquet', memory_map = True)
            regions = regions[['country_code','geometry']]

            urban_boundary = pd.merge(blocks, full_xwalk[['block_id','urban_layer_code','area_type']], on = 'block_id', how = 'left')
            urban_boundary = urban_boundary[urban_boundary['area_type'] != 'Non-urban']
            urban_boundary = urban_boundary[['urban_layer_code','geometry']]
            urban_boundary['geometry'] = urban_boundary['geometry'].to_crs(3395).buffer(0.0001).to_crs(4326)
            logging.info(f'Dissolving urban boundaries')
            urban_boundary = dask_geopandas.from_geopandas(urban_boundary, npartitions = 50)
            urban_boundary = urban_boundary.dissolve(by = 'urban_layer_code').compute()
            urban_boundary.reset_index(inplace = True)

            regions = gpd.overlay(df1 = regions, df2 = urban_boundary[['geometry']], how = 'difference')
            logging.info(f'Dissolving regions')
            regions = dask_geopandas.from_geopandas(regions, npartitions = 50)
            regions = regions.dissolve(by = 'country_code').compute()
            regions.reset_index(inplace = True)
            regions['urban_layer_code'] = 'nonurban' + '_' + regions['country_code']
            regions = regions[['urban_layer_code','geometry']]
            urban_boundary = pd.concat([urban_boundary, regions], ignore_index=True)
            urban_boundary = urban_boundary.merge(urban_delineations_labels, on = 'urban_layer_code', how = 'left')
            urban_delineations = pd.concat([urban_delineations, urban_boundary], ignore_index=True)

        logging.info(f'Writing boundaries')
        urban_delineations.to_parquet(path = Path(crosswalk_dir) / f'urban_boundaries.parquet')
        urban_delineations.to_file(filename = Path(crosswalk_dir) / f'urban_boundaries.gpkg')

    logging.info(f'Finished.')
    print('Finished')


def setup(args=None):
    parser = argparse.ArgumentParser(description='Compute crosswalk.')
    parser.add_argument('--log_file', required=False, type=Path, dest="log_file", help="Path to write log file.") 
    parser.add_argument('--country_chunk', required=False, type=str, dest="country_chunk", nargs='+', help="List of country codes following ISO 3166-1 alpha-3 format.")
    parser.add_argument('--africapolis_file', required=True, type=Path, dest="africapolis_file", help="Path to Africapolis file.")
    parser.add_argument('--ghsl_file', required=True, type=Path, dest="ghsl_file", help="Path to GHSL file.")
    parser.add_argument('--gadm_dir', required=True, type=Path, dest="gadm_dir", help="Path to GADM directory.")
    parser.add_argument('--blocks_dir', required=True, type=Path, dest="blocks_dir", help="Path to blocks directory.")
    parser.add_argument('--output_dir', required=True, type=Path, dest="output_dir", help="Path to outputs directory.")
    return parser.parse_args(args)

if __name__ == "__main__":
    main(**vars(setup()))
