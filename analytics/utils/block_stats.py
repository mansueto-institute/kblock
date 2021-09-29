import numpy as np 
import geopandas as gpd 
import pandas as pd 
from shapely.wkt import loads
from shapely.geometry import Polygon, MultiPolygon
from typing import Tuple, Union, Optional
from pathlib import Path 
from momepy import Rectangularity, BuildingAdjacency, Gini
import utils
import libpysal
from libpysal.weights import W

'''
FILE DESCRIPTION:
Provides capacity to generate block-level metrics.
Structure is to:
    1. Load a building geomtry w/ bldg level pop allocation
    2. If needed, add the block_id
    3. Then there are functions to add additional columns 
       to the GeoDataFrame including
        - block_area
        - block_bldg_count
        - block_bldg_density
        - block_pop_total
        - block_pop_density
    4. Save this out, either maintaining the bldg-level detail
       or reducing to block-level
'''

def flex_load(block: Union[gpd.GeoDataFrame, str]) -> gpd.GeoDataFrame:
    """
    flex_load
    Helper function to allow downstream fns to accept 
    either a path to a GeoDataFrame or the dataframe itself,
    depending on whether you've already loaded it or not. Checks if 
    the object named block being passed in is a string or Path, and if so
    reads it in and returns it.
    """
    if isinstance(block, str) or isinstance(block, Path):
        block = utils.load_csv_to_geo(block)
    return block     


def load_bldg_pop(bldg_pop_path: str,
                  pop_variable: str = 'bldg_pop',
                  ) -> gpd.GeoDataFrame:
    """
    load_bldg_pop
    Loads in a building's geometry with building level population info
    """
    bldg_pop = gpd.read_file(bldg_pop_path)
    assert pop_variable in bldg_pop.columns, "ERROR - loading the building level pop file but looking for pop column |{}| which is not in file {}".format(pop_variable, bldg_pop_path)
    return bldg_pop


def add_block_id(bldg_pop: gpd.GeoDataFrame,
                 block: Union[gpd.GeoDataFrame, str],
                 ) -> gpd.GeoDataFrame:
    """
    add_block_id()
    Step 2: some bldg files don't have the block_id so that may need 
    to be joined on
    NOTE: block can be a path to the block GeoDataFrame, or the already loaded GeoDataFrame
    Joins block_id column on to the builing geodf.
    """
    block = flex_load(block)
    bldg_pop = utils.join_block_building(block, bldg_pop)
    if 'index_right' in bldg_pop.columns:
        bldg_pop.drop(columns=['index_right'], inplace=True)
    return bldg_pop


#######################################
# BASIC BLOCK-LEVEL STATISTICS TO ADD #
#######################################

def add_block_area(bldg_pop: gpd.GeoDataFrame,
                   block: gpd.GeoDataFrame,
                   ) -> gpd.GeoDataFrame:
    """
    Calculates the area of each block and adds that to the bldg_pop geodf
    """    
    block = block.to_crs("EPSG:3395")
    block['block_area'] = block.area * 1e-6
    block = block.to_crs("EPSG:4326")

    bldg_pop = bldg_pop.merge(block[['block_id', 'block_area']],
                              how='left', on='block_id')
    return bldg_pop


def add_block_bldg_count(bldg_pop: gpd.GeoDataFrame
                         ) -> gpd.GeoDataFrame:
    """
    Calculates the number of buildings in a block and adds that to the bldg_pop geodf
    """
    counts = bldg_pop[['block_id', 'bldg_id']].groupby('block_id').count().reset_index()
    counts.rename(columns={'bldg_id': 'block_bldg_count'}, inplace=True)
    bldg_pop = bldg_pop.merge(counts, how='left', on='block_id')
    return bldg_pop


def add_block_bldg_area(bldg_pop: gpd.GeoDataFrame,
                        block: gpd.GeoDataFrame,
                        ) -> gpd.GeoDataFrame:
    """
    Calculates the number of buildings in a block and adds that to the bldg_pop geodf
    """
    bldg_pop = bldg_pop.to_crs("EPSG:3395")
    bldg_pop['bldg_area'] = (bldg_pop.area * 1e-6)
    block_bldg_area = bldg_pop[['block_id', 'bldg_area']].groupby('block_id').sum().reset_index()
    block_bldg_area.rename(columns={'bldg_area': 'block_bldg_area'}, inplace=True)

    bldg_pop = bldg_pop.merge(block_bldg_area, how='left', on='block_id')
    bldg_pop = bldg_pop.to_crs("EPSG:4326")
    bldg_pop.drop(columns=["bldg_area"], inplace=True)
    return bldg_pop


def add_block_bldg_area_density(bldg_pop: gpd.GeoDataFrame,
                                block: gpd.GeoDataFrame,
                                ) -> gpd.GeoDataFrame:
    """
    Calculates the ratio of building density to block area and adds that to the bldg_pop geodf
    """
    if 'block_bldg_area' not in bldg_pop.columns:
        bldg_pop = add_block_bldg_area(bldg_pop, block)

    if 'block_area' not in bldg_pop.columns:
        bldg_pop = add_block_area(bldg_pop, block)

    bldg_pop['block_bldg_area_density'] = bldg_pop['block_bldg_area'] / bldg_pop['block_area']
    return bldg_pop


def add_block_bldg_count_density(bldg_pop: gpd.GeoDataFrame,
                                 block: gpd.GeoDataFrame,
                                 ) -> gpd.GeoDataFrame:
    """
    Calculates the ratio of number of buildings in a block to the block's area and adds that
    to the bldg_pop geodf
    """
    if 'block_bldg_count' not in bldg_pop.columns:
        bldg_pop = add_block_bldg_count(bldg_pop, block)

    if 'block_area' not in bldg_pop.columns:
        bldg_pop = add_block_area(bldg_pop, block) 

    bldg_pop['block_bldg_count_density'] = bldg_pop['block_bldg_count'] / bldg_pop['block_area']
    return bldg_pop


def add_block_pop(bldg_pop: gpd.GeoDataFrame,
                  ) -> gpd.GeoDataFrame:
    """
    Calculates the population for the block and adds that to the bldg_pop geodf
    """
    block_pop = bldg_pop[['block_id', 'bldg_pop']].groupby('block_id').sum()
    block_pop.rename(columns={'bldg_pop': 'block_pop'}, inplace=True)
    bldg_pop = bldg_pop.merge(block_pop, how='left', on='block_id')
    return bldg_pop


def add_block_pop_density(bldg_pop: gpd.GeoDataFrame,
                          block: gpd.GeoDataFrame,
                          ) -> gpd.GeoDataFrame:
    """
    Calculates the ratio of block population to block area and adds that to the bldg_pop geodf
    """    
    if 'block_area' not in bldg_pop.columns:
        bldg_pop = add_block_area(bldg_pop, block)

    if 'block_pop' not in bldg_pop.columns:
        bldg_pop = add_block_pop(bldg_pop, block)

    bldg_pop['block_pop_density'] = bldg_pop['block_pop'] / bldg_pop['block_area']
    return bldg_pop  


def add_block_gini(bldg_pop: gpd.GeoDataFrame,
                   values_col: str,
                   weights: Optional[libpysal.weights.weights.W] = None
                   ) -> gpd.GeoDataFrame:
    """
    Calculates the gini coefficient for the passed in values_col based on the passed in gpd.GeoDataFrame.
    """
    if values_col not in bldg_pop.columns:
        raise IndexError(f'{values_col} was not found in bldg_pop')

    if weights is None:
        weights_dict = {bldg_pop.loc[index_val]['bldg_id']: bldg_pop[bldg_pop['block_id'] == bldg_pop.loc[index_val]['block_id']]['bldg_id'].tolist() for index_val in bldg_pop.index.tolist()}
        weights = W(weights_dict, ids=bldg_pop['bldg_id'].tolist())

    bldg_pop['gini_'+values_col] = Gini(bldg_pop, values=values_col, spatial_weights=weights, unique_id='bldg_id').series

    return bldg_pop


def add_block_rectangularity(bldg_pop: gpd.GeoDataFrame,
                             block: gpd.GeoDataFrame,
                             ) -> gpd.GeoDataFrame:
    """
    Adds rectangularity for each block in the data frame and adds that to the bldg_pop geodf.
    Calculated as (area)/(minimum bounding rotated rectangle area)
    See: http://docs.momepy.org/en/stable/generated/momepy.Rectangularity.html#momepy.Rectangularity
    """
    if 'block_area' not in bldg_pop.columns:
        bldg_pop = add_block_area(bldg_pop, block)

    bldg_pop['rectangularity'] = Rectangularity(bldg_pop, areas='block_area').series

    return bldg_pop


def add_building_adjacency(bldg_pop: gpd.GeoDataFrame, 
                           weights: Optional[libpysal.weights.weights.W] = None
                           ) -> gpd.GeoDataFrame:
    """
    Calculates the building agency, roughly a ratio of number of buildings to number of "built-up patches".
    The "spatial_weights_higher" is created by defining every building in the same block to be connected to one another.

    The intention is to figure out whether the block is composed of adjoining buildings or whether it's primarily freestanding
    buildings. 
    Defined here: https://www-sciencedirect-com.proxy.uchicago.edu/science/article/pii/S0169204617301275?via%3Dihub
    Documentation here: http://docs.momepy.org/en/stable/generated/momepy.BuildingAdjacency.html#momepy.BuildingAdjacency
    """
    if weights is None:
        weights_dict = {bldg_pop.loc[index_val]['bldg_id']: bldg_pop[bldg_pop['block_id'] == bldg_pop.loc[index_val]['block_id']]['bldg_id'].tolist() for index_val in bldg_pop.index.tolist()}
        weights = W(weights_dict, ids=bldg_pop['bldg_id'].tolist())

    bldg_pop['building_adjacency'] = BuildingAdjacency(bldg_pop, spatial_weights_higher=weights, unique_id='bldg_id').series

    return bldg_pop


######################################
# COMMANDS FOR GENERAL AOI SUMMARIES #
######################################
def make_aoi_summary(bldg_pop_data: Union[str, gpd.GeoDataFrame], 
                     block_data: Union[str, gpd.GeoDataFrame],
                     aoi_out_path: str = None,
                     ) -> None:
    '''
    Calculates all statistics given:
        1. bldg-level pop allocation
        2. block geometry
        3. path to save output to
    '''

    if isinstance(bldg_pop_data, gpd.GeoDataFrame):
        bldg_pop = bldg_pop_data
    else:
        bldg_pop = load_bldg_pop(bldg_pop_data)
    block = flex_load(block_data)
    if 'block_id' not in bldg_pop.columns:
        bldg_pop = add_block_id(bldg_pop, block)
    # bldg_pop = add_block_gini(bldg_pop, 'bldg_pop')
    # bldg_pop = add_building_adjacency(bldg_pop)
    bldg_pop = add_block_area(bldg_pop, block)
    bldg_pop = add_block_rectangularity(bldg_pop, block)
    bldg_pop = add_block_bldg_count(bldg_pop)
    bldg_pop = add_block_bldg_area(bldg_pop, block)
    bldg_pop = add_block_bldg_area_density(bldg_pop, block)
    bldg_pop = add_block_bldg_count_density(bldg_pop, block)
    bldg_pop = add_block_pop(bldg_pop)
    bldg_pop = add_block_pop_density(bldg_pop, block)

    if aoi_out_path is not None:
        aoi_out_path = Path(aoi_out_path)
        aoi_out_path.parent.mkdir(parents=True, exist_ok=True)
        bldg_pop.to_file(str(aoi_out_path), driver='GeoJSON')
    return bldg_pop
