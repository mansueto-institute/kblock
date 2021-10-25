import numpy as np 
import geopandas as gpd 
import pandas as pd 
from shapely.wkt import loads
from typing import Tuple, Union, List
from pathlib import Path 
import argparse
from pygeos import GEOSException
import time
import os
import pyarrow

#from . import utils, block_stats
import utils
import block_stats
from raster_tools import extract_aoi_data_from_raster, allocate_population


def load_gadm_file(gadm_dir: str) -> gpd.GeoDataFrame:
    """
    Loads in a GADM gep df from file, sorts the index and returns
    just the index and geometry. 
    """
    sort_fn = lambda s: -int(s.stem.split("_")[-1])
    gadm_dir = Path(gadm_dir)
    shp_files = [p for p in gadm_dir.iterdir() if p.suffix == '.shp']
    shp_files.sort(key=sort_fn)
    gdf = gpd.read_file(str(shp_files[0]))
    n = abs(sort_fn(shp_files[0]))
    s = 'GID_{}'.format(n)
    gdf.rename(columns={s: 'gadm'}, inplace=True)
    gdf = gdf[['gadm', 'geometry']]
    return gdf 

# def make_summary_in_memory(superblock_blocks: gpd.GeoDataFrame, 
#                            landscan_path: Union[str, Path],
#                            superblock_buildings: gpd.GeoDataFrame,
#                            summary_out_path: Union[str, Path],
#                            ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
#     """
#     make_summary_in_memory():

#     Takes in all of the superblock data and the superblocks data a
#     """
#     gadm_list = list(set(superblock_blocks['gadm_code']))
#     _, superblock_ls = extract_aoi_data_from_raster(superblock_blocks, landscan_path, save_geojson=False, save_tif=False)

#     bldg_pop_alloc = allocate_population(superblock_buildings, superblock_ls, 'pop')

#     superblock_bldg_summary = block_stats.make_superblock_summary(bldg_pop_alloc, superblock_blocks)
#     block_cols = [x for x in superblock_bldg_summary.columns if "block" in x]
#     superblock_stats = superblock_bldg_summary[block_cols].drop_duplicates()
#     superblock_summary = superblock_blocks.merge(superblock_stats, how='left', on='block_id')

#     bldg_pop_alloc = allocate_population(superblock_buildings, superblock_ls, 'pop')

#     # (2) Now assemble the other data
#     superblock_bldg_summary = block_stats.make_superblock_summary(bldg_pop_alloc, superblock_blocks)
#     block_cols = [x for x in superblock_bldg_summary.columns if "block" in x]
#     superblock_stats = superblock_bldg_summary[block_cols].drop_duplicates()
#     superblock_summary = superblock_blocks.merge(superblock_stats, how='left', on='block_id')

#     # (3) Save
#     summary_out_path = Path(summary_out_path)
#     fname = summary_out_path.stem
#     outdir = summary_out_path.parent
#     outdir.mkdir(exist_ok=True, parents=True)

#     superblock_summary = utils.remove_duplicated_cols_from_merge(superblock_summary)
#     superblock_buildings_out_path = outdir / (fname + "-bldgs" + summary_out_path.suffix)
#     if summary_out_path.suffix == '.geojson':
#         superblock_summary.to_file(str(summary_out_path), driver='GeoJSON')
#         superblock_bldg_summary.to_file(str(superblock_buildings_out_path), driver='GeoJSON')
#     elif summary_out_path.name.ends_with('parquet'):
#         superblock_summary.to_parquet(str(summary_out_path))
#         superblock_buildings_out_path.to_parquet(str(superblock_buildings_out_path))
#     print("Saved to: {}".format(str(summary_out_path)))
#     print("Saved to: {}".format(str(superblock_buildings_out_path)))

#     return superblock_summary, superblock_bldg_summary



def make_summary(superblock_blocks: Union[str, Path, gpd.GeoDataFrame],
                 landscan_path: Union[str, Path],
                 superblock_buildings: Union[str, Path, gpd.GeoDataFrame],
                 summary_out_path: Union[str, Path],
                 ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Creates a summary of the given area using its population, blocks,
    and boundary lines. 
    """
    # (1) Read in the superblock data if provided a path, or pass over this
    #     logic if the superblock data is already in memory
    if not isinstance(superblock_blocks, gpd.GeoDataFrame):
        if isinstance(superblock_blocks, str):
            superblock_path = Path(superblock_blocks)
        elif isinstance(superblock_blocks, Path):
            superblock_path = superblock_blocks
        else: 
            raise Exception(f'Unknown value of type {type(superblock_blocks)} passed as superblock_blocks')
        superblock_blocks = gpd.read_file(superblock_path)

    if not isinstance(superblock_buildings, gpd.GeoDataFrame):
        if isinstance(superblock_buildings, str):
            superblock_buildings_path = Path(superblock_buildings)
        elif isinstance(superblock_buldings, Path):
            superblock_buildings_path = superblock_buildings
        else: 
            raise Exception(f'Unknown value of type {type(superblock_buildings)} passed as superblock_buildings')
        superblock_buildings = gpd.read_file(superblock_buildings_path)


    # (2) Allocate Landscan
    gadm_list = list(set(superblock_blocks['gadm_code']))
    _, superblock_ls = extract_aoi_data_from_raster(superblock_blocks, landscan_path, save_geojson=False, save_tif=False)

    ### fiona error ###
    if superblock_buildings is None:
        summary_out_path = Path(summary_out_path)
        fname = summary_out_path.stem
        outdir = summary_out_path.parent
        outdir.mkdir(exist_ok=True, parents=True)

        # create empty files
        with open(fname + ".geojson", 'w') as fp:
            pass
        with open(fname + "-bldgs.geojson", 'w') as fp:
            pass

        print("No building file: {}".format(fname))

        return None
    ### --- ###
    bldg_pop_alloc = allocate_population(superblock_buildings, superblock_ls, 'pop')

    # (2) Now assemble the other data
    superblock_bldg_summary = block_stats.make_superblock_summary(bldg_pop_alloc, superblock_blocks)
    block_cols = [x for x in superblock_bldg_summary.columns if "block" in x]
    superblock_stats = superblock_bldg_summary[block_cols].drop_duplicates()
    superblock_summary = superblock_blocks.merge(superblock_stats, how='left', on='block_id')

    # (3) Save
    summary_out_path = Path(summary_out_path)
    fname = summary_out_path.stem
    outdir = summary_out_path.parent
    outdir.mkdir(exist_ok=True, parents=True)

    superblock_summary = utils.remove_duplicated_cols_from_merge(superblock_summary)
    superblock_buildings_out_path = outdir / (fname + "-bldgs" + summary_out_path.suffix)
    if summary_out_path.suffix == '.geojson':
        superblock_summary.to_file(summary_out_path, driver='GeoJSON')
        superblock_bldg_summary.to_file(superblock_buildings_out_path, driver='GeoJSON')
    elif summary_out_path.name.ends_with('parquet'):
        utils.parquet_write(superblock_blocks, summary_out_path)
        utils.parquet_write(superblock_buildings, superblock_buildings_out_path)
    else:
        raise Exception(f'Out path file format {superblock_buildings_out_path.suffix} not recognized')
    print("Saved to: {}".format(str(summary_out_path)))
    print("Saved to: {}".format(str(superblock_buildings_out_path)))

    return superblock_summary, superblock_bldg_summary


if __name__ == "__main__":
    t0 = time.time()
    parser = argparse.ArgumentParser(description='Make block-level and building-level summary for Area of Interest')
    parser.add_argument('--superblock_path', required=True, type=str, help='Path to geometry which defines AoI')
    parser.add_argument('--landscan_path', required=True, type=str, help='Path to Landscan tif file')
    parser.add_argument('--buildings_path', required=True, type=str, help='Path to buildings geomatries')
    parser.add_argument('--summary_out_path', required=True, type=str, help='Path to save block summary')
    args = parser.parse_args()
    make_summary(**vars(args))
    t1 = time.time()
    print(f"block summary took {t1-t0}")
