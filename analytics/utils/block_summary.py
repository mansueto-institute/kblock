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


def make_summary(superblock_path: Union[str, Path],
                 landscan_path: Union[str, Path],
                 buildings_path: Union[str, Path],
                 summary_out_path: Union[str, Path],
                 pop_only: bool
                 ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Creates a summary of the given area using its population, blocks,
    and boundary lines. 
    """
    # (1) Allocate Landscan
    if not isinstance(superblock_path, Path):
        superblock_path = Path(superblock_path)
    if superblock_path.suffix == '.csv':
        superblock_blocks = utils.load_csv_to_geo(superblock_path)
    else:
        superblock_blocks = gpd.read_file(str(superblock_path))

    gadm_list = list(set(superblock_blocks['gadm_code']))
    _, superblock_ls = extract_aoi_data_from_raster(superblock_blocks, landscan_path, save_geojson=False, save_tif=False)

    superblock_buildings = gpd.read_file(buildings_path)

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
    if not pop_only:
        superblock_bldg_summary = block_stats.make_superblock_summary(bldg_pop_alloc, superblock_blocks)
        block_cols = [x for x in superblock_bldg_summary.columns if "block" in x]
        superblock_stats = superblock_bldg_summary[block_cols].drop_duplicates()
        superblock_summary = superblock_blocks.merge(superblock_stats, how='left', on='block_id')

    # (3) Save
    summary_out_path = Path(summary_out_path)
    fname = summary_out_path.stem
    outdir = summary_out_path.parent
    outdir.mkdir(exist_ok=True, parents=True)

    superblock_summary.to_file(str(summary_out_path), driver='GeoJSON')
    print("Saved to: {}".format(str(summary_out_path)))
    
    superblock_buildings_out_path = outdir / (fname + "-bldgs.geojson")
    superblock_bldg_summary.to_file(str(superblock_buildings_out_path), driver='GeoJSON')
    print("Saved to: {}".format(str(superblock_buildings_out_path)))

    return superblock_summary, superblock_bldg_summary


if __name__ == "__main__":
    t0 = time.time()
    parser = argparse.ArgumentParser(description='Make block-level and building-level summary for Area of Interest')
    parser.add_argument('--superblock_path', required=True, type=str, help='Path to geometry which defines AoI')
    parser.add_argument('--landscan_path', required=True, type=str, help='Path to Landscan tif file')
    parser.add_argument('--buildings_path', required=True, type=str, help='Path to buildings geomatries')
    parser.add_argument('--summary_out_path', required=True, type=str, help='Path to save block summary')
    parser.add_argument('--pop_only', required=False, default=False, action='store_true', help='Should only the population be computed, not the analytics?')
    args = parser.parse_args()
    make_summary(**vars(args))
    t1 = time.time()
    print(f"block summary took {t1-t0}")
