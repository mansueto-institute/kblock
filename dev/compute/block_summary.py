import numpy as np
import geopandas as gpd 
import pandas as pd 
from shapely.wkt import loads
from shapely.geometry import (MultiPolygon, Point, Polygon, 
                              MultiLineString, LineString, MultiPoint)
from typing import Tuple, Union, Optional, List
from pathlib import Path 
from pygeos import GEOSException
import rasterio 
from shapely.ops import unary_union
from rasterio.crs import CRS 
from rasterio import features 
import affine 
from tobler.area_weighted import _area_tables, area_interpolate
from tobler.util.util import _check_crs, _nan_check, _check_presence_of_crs
import rasterio 
import rasterio.features


ShapelyGeom = Union[MultiPolygon, Polygon, MultiLineString, LineString]


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
        block = load_csv_to_geo(block)
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
    bldg_pop = join_block_building(block, bldg_pop)
    if 'index_right' in bldg_pop.columns:
        bldg_pop.drop(columns=['index_right'], inplace=True)
    return bldg_pop


#######################################
# BASIC BLOCK-LEVEL STATISTICS TO ADD #
#######################################

def set_dtypes(bldg_pop: gpd.GeoDataFrame,
               block: gpd.GeoDataFrame,
              ) -> gpd.GeoDataFrame:
    """
    Mandate dtypes for each col
    """    
    dtypes = {'block_id': str, 'block_area': float, 'building_count': float, 'building_area': float}

    for col in dtypes.keys():
        bldg_pop[col] = bldg_pop[col].apply(lambda x: dtypes[col](x))
    return bldg_pop


def add_block_bldg_area_density(bldg_pop: gpd.GeoDataFrame,
                                block: gpd.GeoDataFrame,
                                ) -> gpd.GeoDataFrame:
    """
    Calculates the ratio of building density to block area and adds that to the bldg_pop geodf
    """
    bldg_pop['block_bldg_area_density'] = bldg_pop['building_area'] / bldg_pop['block_area']
    return bldg_pop


def add_block_bldg_count_density(bldg_pop: gpd.GeoDataFrame,
                                 block: gpd.GeoDataFrame,
                                 ) -> gpd.GeoDataFrame:
    """
    Calculates the ratio of number of buildings in a block to the block's area and adds that
    to the bldg_pop geodf
    """
    bldg_pop['block_building_count_density'] = bldg_pop['building_count'] / bldg_pop['block_area']
    return bldg_pop


def add_block_pop(bldg_pop_ls: gpd.GeoDataFrame,
                  bldg_pop_wp: gpd.GeoDataFrame,
                  ) -> gpd.GeoDataFrame:
    """
    Calculates the population for the block and adds that to the bldg_pop geodf
    """
    block_pop_ls = bldg_pop_ls[['block_id', 'bldg_pop']].groupby('block_id').sum()
    block_pop_wp = bldg_pop_wp[['block_id', 'bldg_pop']].groupby('block_id').sum()
    block_pop_ls.rename(columns={'bldg_pop': 'block_pop_ls'}, inplace=True)
    block_pop_wp.rename(columns={'bldg_pop': 'block_pop_wp'}, inplace=True)
    bldg_pop = bldg_pop_ls.merge(block_pop_ls, how='left', on='block_id')
    bldg_pop = bldg_pop.merge(block_pop_wp, how='left', on='block_id')
    return bldg_pop


def add_block_pop_density(bldg_pop: gpd.GeoDataFrame,
                          block: gpd.GeoDataFrame,
                          ) -> gpd.GeoDataFrame:
    """
    Calculates the ratio of block population to block area and adds that to the bldg_pop geodf
    """    
    bldg_pop['block_pop_density'] = bldg_pop['block_pop'] / bldg_pop['block_area']
    return bldg_pop  


######################################
# COMMANDS FOR GENERAL AOI SUMMARIES #
######################################
def make_superblock_summary(bldg_pop_data_ls: gpd.GeoDataFrame, 
                            bldg_pop_data_wp: gpd.GeoDataFrame, 
                            block_data: gpd.GeoDataFrame,
                            aoi_out_path: str = None,
                            ) -> None:
    '''
    Calculates all statistics given:
        1. bldg-level pop allocation
        2. block geometry
        3. path to save output to
    '''

    if isinstance(bldg_pop_data_ls, gpd.GeoDataFrame):
        bldg_pop_ls = bldg_pop_data_ls
    else:
        bldg_pop_ls = load_bldg_pop(bldg_pop_data_ls)

    if isinstance(bldg_pop_data_wp, gpd.GeoDataFrame):
        bldg_pop_wp = bldg_pop_data_wp
    else:
        bldg_pop_wp = load_bldg_pop(bldg_pop_data_wp)


    if 'block_id' not in bldg_pop_ls.columns:
         bldg_pop_ls = add_block_id(bldg_pop_ls, block_data)
    
    if 'block_id' not in bldg_pop_wp.columns:
         bldg_pop_wp = add_block_id(bldg_pop_wp, block_data)

    bldg_pop_ls = set_dtypes(bldg_pop_ls, block_data)
    bldg_pop = add_block_pop(bldg_pop_ls, bldg_pop_wp)
    bldg_pop['block_pop_ls'] = bldg_pop['block_pop_ls'].apply(lambda x: np.nan if x == 0 else x)
    bldg_pop['block_pop_wp'] = bldg_pop['block_pop_wp'].apply(lambda x: np.nan if x == 0 else x)


    if aoi_out_path is not None:
        aoi_out_path = Path(aoi_out_path)
        aoi_out_path.parent.mkdir(parents=True, exist_ok=True)
        bldg_pop.to_file(str(aoi_out_path), driver='GeoJSON')
    return bldg_pop



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



def make_summary(superblock: Union[str, Path, gpd.GeoDataFrame],
                 landscan_path: Union[str, Path],
                 superblock_buildings: Union[str, Path, gpd.GeoDataFrame],
                 summary_out_path: Union[str, Path] = None,
                 ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Creates a summary of the given area using its population, blocks,
    and boundary lines. 
    """
    # (1) Read in the superblock data if provided a path, or pass over this
    #     logic if the superblock data is already in memory
    if not isinstance(superblock, gpd.GeoDataFrame):
        if isinstance(superblock, str):
            superblock_path = Path(superblock)
        elif isinstance(superblock, Path):
            superblock_path = superblock
        else: 
            raise Exception(f'Unknown value of type {type(superblock)} passed as superblock')
        superblock = gpd.read_file(superblock_path)

    if not isinstance(superblock_buildings, gpd.GeoDataFrame):
        if isinstance(superblock_buildings, str):
            superblock_buildings_path = Path(superblock_buildings)
        elif isinstance(superblock_buldings, Path):
            superblock_buildings_path = superblock_buildings
        else: 
            raise Exception(f'Unknown value of type {type(superblock_buildings)} passed as superblock_buildings')
        superblock_buildings = gpd.read_file(superblock_buildings_path)


    # (2) Allocate Landscan
    country_code = superblock['block_id'][0].split('.')[0]
    gadm_list = list(set(superblock['gadm_code']))
    world_pop_path = '/project2/bettencourt/mnp/analytics/data/population/WorldPop_tifs/' + country_code +'.tif'
    _, superblock_ls = extract_aoi_data_from_raster(superblock, str(landscan_path), save_geojson=False, save_tif=False)
    _, superblock_wp = extract_aoi_data_from_raster(superblock, world_pop_path, save_geojson=False, save_tif=False)
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
    bldg_pop_alloc_ls = allocate_population(superblock_buildings, superblock_ls, 'pop')
    bldg_pop_alloc_wp = allocate_population(superblock_buildings, superblock_wp, 'pop')

    # (2) Now assemble the other data
    superblock_bldg_summary = make_superblock_summary(bldg_pop_alloc_ls, bldg_pop_alloc_wp, superblock)
    block_cols = [x for x in superblock_bldg_summary.columns if "block" in x]
    superblock_stats = superblock_bldg_summary[block_cols].drop_duplicates()
    superblock_summary = superblock.merge(superblock_stats, how='left', on='block_id')

    # (3) Save
    summary_out_path = Path(summary_out_path)
    fname = summary_out_path.stem
    outdir = summary_out_path.parent
    outdir.mkdir(exist_ok=True, parents=True)

    superblock_summary = remove_duplicated_cols_from_merge(superblock_summary)
    superblock_buildings_out_path = outdir / (fname + "-bldgs" + summary_out_path.suffix)

    # Not writing because it gets written in batch_block
    # if summary_out_path.suffix == '.geojson':
    #     superblock_summary.to_file(summary_out_path, driver='GeoJSON')
    #     superblock_bldg_summary.to_file(superblock_buildings_out_path, driver='GeoJSON')
    # elif summary_out_path.name.endswith('parquet'):
    #     parquet_write(superblock_summary, summary_out_path)
    #     parquet_write(superblock_bldg_summary, superblock_buildings_out_path)
    # else:
    #     raise Exception(f'Out path file format {superblock_buildings_out_path.suffix} not recognized')
    # print("Saved to: {}".format(str(summary_out_path)))
    # print("Saved to: {}".format(str(superblock_buildings_out_path)))

    return superblock_summary, superblock_bldg_summary




def load_raster_selection(raster_io: rasterio.io.DatasetReader,
                          geom_list: List[ShapelyGeom],
                          ) -> np.ndarray:
    '''
    rasterio allows loading of subselection of a tiff file, so 
    given a list of geometries and an open raster DatasetReader,
    loads in the values within the geom_list.
    This allows for loading of small selections from otherwise 
        huge tiff's
    '''

    if not isinstance(geom_list, gpd.GeoSeries):
        geom_list = [geom_list]

    # Find the window around the geom_list
    geom = [geom_list.unary_union]
    window = rasterio.features.geometry_window(raster_io, geom)
    # window = rasterio.features.geometry_window(raster_io, geom_list)
    transform = raster_io.window_transform(window)

    # Perform the windowed read
    sub_data = raster_io.read(window=window)
    return sub_data, transform, window 


def save_np_as_geotiff(np_array: np.ndarray, 
                       transform: affine.Affine,
                       out_file: Path,
                       crs: CRS = CRS.from_epsg(4326),
                       ) -> None:
    '''
    Given a numpy array of data, and transform and crs defining
    the coordinate system, save as a geotiff
    '''

    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    c, h, w = np_array.shape

    new_file = rasterio.open(
        out_file,
        'w',
        driver='GTiff',
        height=h,
        width=w,
        count=c,
        dtype=np_array.dtype,
        crs=crs,
        transform=transform)

    #print("shape = {}".format(np_array.shape))
    new_file.write(np_array)
    print("Saved GeoTiff at: {}".format(out_file))


def raster_to_geodataframe(raster_data: np.ndarray,
                           transform: affine.Affine,
                           window: rasterio.windows.Window,
                           raster_io: rasterio.io.DatasetReader,
                           ) -> gpd.GeoDataFrame:
    """
    Takes in raster data as a numpy array and converts it to a geo df
    """
    if raster_data.ndim == 3:
        assert raster_data.shape[0] == 1, "ERROR - can't handle multichannel yet"
        raster_data = raster_data[0]

    # To convert when our raster_data is from a windowed
    # read we need to go from the windowed raster_data 
    # coords -> orig dataset coords -> x,y spatial
    data_h, data_w = raster_data.shape 
    all_geoms = []
    data = []
    for h in range(data_h):
        for w in range(data_w):
            # Convert to the coord sys of the raster_io dataset
            abs_h = h + window.row_off 
            abs_w = w + window.col_off 

            # And use its xy method now
            # abs_x1, abs_y1 = raster_io.xy(abs_h, abs_w)
            # abs_x0, abs_y0 = raster_io.xy(abs_h-1, abs_w-1)
            p_ul = Point(raster_io.xy(abs_h, abs_w, offset='ul'))
            p_lr = Point(raster_io.xy(abs_h, abs_w, offset='lr'))
            geom = MultiPoint([p_ul, p_lr]).envelope
            all_geoms.append(geom)
            data.append(raster_data[h,w])
    return gpd.GeoDataFrame.from_dict({'geometry': all_geoms,
                                       'pop': data
                                       })


def extract_aoi_data_from_raster(geometry_data: Union[str, gpd.GeoDataFrame], 
                                 raster_path: str, 
                                 out_path: str = None, 
                                 save_geojson: bool = True, 
                                 save_tif: bool = True,
                                 ) -> Tuple:
    '''
    Extracts the relevant data from a larger raster tiff, per the geometries
    in geometry_data file and saves out.

    Inputs:
        - geometry_data (str) geojson file of vector geometries
        - raster_path (str) tiff file of raster data
        - out_path (str) path to save extracted data
    '''
    out_path = Path(out_path) if out_path is not None else out_path

    raster_io = rasterio.open(raster_path)
    if isinstance(geometry_data, gpd.GeoDataFrame):
        aoi_geoms = geometry_data
    else:
        if geometry_data.suffix == ".csv":
            aoi_geoms = load_csv_to_geo(geometry_data)
        else:
            aoi_geoms = gpd.read_file(geometry_data)

    raster_data_selection, transform, window = load_raster_selection(raster_io, 
                                                    aoi_geoms['geometry'])
    
    gdf_data = raster_to_geodataframe(raster_data_selection, 
                                      transform,
                                      window,
                                      raster_io)

    gdf_data['pop'] = gdf_data['pop'].apply(lambda x: max(0, x))

    if save_tif:
        tif_path = out_path.with_suffix('.tiff')
        save_np_as_geotiff(raster_data_selection, transform, str(tif_path))
    if save_geojson:
        geojson_path = out_path.with_suffix('.geojson')
        gdf_data.to_file(geojson_path, driver='GeoJSON')
        print("Saved GeoJSON at: {}".format(geojson_path))
    return raster_data_selection, gdf_data


def fix_invalid_polygons(geom: Union[Polygon, MultiPolygon]) -> Polygon:
    """
    Fix self-intersection polygons
    """
    if geom.is_valid:
        return geom 
    else:
        return geom.buffer(0)


def allocate_population(buildings_gdf: gpd.GeoDataFrame,
                        population_gdf: gpd.GeoDataFrame,
                        pop_variable: str,
                        ) -> gpd.GeoDataFrame:
    """

    """
    for k in ['index_right', 'index_left']:
        for df in [buildings_gdf, population_gdf]:
            if k in df.columns:
                df.drop(columns=[k], inplace=True)

    # Map each building to the pop geom it is in
    in_cols = list(buildings_gdf.columns)
    population_gdf['pop_id'] = np.arange(population_gdf.shape[0])
    buildings_gdf['bldg_id'] = np.arange(buildings_gdf.shape[0])
    geo = gpd.sjoin(buildings_gdf, population_gdf,
                    how='left', op='intersects')[['bldg_id', 'geometry', 'pop_id', pop_variable]]
    
    # Handle building geoms intersecting multiple pop geoms
    # Split building into each piece, allocate, then sum by bldg_id
    geo = geo.join(population_gdf[['pop_id', 'geometry']], on='pop_id', how='left', rsuffix='_pop')
    geo['obs_count'] = 1
    bldg_count = geo[['bldg_id', 'obs_count']].groupby('bldg_id').sum()
    bldg_count.reset_index(inplace=True)
    geo.drop(columns=['obs_count'], inplace=True)
    geo = geo.merge(bldg_count, on='bldg_id', how='left')

    null_geo = geo[geo['geometry_pop'].isna()]
    print(f'Superblock had {len(null_geo)} null geometries.')
    geo = geo[~geo['geometry_pop'].isna()]

    fn = lambda s: s['geometry'].intersection(s['geometry_pop'])
    geo['unique_geom'] = geo.apply(fn, axis=1)
    geo.set_geometry('unique_geom', inplace=True)

    # Numerator is the buildings area
    geo['num_area'] = geo.geometry.area

    # Denom is the area of all buildings in that pop_id
    geo_by_pop_id = geo[['pop_id', 'num_area']].groupby('pop_id').sum()
    geo_by_pop_id.rename(columns={'num_area': 'den_area'}, inplace=True)
    geo_by_pop_id.reset_index(inplace=True)
    
    # Merge the denom and generate the factor and allocation pop
    geo = geo.merge(geo_by_pop_id, on='pop_id', how='left')
    geo['alloc_factor'] = geo['num_area'] / geo['den_area']
    geo['bldg_pop'] = geo['alloc_factor'] * geo[pop_variable]
    
    geo = gpd.GeoDataFrame(pd.concat([geo, null_geo], ignore_index=True), crs=geo.crs)

    # Because we split the bldg geoms w > 1 intersection, sum
    # up by bldg_id to reassemble
    bldg_pop = geo[['bldg_pop', 'bldg_id']].groupby('bldg_id').sum().reset_index()
    assert bldg_pop.shape[0] == buildings_gdf.shape[0], "ERROR - bldg count IN = {} but bldg count OUT = {}".format(buildings_gdf.shape[0], bldg_pop.shape[0] )
    if 'bldg_id' not in in_cols:
        in_cols.append('bldg_id')
    buildings_gdf = buildings_gdf[in_cols].merge(bldg_pop, on='bldg_id', how='left')

    return buildings_gdf


def alloc_pop_to_buildings(buildings_gdf_path: str,
                           master_pop_path: str,
                           alloc_out_path: str,
                           pop_out_path: str,
                           ) -> None:
    '''
    Extracts the relevant AoI in the master pop file based on the
    buildings gdf, saves that out if you want, then does a building-level
    allocation based on the building area. Saves out that building-level
    population allocation

    Inputs:
        - buildings_gdf_path: path to buildings gdf
        - master_pop_path: path to Landscan raster pop data
        - alloc_out_path: path to save final building-level pop allocation
        - pop_out_path: path to save extracted data from global Landscan pop
    '''
    buildings_gdf_path = Path(buildings_gdf_path)
    master_pop_path = Path(master_pop_path)
    pop_out_path = Path(pop_out_path)

    extract_aoi_data_from_raster(buildings_gdf_path, 
                                 master_pop_path, 
                                 pop_out_path,
                                 ) 
    pop_out_path = pop_out_path.with_suffix(".geojson")
    alloc_out_path = Path(alloc_out_path)
    population_gdf = gpd.read_file(pop_out_path)
    pop_variable = 'pop'
    buildings_gdf = gpd.read_file(buildings_gdf_path)

    bldg_pop = allocate_population(buildings_gdf, population_gdf, pop_variable)
    alloc_out_path.parent.mkdir(exist_ok=True, parents=True)
    bldg_pop.to_file(alloc_out_path, driver='GeoJSON')




def join_block_building(block_gdf: gpd.GeoDataFrame,
                        buildings_gdf: gpd.GeoDataFrame,
                        ) -> gpd.GeoDataFrame:
    """
    A helper function that merges a blocks geodataframe
    with a buildings geodataframe. Returns just the buildings
    geo df.
    """
    buildings_gdf = gpd.sjoin(buildings_gdf, block_gdf,
                              how='left', op='intersects')
    return buildings_gdf


def load_blocks_buildings(block_path: str, 
                          building_path: str,
                          merge_bldgs=False,
                          ) -> Tuple[gpd.GeoDataFrame]:
    """
    Loads in the blocks and buildings geopandas dataframes, optionally
    adding the blocks dataframe onto the buildings dataframe. Returns
    both no matter what. 
    """    
    buildings_gdf = gpd.read_file(building_path)
    blocks_gdf = gpd.read_file(block_path)

    if merge_bldgs:
        buildings_gdf = join_block_building(blocks_gdf, buildings_gdf)
    return blocks_gdf, buildings_gdf


def load_csv_to_geo(csv_path: str, 
                    add_file_col: bool = False,
                    crs: str = "EPSG:4326",
                    ) -> gpd.GeoDataFrame:
    """
    Loads a csv with geometries, turns it into a geo df, and applies
    the proper CRS. Returns the geo df.
    """
    df = pd.read_csv(csv_path, usecols=['block_id', 'geometry'])
    df['geometry'] = df['geometry'].apply(loads)
    gdf = gpd.GeoDataFrame(df, geometry='geometry')
    gdf.crs = crs
    gdf['geometry'].crs = crs
    return gdf


def remove_duplicated_cols_from_merge(block_data: gpd.GeoDataFrame):
    block_data['block_area'] = block_data['block_area_x']
    block_data = block_data.drop(['block_area_x', 'block_area_y'], axis='columns')
    return block_data


def parquet_write(block_data: gpd.GeoDataFrame,
                  output_path: Path):
    assert block_data.crs == 4326
    block_data.to_parquet(output_path)


def parquet_read(input_path: Union[str, Path]):
    gdf = gpd.read_parquet(input_path)
    return gdf
