
## kblock
Python tools for creating residential street block delineations, estimating population, and measuring access of buildings to street networks. See [https://www.millionneighborhoods.africa](https://www.millionneighborhoods.africa) for data and [https://github.com/mansueto-institute/kblock-analysis](https://github.com/mansueto-institute/kblock-analysis) for analysis.

### Set up data and environment

#### Create conda environment
##### Build environment from scratch.
```
conda update conda
conda update -n base -c defaults conda

conda env remove -n kblock_env
conda create --name kblock_env python=3.9.7 --yes

source activate kblock_env

conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
conda install -c conda-forge geopandas=0.12.1 --yes
conda install -c conda-forge pygeos=0.12.0 --yes
conda install -c conda-forge urlpath --yes
conda install -c conda-forge dask --yes
conda install -c conda-forge dask-geopandas --yes
conda install -c conda-forge pyarrow --yes
conda install -c conda-forge rasterio --yes
conda install -c conda-forge xarray --yes
conda install -c conda-forge rioxarray --yes
conda install -c conda-forge momepy --yes
conda install -c conda-forge pygeohash --yes
conda install -c conda-forge pyogrio --yes
conda install -c conda-forge ipykernel --yes
```

##### Or install environment directly from repo
```
mkdir -p /users/project/repo
cd /users/project/repo
git clone git@github.com:mansueto-institute/kblock.git
cd /users/project/repo/kblock

conda env create --name kblock_env --file geospatial.yml --force
conda activate kblock_env
```

## Download data
```
brew install wget
mkdir -p /users/project/inputs
```
##### Download [coastlines](https://daylightmap.org/coastlines.html) from daylightmap.org
```
wget -O users/project/inputs/coastlines.tgz https://daylight-map-distribution.s3.us-west-1.amazonaws.com/release/v1.20/coastlines-v1.20.tgz
```
##### Download [GHLS Urban Centre Database](https://data.jrc.ec.europa.eu/dataset/53473144-b88c-44bc-b4a3-4583ed1f547e#dataaccess) from the EU Joint Research Centre
```
wget -O users/project/inputs/ghsl.zip https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_STAT_UCDB2015MT_GLOBE_R2019A/V1-2/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.zip
```
##### Run OSM and GADM download job (from geofabrik.de and gadm.org)
```
brew install osmium-tool
cd /users/project/kblock/deploy/0-ingest-data/
bash deploy_download.sh -o users/project/inputs
```

##### Download [LandScan 2020 1 kilometer grid](https://landscan.ornl.gov/) from Oak Ridge National Laboratory
`mkdir -p users/project/inputs/rasters/landscan` and add file `/landscan-global-2020.tif`
##### Download [WorldPop 2020 Constrained 100 meter grid](https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/)
```
mkdir -p users/project/inputs/rasters/worldpop
outpath=users/project/inputs/rasters/worldpop

iso_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)         
for c in ${iso_list[@]}; do
	echo "$c"
	d=$(echo "$c" | tr '[:upper:]' '[:lower:]') 
	wget -O $outpath/${c}.tif https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/${c}/${d}_ppp_2020_UNadj_constrained.tif
done
```
##### Download [Ecopia Landbase Africa powered by Maxar](https://platform.ecopiatech.com/login) (requires gaining access to DigitizeAfrica and manually downloading files)
```
mkdir -p users/project/inputs/ecopia
```

##### Download [UN Population](https://population.un.org/wpp/Download/Standard/CSV/) data 
```
wget "https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2022_Demographic_Indicators_Medium.zip" -O temp.zip
unzip temp.zip
rm temp.zip
```

## Local deployment

### Expected output from each job
* `deploy_1a_prepare_gadm.sh` prepares the GADM delineations (cleaning the boundaries to align with coastlines and remove water features from land area.
* `deploy_1b_generate_blocks.sh` generates the block level geometries from OSM street networks and natural features (and GADM boundaries as well).
* `deploy_1c_regions_crosswalk.sh` aligns block geometries with GHSL and Africapolis urban boundaries
* `deploy_2_prepare_buildings.sh` takes the Ecopia files and converts the building polygons to points and calculates building areas.
* `deploy_3_model_population.sh` implements a dasymetric population model to estimate block level population using raster grids (LandScan and WorldPop) and building area.
* `deploy_4_compute_k.sh` computes block complexity statistics (and other block level properties)
* `deploy_5_combine_data.sh` combines data and computes additional block statistics


### Minimal reproducible example
#### Download data [here](dsbprylw7ncuq.cloudfront.net/_sampledata/sample-data.zip) (recommended) (or [here](https://drive.google.com/drive/folders/1Cs9RK01hltsw9ZK-y3dDU_5LXclXNx0K?usp=sharing))
* Unzip `sample-data.zip` into `/users/downloads/sample-data`

#### Copy prepared buildings and prepared land polygons (this uses OSM buildings for demonstration purposes)
```
cd /users/downloads/sample-data
cp -r _minreprex/gadm/parquet inputs-reprex/gadm
cp -r _minreprex/buildings/osm outputs-reprex/buildings
```

#### Clone repo and cd into it
```
cd /users/desktop
git clone https://github.com/mansueto-institute/kblock.git
```

#### Install environment
```
cd /users/desktop/kblock
conda env create --name kblock_env --file geospatial.yml --force
conda activate kblock_env
```

#### Generate blocks, model population, compute block complexity
```
bash ./deploy/1-prepare-blocks/deploy_1b_generate_blocks.sh /users/downloads/sample-data
bash ./deploy/3-model-population/deploy_3_model_population.sh /users/downloads/sample-data
bash ./deploy/4-compute-k/deploy_4_compute_k.sh /users/downloads/sample-data
```

A more extended reprex with additional steps is available [here](https://gist.github.com/nmarchio/8dc8aa5736202da32cf02c4c2da0e85b).

## Midway HPC deployment

#### Transfer files to HPC
```
rsync -avz /users/project/inputs/daylight/coastlines-v1.19 <cnetid>@midway.rcc.uchicago.edu:/project2/<pi-cnetid>/projects/mnp/inputs/daylight
rsync -avz /users/project/inputs/gadm/geojson <cnetid>@midway.rcc.uchicago.edu:/project2/<pi-cnetid>/projects/mnp/inputs/gadm
rsync -avz /users/project/inputs/osm/parquet <cnetid>@midway.rcc.uchicago.edu:/project2/<pi-cnetid>/projects/mnp/inputs/osm
rsync -avz /users/project/inputs/ghsl <cnetid>@midway.rcc.uchicago.edu:/project2/bett<pi-cnetid>/projects/mnp/inputs/
rsync -avz /users/project/inputs/rasters <cnetid>@midway.rcc.uchicago.edu:/project2/<pi-cnetid>/projects/mnp/inputs/
rsync -avz /users/project/inputs/ecopia <cnetid>@midway.rcc.uchicago.edu:/project2/<pi-cnetid>/projects/mnp/inputs/
```

#### Deploy jobs on SLURM

##### Clone repo and set up environment.
```
ssh <cnetid>@midway2.rcc.uchicago.edu
cd /project2/<pi-cnetid>/projects/mnp/repo/
git clone git@github.com:mansueto-institute/kblock.git
module load python/anaconda-2021.05
conda env update --name geospatial --file ./kblock/geospatial.yml
source activate geospatial
```

##### Run each sbatch files serially from terminal on login node (takes a few days for each job).
```
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/1-prepare-blocks/deploy_1a_prepare_gadm.sbatch
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/1-prepare-blocks/deploy_1b_generate_blocks.sbatch
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/1-prepare-blocks/deploy_1c_regions_crosswalk.sbatch
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings_bigmem.sbatch
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings.sbatch
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/3-model-population/deploy_3_model_population.sbatch
bash project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/4-compute-k/deploy_4_job_trigger.sh
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/5-combine-data/deploy_5_combine_data.sbatch
```

```
cd /project2/bettencourt/mnp/update
cd /Users/nm/Desktop/kblock/kblock
source activate geospatial
bash /Users/nm/Desktop/kblock/deploy/1-prepare-blocks/deploy_1a_prepare_gadm.sh
bash /Users/nm/Desktop/kblock/deploy/1-prepare-blocks/deploy_1b_generate_blocks.sh
bash /Users/nm/Desktop/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings.sh
bash /Users/nm/Desktop/kblock/deploy/3-model-population/deploy_3_model_population.sh
bash /Users/nm/Desktop/kblock/deploy/4-compute-k/deploy_4_compute_k.sh
bash /users/project/repo/kblock/deploy/5-combine-data/deploy_5_combine_data.sh
```

```
cd /project2/bettencourt/mnp/update/repo
git clone git@github.com:mansueto-institute/kblock.git
module load python/anaconda-2021.05
source activate geospatial

cd /project2/bettencourt/mnp/update/repo/kblock/kblock
source activate geospatial

sbatch project2/bettencourt/mnp/update/repo/kblock/deploy/1-prepare-blocks/deploy_1a_prepare_gadm.sbatch
sbatch project2/bettencourt/mnp/update/repo/kblock/deploy/1-prepare-blocks/deploy_1b_generate_blocks.sbatch
sbatch project2/bettencourt/mnp/update/repo/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings_bigmem.sbatch
sbatch project2/bettencourt/mnp/update/repo/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings.sbatch
sbatch project2/bettencourt/mnp/update/repo/kblock/deploy/3-model-population/deploy_3_model_population.sbatch
bash project2/bettencourt/mnp/update/repo/kblock/deploy/4-compute-k/deploy_4_job_trigger.sh
```



