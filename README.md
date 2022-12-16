
## kblock
Python tools for creating residential street block delineations, estimating population, and measuring access of buildings to street networks.

### Set up data and environment

#### Clone repository
```
mkdir -p /users/project/repo
cd /users/project/repo
git clone git@github.com:mansueto-institute/kblock.git
```

#### Create conda environment
##### Build environment from scratch.
```
conda update conda
conda update -n base -c defaults conda

conda env remove -n geospatial
conda create --name geospatial python=3.9.7 --yes

source activate geospatial

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

##### Or install directly from conda yml.
```
conda env update --name geospatial --file /users/project/repo/kblock/geospatial.yml
```

#### Download data
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

##### Download the Natural Earth Countries
https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip


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



### Local deployment

#### Run each bash file serially from terminal.
```
cd /users/project/repo/kblock/kblock
source activate geospatial
bash /users/project/repo/kblock/deploy/1-prepare-blocks/deploy_1a_prepare_gadm.sh
bash /users/project/repo/kblock/deploy/1-prepare-blocks/deploy_1b_generate_blocks.sh
bash /users/project/repo/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings.sh
bash /users/project/repo/kblock/deploy/3-model-population/deploy_3_model_population.sh
bash /users/project/repo/kblock/deploy/4-compute-k/deploy_4_compute_k.sh
```

### Midway HPC deployment

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
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings_bigmem.sbatch
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/2-centroid-buildings/deploy_2_prepare_buildings.sbatch
sbatch project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/3-model-population/deploy_3_model_population.sbatch
bash project2/<pi-cnetid>/projects/mnp/repos/kblock/deploy/4-compute-k/deploy_4_job_trigger.sh
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



