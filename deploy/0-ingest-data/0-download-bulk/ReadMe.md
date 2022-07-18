# OSM and GADM Data Ingestion 

## Dependencies for download script (MacOS)
#### OSM extraction setup (requires homebrew at https://brew.sh/)
```
brew install wget
brew install osmium-tool
```
#### Python environment setup (requires miniconda at https://docs.conda.io/en/latest/miniconda.html)
```
conda create --name geospatial python=3.9.7 --yes
source activate geospatial

conda install -c conda-forge pygeos=0.10.2 --yes
conda install -c conda-forge geopandas=0.10.2 --yes 
conda install -c conda-forge urlpath --yes
conda install -c conda-forge requests --yes

conda install -c conda-forge dask --yes
conda install -c conda-forge dask-geopandas --yes
conda install -c conda-forge pyarrow --yes

conda install -c conda-forge rasterio --yes
conda install -c conda-forge xarray --yes
conda install -c conda-forge rioxarray --yes
conda install -c conda-forge momepy --yes
```

## Steps to download OSM linestrings and GADMs boundaries and write to local directory

#### 1. Activate environment, clone repo, cd into repo, make output directory
```
source activate geospatial
cd /users/repos/
git clone git@github.com:mansueto-institute/kblock.git
cd ./kblock/deploy/0-ingest-data/0-download-bulk
mkdir -p /users/downloads/folder
```
#### 2. Deploy download job (pass in output directory with -o argument)
> **Caution:** Running this script will download 20GB for all Sub-Saharan countries in the `gadm_list=()` in the `geofabrik_list=()` in `download_job.sh`.
```
bash download_job.sh -o /users/downloads/folder
```
