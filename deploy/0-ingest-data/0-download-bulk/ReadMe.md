# OSM and GADM Data Ingestion 

## Dependencies for download script (MacOS)
#### OSM extraction setup (requires homebrew at https://brew.sh/)
```
brew install wget
brew install osmium-tool
```

## Steps to download OSM linestrings and GADMs boundaries and write to local directory

#### 1. Clone repo
```
cd /users/repos/
git clone git@github.com:mansueto-institute/kblock.git
```
#### 2. Activate environment and deploy job (pass in output directory with -o argument)
> **Caution:** Running this script will download 20GB for all Sub-Saharan countries in the `gadm_list=()` in the `geofabrik_list=()` in `download_job.sh`.
```
source activate geospatial
cd ./kblock/deploy/0-ingest-data/0-download-bulk
mkdir -p /users/downloads/folder
bash download_job.sh -o /users/downloads/folder
```
