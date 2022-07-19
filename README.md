# kblock
Python tools for building and measuring street blocks using geographic data

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
