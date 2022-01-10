
### PyrOSM Environment 

Due to performance issues the [PyrOSM](https://pyrosm.readthedocs.io/en/latest/) library utilized in the `download.py` module is not recommended for large countries or large scale data processing.
However, this may be convenient for local development use cases. 

```
conda create --name pyrosm_env python=3.9.7 --yes
source activate pyrosm_env
conda install -c conda-forge pygeos=0.12.0 --yes
conda install -c conda-forge geopandas=0.10.2 --yes
conda install -c conda-forge cykhash=1.0.2 --yes
conda install -c conda-forge pyrosm=0.6.1 --yes
conda install -c conda-forge urlpath --yes
conda install -c conda-forge pyarrow --yes
conda install -c conda-forge dask --yes
conda install -c conda-forge dask-geopandas --yes
```
