### Build on Midway
```
module load python/anaconda-2021.05
conda create --name build_block_env python=3.9.7 --yes
source activate build_block_env
conda install -c conda-forge pygeos=0.12.0 --yes
conda install -c conda-forge geopandas=0.10.2 --yes
conda install -c conda-forge cykhash=1.0.2 --yes
conda install -c conda-forge pyrosm=0.6.1 --yes
conda install -c conda-forge urlpath --yes
conda install -c conda-forge pyarrow --yes
conda install -c conda-forge dask --yes
conda install -c conda-forge dask-geopandas --yes
```

### Compute on Midway
```
module load python/anaconda-2021.05
conda create --name kcompute_env python=3.9.7 --yes
source activate kcompute_env
conda install -c conda-forge pygeos=0.12.0 --yes
conda install -c conda-forge geopandas=0.10.2 --yes
conda install -c conda-forge urlpath --yes
conda install -c conda-forge pyarrow --yes
conda install -c conda-forge dask --yes
conda install -c conda-forge dask-geopandas --yes
```

### Local deploy
```
cd /Users/nm/Desktop/Projects/work/mnp/kblock/dev/compute
source activate py395_pygeos
bash /Users/nm/Desktop/Projects/work/mnp/kblock/dev/compute/deploy_local.sh	
conda deactivate 
```

### Login and update repo
```
ssh <cnetid>@midway2.rcc.uchicago.edu
cd /project2/bettencourt/mnp/analytics/scripts/kblock
git pull
```

### Midway deploy
```
cd /project2/bettencourt/mnp/analytics/scripts/kblock/dev/compute/
module load python/anaconda-2021.05
source activate py395_pygeos
sbatch /project2/bettencourt/mnp/analytics/scripts/kblock/dev/compute/deploy_midway.sbatch 
```

### Check logs and outputs
```
less /project2/bettencourt/mnp/analytics/deployments/log_output_sle.log
less /project2/bettencourt/mnp/analytics/deployments/sle_job.err 
ls /project2/bettencourt/mnp/analytics/outputs/SLE

conda deactivate 
```
