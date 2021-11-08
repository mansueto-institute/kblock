### Build on Midway
```
module load python/anaconda-2021.05
conda create --name py395_pygeos python=3.9.5 --channel conda-forge
source activate py395_pygeos
conda install pygeos=0.11.1 --channel conda-forge --yes
conda install geopandas=0.10.2 --channel conda-forge --yes
conda install pyarrow --channel conda-forge --yes
conda install dask --channel conda-forge --yes
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
