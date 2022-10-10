#!/bin/bash
module purge
module use -a /app/spack/2022.01/modulefiles
module load gcc/11.2.0
#module load netcdf/4.2  
module load netcdf-c/4.8.1
module load netcdf-fortran/4.5.3

./examples/rfmip-clear-sky/clearsky_rrtmgp test.nc\
                           ./rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc \
                           ./rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
