#!/bin/bash

module purge
module use -a /app/spack/2022.01/modulefiles
module load gcc/11.2.0
module load netcdf-c/4.8.1
module load netcdf-fortran/4.5.3

export NCHOME="/app/spack/v0.17/linux-rhel7-x86_64/gcc-11.2.0/netcdf-c/4.8.1-kd2fizpbah22hkghjxhhzkwh2ekxetdj"
export NFHOME="/app/spack/v0.17/linux-rhel7-x86_64/gcc-11.2.0/netcdf-fortran/4.5.3-m4arkb34bmowqtb4br3xgceig3x7tovf"
export RRTMGP_DIR="$PWD/build"

export FC=gfortran
export FCFLAGS="-g -O3 -I$NFHOME -I$RRTMGP_DIR -ffree-line-length-none"
export LD=gfortran

cd build
make
cd -

export RRTMGP_ROOT="/local2/home/projects/rte-rrtmgp/"

cd examples/rfmip-clear-sky
make clean; make
