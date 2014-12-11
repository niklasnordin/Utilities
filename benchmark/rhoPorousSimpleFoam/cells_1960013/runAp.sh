#!/bin/bash

#PBS -N Rho_1.9M
#PBS -l walltime=0:30:00
#PBS -l mppwidth=1440
NSLOTS=1440

# set case and code
# =================================
FOAM_ROOT=/cfs/scania/home/n/nnordi/OpenFOAM/nnordi-2.2.0/run
CASE_DIR=$FOAM_ROOT/benchmark/rhoPorousSimpleFoam/cells_1960013/ref
FOAM_EXE=rhoPorousSimpleFoam

# *********************************************************
# ******** usually no need to modify anything below *******
# *********************************************************

# set the environment
# =================================
. /opt/modules/default/etc/modules.sh
module load PrgEnv-gnu
. /cfs/scania/home/n/nnordi/OpenFOAM/OpenFOAM-2.2.0/etc/bashrc

LOGFILE=$CASE_DIR/$FOAM_EXE.log
FOAM_CODE=$FOAM_APPBIN/$FOAM_EXE

MYSYS=/cfs/scania/home/n/nnordi/OpenFOAM/ThirdParty-2.0.1/system
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MYSYS/alps/lib:$MYSYS/gcc-4.5.1/lib64:$MYSYS/mpich/lib:$MYSYS/pmi/lib64:$MYSYS/udreg/lib64:$MYSYS/ugni/lib64:$MYSYS/xpmem/lib64

# run the case
# =================================
aprun -n $NSLOTS $FOAM_EXE -case $CASE_DIR -parallel >& $CASE_DIR/$FOAM_EXE.log.$NSLOTS
