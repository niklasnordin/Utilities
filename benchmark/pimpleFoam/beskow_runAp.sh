#!/bin/bash -l
# The -l above is required to get the full environment with modules

#The name of the script
#SBATCH -J nPimple_916k

# 24 hour wall-clock time will be given to this job
#SBATCH -t 1:00:00

# Number of nodes
#SBATCH -N 3

# Number of MPI processes per node (default 40)
##SBATCH --ntasks-per-node=40
#SBATCH --ntasks-per-node=32


NSLOTS=96

# set case and code
# =================================
FOAM_ROOT=/cfs/scania/home/n/nnordi/OpenFOAM/nnordi-2.3.x/run
CASE_DIR=$FOAM_ROOT/benchmark/pimpleFoam/cells_916801
FOAM_EXE=pimpleFoam

# *********************************************************
# ******** usually no need to modify anything below *******
# *********************************************************

# set the environment
# =================================
##. /opt/modules/default/etc/modules.sh
#module add PrgEnv-intel
module add PrgEnv-gnu

LOGFILE=$CASE_DIR/$FOAM_EXE.log
FOAM_CODE=$FOAM_APPBIN/$FOAM_EXE

##MYSYS=/cfs/milner/scratch/n/nnordi/OpenFOAM/ThirdParty-2.0.1/system
##export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MYSYS/alps/lib:$MYSYS/gcc-4.5.1/lib64:$MYSYS/mpich/lib:$MYSYS/pmi/lib64:$MYSYS/udreg/lib64:$MYSYS/ugni/lib64:$MYSYS/xpmem/lib64

# run the case
# =================================
aprun -n $NSLOTS $FOAM_EXE -case $CASE_DIR -parallel >& $CASE_DIR/$FOAM_EXE.log
