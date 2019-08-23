#!/bin/bash

# Begin LSF Directives
#BSUB -P CPH102
#BSUB -W 2:00
#BSUB -nnodes NODES
#BSUB -alloc_flags smt4
#BSUB -J dca_cooldown
#BSUB -o dca_cooldown.%J
#BSUB -e dca_cooldown.%J

# Change into scratch filesystem
cd $LS_SUBCWD

# Load required modules
source $MODULESHOME/init/bash

module reset
module load gcc/6.4.0
module load cuda
module load hdf5
module load fftw
module load cmake
module load magma
module load netlib-lapack
module load essl

date

JOBS

date
