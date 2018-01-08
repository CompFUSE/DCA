#!/bin/bash
#
#PBS -A cph102
#PBS -N DCA++_APPLICATION_U=HUBBARDU_d=DENS_Nc=SIZE
#PBS -o out.APPLICATION_U=HUBBARDU_d=DENS_Nc=SIZE.$PBS_JOBID.txt
#PBS -e err.APPLICATION_U=HUBBARDU_d=DENS_Nc=SIZE.$PBS_JOBID.txt
#PBS -l nodes=1
#PBS -l walltime=01:00:00
#PBS -l gres=atlas1%atlas2

export OMP_NUM_THREADS=8
export CRAY_CUDA_MPS=1  # Required for running multiple tasks per node/gpu.

cd $PBS_O_WORKDIR

date

JOBS

date
