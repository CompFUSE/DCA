#!/bin/bash
#
# Loads all modules that are required to build DCA++ on  ORNL's CADES CONDO cluster
#
# Usage: source cades-condo_load_modules.sh

# uncomment if your default module loadout clashes with this
# module purge

module load env/cades-cnms
. $SOFTWARECNMS/spack/share/spack/setup-env.sh
spack load gcc/egooyqw
spack load openmpi@3.1.3%gcc@6.5.0
spack load ninja/gzwd46m
spack load emacs/qp7x25b
spack load hdf5/4gmsnjn
spack load cmake/q76ndqk
spack load git/kibjjo6
spack load fftw/kpdartc
module load cuda/9.2
spack load magma/y3qdgaz
export FFTW_DIR=/software/user_tools/centos-7.2.1511/cades-cnms/spack/opt/spack/linux-centos7-x86_64/gcc-6.5.0/fftw-3.3.8-kpdartcqxfk2kdsbcfdtwin75s24z5uu
