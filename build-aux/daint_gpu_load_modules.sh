#echo "Clearing modules and loading DCA++ modules"
#module purge

module load daint-gpu
module load gcc/8.3.0
module load craype
module load craype-network-aries
module load cray-libsci
module load udreg
module load ugni
module load pmi
module load dmapp
module load gni-headers
module load xpmem
module load job
module load dvs
module load alps
module load rca
module load perftools-base
module load cudatoolkit
module load cray-mpich
module load cray-hdf5
module load cray-fftw
module unload cray-libsci/*
module load intel
module load cudatoolkit
module load cray-mpich
module load cray-hdf5
module load cray-fftw

export MAGMAROOT=/project/s299/dca_ext_libs
export PATH=/project/s299/easybuild/daint/haswell/software/cmake/cmake-3.18.2/bin/:${PATH}
