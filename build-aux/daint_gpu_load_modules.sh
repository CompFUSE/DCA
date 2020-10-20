#echo "Clearing modules and loading DCA++ modules"
#module purge

module load daint-gpu
module swich PrgEnv-cray PrgEnv-gnu

module load cudatoolkit
module load cray-mpich
module load cray-hdf5
module load cray-fftw

# Use intel linear algebra libraries
module unload cray-libsci/*
module load intel

export MAGMAROOT=/project/s299/easybuild/daint/haswell/software/magma/2.5.3-gcc-8.3-cuda-10.2
export PATH=/project/s299/easybuild/daint/haswell/software/cmake/cmake-3.18.2/bin/:${PATH}
