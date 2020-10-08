# Spack modules to get DCA build on Cades

# If you aren't using the suggested CNMS environment you need to uncomment the following two lines.

# module load env/cades-cnms
#. $SOFTWARECNMS/spack/share/spack/setup-env.sh

module load PE-gnu/3.0
spack load emacs@26.3
spack load git
spack load gcc@8.2.0
spack load openmpi/qnfab5m
spack load fftw%gcc@8.2.0
spack load ninja/v2bqky4
spack load cmake/g4ybxxf
spack load openblas@0.3.9

export HDF5_DIR=/software/user_tools/current/cades-cnms/for_nti/hdf5
export MAGMA_DIR=/lustre/or-hydra/cades-cnms/epd/dev/magma
export CUDA_DIR=/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/cuda/11.0/centos7.8_binary
export CUDADIR=/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/cuda/11.0/centos7.8_binary
export CMAKE_PREFIX_PATH=${HDF5_DIR}:${MAGMA_DIR}:$CMAKE_PREFIX_PATH

export FFTW_DIR=`spack find --loaded -p fftw | awk -e '/fftw/ {print $2}'`

export CC=$(which mpicc)
export CXX=$(which mpic++)

# cmake like this can work if you don't want to use cades.cmake
#rm -rf *; CXX=$(which mpic++) CC=$(which mpicc) cmake -DCUDA_TOOLKIT_ROOT_DIR=${CUDA_DIR} -DMAGMA_DIR=${MAGMA_DIR} -DDCA_WITH_CUDA=True -DCUDA_GPU_ARCH=sm_60 -DHDF5_ROOT=${HDF5_DIR} -DHDF5_INCLUDE_DIRS=${HDF5_DIR}/include -DHDF5_LIBRARIES="${HDF5_DIR}/lib/libhdf5_cpp.a;${HDF5_DIR}/lib/libhdf5.a" -DDCA_WITH_TESTS_FAST=True -DTEST_RUNNER=srun -DMPIEXEC_NUMPROC_FLAG="-n" -DMPIEXEC_PREFLAGS="--mem=64G --gpus-per-task=1" -DCMAKE_EXE_LINKER_FLAGS="-ldl -fopenmp" -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DFFTW_DIR=${FFTW_DIR} -DFFTW_INCLUDE_DIR=${FFTW_DIR}/include -DFFTW_LIBRARY="${FFTW_DIR}/lib/libfftw3.a;${FFTW_DIR}/lib/libfftw3f.a" -GNinja ..
