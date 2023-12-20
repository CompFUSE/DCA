
#export CXXFLAGS="-I/opt/cray/pe/papi/6.0.0.17/include"
#export LDFLAGS="-L/opt/cray/pe/papi/6.0.0.17/lib -lpapi"
#-DDCA_PROFILER=PAPI

FFTW_PATH=/sw/frontier/spack-envs/base/opt/cray-sles15-zen3/gcc-11.2.0/fftw-3.3.10-njtwalw5xctv6e3fshucadlgu32jb4k6
MAGMA_ROOT=/lustre/orion/world-shared/cph102/epd/gcc-11.2.0/magma-2.7.2-c5m5kzaz7irix5hk5zzf3mrwwlij43is
OPENBLAS_ROOT=/lustre/orion/world-shared/cph102/epd/gcc-11.2.0/openblas-0.3.25-scaywvuh5zsm5u7smg54plj2oyf7nekv
HDF5_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/gcc-11.2.0/hdf5-1.14.3-rif5452vfxb2jzpa7nywj7rokpaofyde

cmake -DDCA_WITH_CUDA=off -DDCA_WITH_HIP=ON \
      -DFFTW_ROOT=$FFTW_PATH \
      -DDCA_FIX_BROKEN_MPICH=ON \
      -DROCM_ROOT=${OLCF_ROCM_ROOT} \
      -DMAGMA_ROOT=${MAGMA_ROOT} \
      -DLAPACK_ROOT=${OPENBLAS_ROOT} \
      -DBLAS_ROOT=${OPENBLAS_ROOT} \
      -DDCA_WITH_TESTS_FAST=ON \
      -DTEST_RUNNER="srun" \
      -DGPU_TARGETS=gfx90a \
      -DAMDGPU_TARGETS=gfx90a \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DCMAKE_INSTALL_PREFIX=$INST \
      -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH \
      -DHDF5_ROOT=${HDF5_ROOT} \
      -GNinja \
      ..
