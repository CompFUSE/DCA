
#export CXXFLAGS="-I/opt/cray/pe/papi/6.0.0.17/include"
#export LDFLAGS="-L/opt/cray/pe/papi/6.0.0.17/lib -lpapi"
#-DDCA_PROFILER=PAPI
export FFTW_PATH=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/gcc-12.2.0/fftw-3.3.10-tajdtzkealhold4bmpuq7wiwzurnclr4
export MAGMA_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/gcc-12.2.0/magma-2.7.2-gbjcrprqdw7y5uplm5upmqbi65zqwubb
export OPENBLAS_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/gcc-12.2.0/openblas-0.3.25-t62dxdtaqba6lzrwoy4uddswlprgma6n
export HDF5_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/gcc-12.2.0/hdf5-1.14.3-3so3g5x2roywum3edvjun7jbhwisei6p
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:/opt/rocm-5.7.0/hip/bin:${HDF5_ROOT}
export PATH=/sw/frontier/spack-envs/base/opt/linux-sles15-x86_64/gcc-7.5.0/cmake-3.23.2-4r4mpiba7cwdw2hlakh5i7tchi64s3qd/bin:${PATH}

cmake -DDCA_WITH_CUDA=off -DDCA_WITH_HIP=ON \
      -DFFTW_ROOT=$FFTW_PATH \
      -DDCA_FIX_BROKEN_MPICH=ON \
      -DROCM_ROOT=${ROCM_PATH} \
      -DMAGMA_ROOT=${MAGMA_ROOT} \
      -DLAPACK_ROOT=${OPENBLAS_ROOT} \
      -DBLAS_ROOT=${OPENBLAS_ROOT} \
      -DDCA_WITH_TESTS_FAST=ON \
      -DTEST_RUNNER="srun" \
      -DGPU_TARGETS=gfx90a \
      -DAMDGPU_TARGETS=gfx90a \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DCMAKE_HIP_COMPILER=/opt/rocm-5.7.0/llvm/bin/clang++ \
      -DCMAKE_INSTALL_PREFIX=$INST \
      -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
      -GNinja \
      ..
