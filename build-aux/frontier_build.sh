cmake -DDCA_WITH_CUDA=off -DDCA_WITH_HIP=ON \
      -DFFTW_ROOT=$FFTW_PATH \
      -DDCA_WITH_ADIOS2=OFF \
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
      -DCMAKE_HIP_COMPILER=/opt/rocm-6.3.1/llvm/bin/clang++ \
      -DCMAKE_INSTALL_PREFIX=$INST \
      -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
      -DCMAKE_HIP_LINK_FLAGS=--hip-link \
      -GNinja \
      ..
