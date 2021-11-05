
cmake -DDCA_WITH_CUDA=OFF -DDCA_WITH_HIP=ON \
      -DFFTW_INCLUDE_DIR=${OLCF_FFTW_ROOT}/include \
      -DFFTW_LIBRARY=${OLCF_FFTW_ROOT}/lib/libfftw3.a \
      -GNinja \
      -DROCM_ROOT=${ROCM_PATH} \
      -DDCA_WITH_TESTS_FAST=ON -DTEST_RUNNER="srun" \
      -DGPU_TARGETS=gfx906 \
      -DCMAKE_C_COMPILER=cc \
      -DCMAKE_CXX_COMPILER=CC \
      ..
