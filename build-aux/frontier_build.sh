
export CXXFLAGS="-I/opt/cray/pe/papi/6.0.0.17/include"
export LDFLAGS="-L/opt/cray/pe/papi/6.0.0.17/lib -lpapi"

cmake -DDCA_WITH_CUDA=off -DDCA_WITH_HIP=ON \
	  -DFFTW_INCLUDE_DIR=${OLCF_FFTW_ROOT}/include \
	  -DFFTW_LIBRARY=${OLCF_FFTW_ROOT}/lib/libfftw3.a \
	  -DDCA_FIX_BROKEN_MPICH=ON \
	  -DROCM_ROOT=${ROCM_PATH} \
	  -DDCA_PROFILER=PAPI \
	  -DDCA_WITH_TESTS_FAST=ON \
	  -DTEST_RUNNER="srun" \
	  -DGPU_TARGETS=gfx90a \
	  -DAMDGPU_TARGETS=gfx90a \
	  -DCMAKE_C_COMPILER=mpicc \
	  -DCMAKE_CXX_COMPILER=mpic++ \
	  -DCMAKE_INSTALL_PREFIX=$INST \
	  -DCMAKE_PREFIX_PATH=$INST \
 	  -GNinja \
	  ..
