nvhome=/opt/nvidia/hpc_sdk
target=Linux_ppc64le
version=20.11

nvcudadir=$nvhome/$target/$version/cuda/11.1
nvmathdir=$nvhome/$target/$version/math_libs/11.1
nvcommdir=$nvhome/$target/$version/comm_libs/11.1

export NVHPC=$nvhome
export PATH=/home/epd/local/hdf_cxx_12/bin:/home/epd/opt/ompi12/bin:/home/epd/opt/clang12/bin:/home/epd/opt/rtags/bin:/home/epd/local/bin:$nvcudadir/bin:/home/epd/opt/fftw3_clang12/bin:/home/epd/opt/cmake_3_19_2/bin:$PATH

export CPATH=/home/epd/opt/ompi12/include:/home/epd/opt/clang12/include:$nvcudadir/include:$nvmathdir/include:$nvcommdir/nccl/include:$nvcommdir/nvshmem/include:/home/epd/opt/fftw3_clang12/include:/home/epd/opt/magma_clang12_cuda11/include:$CPATH

export LD_LIBRARY_PATH=/home/epd/opt/clang12/lib:$nvcudadir/lib64:$nvmathdir/lib64:$nvcommdir/nccl/lib:$nvcommdir/nvshmem/lib:/home/epd/local/lib:/home/epd/local/hdf_cxx_12/lib:/home/epd/opt/magma_clang12_cuda11/lib:${LD_LIBRARY_PATH}

export MAN_PATH=/home/epd/opt/rtags/share/man:/home/epd/opt/ompi12/man:${MAN_PATH}

export CMAKE_PREFIX_PATH="/home/epd/opt/clang12:/home/epd/opt/cmake_3_19_2:${CMAKE_PREFIX_PATH}"

export HDF5_DIR="/home/epd/local/hdf_cxx_12"
export CUDA_DIR="${nvcudadir};${nvmathdir}"
export FFTW_DIR="/home/epd/opt/fftw3_clang12"


spack load ninja%clang@10.0.1
#spack load openblas
#spack load netlib-lapack
spack load libxml2/av5cl7l
spack load the-silver-searcher
spack load bash-completion

export CC=$(which mpicc)
export CXX=$(which mpic++)

#rm -rf *; CC=mpicc CXX=mpic++ CXXFLAGS=--gcc-toolchain=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-4.8.5/gcc-7.3.0-bco5a3lq3pzlot65mqywljwofqhsgxim CFLAGS=--gcc-toolchain=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-4.8.5/gcc-7.3.0-bco5a3lq3pzlot65mqywljwofqhsgxim cmake -DADIOS2_USE_MPI=ON ADIOS2_USE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DADIOS2_USE_Python=OFF -DADIOS2_USE_Blosc=OFF -DADIOS2_USE_Bzip2=OFF -DADIOS2_USE_SST=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_ZeroMQ=OFF -DADIOS2_USE_ZFP=OFF -DADIOS2_USE_SZ=OFF -DADIOS2_USE_MGARD=OFF -GNinja -DCMAKE_CXX_STANDARD=17 ..
#rm -rf *; export ADIOS2_DIR=/home/epd/ADIOS2/build_clang17; export HDF5_ROOT=${HDF5_DIR}; CXXFLAGS=--gcc-toolchain=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-4.8.5/gcc-7.3.0-bco5a3lq3pzlot65mqywljwofqhsgxim CFLAGS=--gcc-toolchain=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-4.8.5/gcc-7.3.0-bco5a3lq3pzlot65mqywljwofqhsgxim cmake -DDCA_ESSL_INCLUDES=/opt/ibmmath/essl/5.5/include -DDCA_HAVE_LAPACK=TRUE -DLAPACK_LIBRARIES="/opt/ibmmath/essl/5.5/lib64/libessl.so;/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-7.3.0/netlib-lapack-3.8.0-5vb67z6mxivntovciglsssz54nundlxv/lib64/liblapack.so;/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-7.3.0/netlib-xblas-1.0.248-pd7d4wruu3efdaupjijl57sy37hped73/lib/libblas.a" -DTEST_RUNNER="mpirun" -DMPIEXEC_NUMPROC_FLAG="-n" -DMPIEXEC_PRE_FLAGS="-mca btl self,tcp,smcuda,vader" -DDCA_WITH_CUDA=TRUE -DCUDA_GPU_ARCH=sm_70 -DFFTW_INCLUDE_DIR=/opt/ibmmath/essl/5.5/FFTW3/include -DFFTW_LIBRARY="/usr/local/lib64/libfftw3_essl_gcc.a;/opt/ibm/xlf/16.1.0/lib/libxlf90_r.a;/opt/ibm/xlf/16.1.0/lib/libxlfmath.a;/opt/ibmmath/essl/5.5/lib64/libessl.so" -DMAGMA_DIR=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-7.3.0/magma-2.4.0-s4fji5q5t5esgsvfsk7xm -DHDF5_ROOT=${HDF5_ROOT} -DHDF5_INCLUDE_DIRS=${HDF5_ROOT}/include -DHDF5_LIBRARIES="${HDF5_ROOT}/lib/libhdf5.so;${HDF5_ROOT}/lib/libhdf5_cpp.so" -DDCA_WITH_TESTS_FAST=True -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_PREFIX_PATH="${ADIOS2_DIR}" -GNinja -DCMAKE_CXX_STANDARD=17 ..
#rm -rf *; export ADIOS2_DIR=/home/epd/ADIOS2/build_clang17; export HDF5_ROOT=${HDF5_DIR}; CXXFLAGS=--gcc-toolchain=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-4.8.5/gcc-7.3.0-bco5a3lq3pzlot65mqywljwofqhsgxim CFLAGS=--gcc-toolchain=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-4.8.5/gcc-7.3.0-bco5a3lq3pzlot65mqywljwofqhsgxim cmake -DDCA_ESSL_INCLUDES=/opt/ibmmath/essl/5.5/include -DDCA_HAVE_LAPACK=TRUE -DLAPACK_LIBRARIES="/opt/ibmmath/essl/5.5/lib64/libessl.so;/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-7.3.0/netlib-lapack-3.8.0-5vb67z6mxivntovciglsssz54nundlxv/lib64/liblapack.so;/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-7.3.0/netlib-xblas-1.0.248-pd7d4wruu3efdaupjijl57sy37hped73/lib/libblas.a" -DTEST_RUNNER="mpirun" -DMPIEXEC_NUMPROC_FLAG="-n" -DMPIEXEC_PRE_FLAGS="-mca btl self,tcp,smcuda,vader" -DDCA_WITH_CUDA=TRUE -DCUDA_GPU_ARCH=sm_70 -DFFTW_INCLUDE_DIR=/opt/ibmmath/essl/5.5/FFTW3/include -DFFTW_LIBRARY="/usr/local/lib64/libfftw3_essl_gcc.a;/opt/ibm/xlf/16.1.0/lib/libxlf90_r.a;/opt/ibm/xlf/16.1.0/lib/libxlfmath.a;/opt/ibmmath/essl/5.5/lib64/libessl.so" -DMAGMA_DIR=/home/epd/spack/opt/spack/linux-centos7-ppc64le/gcc-7.3.0/magma-2.4.0-s4fji5q5t5esgsvfsk7xm -DHDF5_ROOT=${HDF5_ROOT} -DHDF5_INCLUDE_DIRS=${HDF5_ROOT}/include -DHDF5_LIBRARIES="${HDF5_ROOT}/lib/libhdf5.so;${HDF5_ROOT}/lib/libhdf5_cpp.so" -DDCA_WITH_TESTS_FAST=True -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_PREFIX_PATH="${ADIOS2_DIR}" -GNinja -DCMAKE_CXX_STANDARD=17 -DDCA_WITH_ADIOS2=TRUE -DCMAKE_CXX_FLAGS_DEBUG="-fstandalone-debug" -DCMAKE_BUILD_TYPE=Debug ..
