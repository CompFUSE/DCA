# Initial cache list for CADES-CONDO
#
# Usage: cmake -C /path/to/this/file /path/to/DCA/source -D<option>=<value> ...

add_compile_options(-fPIC -march=native -mtune=native)

set(FFTW_DIR /software/user_tools/centos-7.2.1511/cades-cnms/spack/opt/spack/linux-centos7-x86_64/clang-7.9.0/fftw-3.3.8-3mx5kwmdcznqlizyos5uipzhpj6jou2w)
set(FFTW_INCLUDE_DIR ${FFTW_DIR}/include CACHE PATH "FFTW include path")
set(FFTW_LIBRARY "${FFTW_DIR}/lib/libfftw3.a;${FFTW_DIR}/lib/libfftw3f.a" CACHE FILEPATH "FFTW libary")

set(TEST_RUNNER "mpirun" CACHE STRING "Command for executing (MPI) programs.")
set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING
  "Flag used by TEST_RUNNER to specify the number of processes.")
set(MPIEXEC_PREFLAGS "" CACHE STRING
  "Flags to pass to TEST_RUNNER directly before the executable to run.")

# Enable the GPU support.
option(DCA_WITH_CUDA "Enable GPU support." ON)
option(DCA_WITH_MPI "Enable MPI support." ON)
# Compile for Volta compute architecture.
set(CUDA_GPU_ARCH "sm_60" CACHE STRING "Name of the *real* architecture to build for.")

# For the GPU support we also need MAGMA.
#set(MAGMA_DIR $ENV{MAGMA_DIR} CACHE PATH)
#  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_DIR})

option(DCA_WITH_TESTS_FAST "Enable fast tests" ON)

# Followed by for "llvm7.9.0" i.e. dev llvm with its version set back for nvcc
#rm -rf *; CXX=$(which mpic++) CC=$(which mpicc) CXXFLAGS="--gcc-toolchain=/software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5" LDFLAGS="-L/software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5/lib64 -Wl,-rpath,/software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5" LD_FLAGS="-L/software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5/lib64 -Wl,-rpath,/software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5 -fuse-ld=lld" cmake -C../build-aux/cades-condo.cmake -DCMAKE_EXE_LINKER_FLAGS="-L/software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5/lib64 -Wl,-rpath,/software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5 -fuse-ld=lld" -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -GNinja ..
