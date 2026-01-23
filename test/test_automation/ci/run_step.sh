#!/bin/bash

set -x
HOST_NAME=$(hostname -s)

case "$1" in

  # Configure DCA++ using cmake out-of-source builds
  configure)

    if [ -d ${GITHUB_WORKSPACE}/../dca-build ]
    then
      echo "Found existing out-of-source build directory ${GITHUB_WORKSPACE}/../dca-build, removing"
      rm -fr ${GITHUB_WORKSPACE}/../dca-build
    fi

    echo "Creating new out-of-source build directory ${GITHUB_WORKSPACE}/../dca-build"
    cd ${GITHUB_WORKSPACE}/..
    mkdir dca-build
    cd dca-build

    # Real or Complex configuration
    case "${GH_JOBNAME}" in
      *"Real"*)
        echo 'Configure for real build'
        IS_COMPLEX=0
      ;;
      *"Complex"*)
        echo 'Configure for complex build'
        IS_COMPLEX=1
      ;;
    esac

    if [[ "${GH_JOBNAME}" =~ (-CUDA) ]]
    then
      echo "Set PATH to cuda-12.9 to be associated with the C and C++ compilers"
        export PATH=/usr/local/cuda-12.9/bin:$PATH
        echo "Set CUDACXX CMake environment variable to nvcc cuda 11.8"
        export CUDACXX=/usr/local/cuda-12.9/bin/nvcc
        # Make current environment variables available to subsequent steps
        echo "PATH=$PATH" >> $GITHUB_ENV
        echo "CUDACXX=$CUDACXX" >> $GITHUB_ENV
    fi

    # Mixed or Non-Mixed (default, full) precision, used with GPU code
    case "${GH_JOBNAME}" in
      *"Mixed"*)
        echo 'Configure for mixed precision build'
        IS_MIXED_PRECISION=1
      ;;
      *)
        IS_MIXED_PRECISION=0
      ;;
    esac

    # Path to QMC_DATA in self-hosted CI system and point at minimum gcc-9
    if [[ "$HOST_NAME" =~ (a30four) ]]
    then
	# use gcc-12
      export PATH=$(spack find -lp llvm | awk '/llvm/{print $3}')/bin:${PATH}
      export LD_LIBRARY_PATH=$(spack find -lp gcc/yel5m4n | awk '/gcc/{print $3}')/lib64:$(spack find -lp hdf5 | awk '/hdf5/{print $3}')/lib64:$(spack find -lp llvm/navafcm | awk '/llvm/{print $3}')/lib:${LD_LIBRARY_PATH}:$(spack find -lp tree-sitter-cmake | awk '/tree-sitter-cmake/{print $3}')/lib:${LD_LIBRARY_PATH}
      export CUDAHOSTCXX=clang++
      export MAGMA_ROOT=/home/epd/opt_a30/magma
      export HDF5_ROOT=$(spack find --loaded -lp hdf5 | awk '/hdf5/{print $3}')
      export OPENBLAS_ROOT=$(spack find --loaded -lp openblas | awk '/openblas/{print $3}')
      export MPI_ROOT=$(spack find --loaded -lp openmpi | awk '/openmpi/{print $3}')
      export FFTW_ROOT=$(spack find --loaded -lp fftw | awk '/fftw/{print $3}')
      # Make current environment variables available to subsequent steps
      echo "PATH=$PATH" >> $GITHUB_ENV
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $GITHUB_ENV
    fi

    case "${GH_JOBNAME}" in
      *"LLVM21-MPI-CUDA-"*)
        echo 'Configure for debug mode to capture asserts with gcc'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=mpicc \
              -DCMAKE_CXX_COMPILER=mpic++ \
              -DDCA_WITH_CUDA=1 \
	      -DCMAKE_CUDA_FLAGS=-allow-unsupported-compiler \
	      -DCMAKE_CUDA_ARCH_LIST=80 \
	      -DCMAKE_FORTRAN_COMPILER=gfortran -DCMAKE_EXE_LINKER_FLAGS_INIT="-lgfortran" -DCMAKE_SHARED_LINKER_FLAGS_INIT="-lgfortran" -DCMAKE_MODULE_LINKER_FLAGS_INIT="-lgfortran" \
              -DDCA_WITH_MPI=1 \
	      -DCMAKE_PREFIX_PATH=${MPI_ROOT}\;${CUDA_ROOT}\;${MAGMA_ROOT}\;${HDF5_ROOT}\;${OPENBLAS_ROOT}\;${ADIOS2_ROOT}\;${FFTW_ROOT} \
              -DCMAKE_BUILD_TYPE=Release \
	      -DTEST_RUNNER="mpiexec" \
	      -DMPIEXEC_NUMPROC_FLAG="-n" -DMPIEXEC_PREFLAGS="-mca btl self,tcp" -DDCA_WITH_CUDA=1 -DDCA_WITH_ADIOS2=1 \
	      -DDCA_WITH_TESTS_FAST=1 \
              ${GITHUB_WORKSPACE}
      ;;
    esac
    ;;

  # Build using ninja (~ 25 minutes on GitHub-hosted runner)
  build)
    # CUDA toolchain can be used implicitly by the compiler. Double check the location.
    if [[ "${GH_JOBNAME}" =~ (CUDA) ]]
    then
      which nvcc
    fi

    cd ${GITHUB_WORKSPACE}/../dca-build
    ninja
    ;;

  # Run deterministic tests
  test)

    # Run only deterministic tests (reasonable for CI) by default
    TEST_LABEL=""
    export MPI_ROOT=$(spack find --loaded -lp openmpi | awk '/openmpi/{print $3}')
    export PATH=${MPI_ROOT}/bin:${PATH}
    cd ${GITHUB_WORKSPACE}/../dca-build

    # Add ctest concurrent parallel jobs
    # Default for Linux GitHub Action runners
    CTEST_JOBS="1"

    ctest --output-on-failure $TEST_LABEL -j $CTEST_JOBS
    ;;

  *)
    echo " Invalid step" "$1"
    exit -1
    ;;
esac
