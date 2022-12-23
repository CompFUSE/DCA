#!/bin/bash

set -x
HOST_NAME=$(hostname -s)

case "$1" in 

  # Configure qmcpack using cmake out-of-source builds 
  configure)
    
    if [ -d ${GITHUB_WORKSPACE}/../dca-build ]
    then
      echo "Found existing out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build, removing"
      rm -fr ${GITHUB_WORKSPACE}/../dca-build
    fi
    
    echo "Creating new out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build"
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
      echo "Set PATH to cuda-11.8 to be associated with the C and C++ compilers"
        export PATH=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/cuda-11.8.0-3fgvnzslqxogwa5djfpvthlj2ebs3w4f/bin:$PATH
        echo "Set CUDACXX CMake environment variable to nvcc cuda 11.8"
        export CUDACXX=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/cuda-11.8.0-3fgvnzslqxogwa5djfpvthlj2ebs3w4f/bin/nvcc
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
    if [[ "$HOST_NAME" =~ (volta-cidev-node) ]]
    then
      # use gcc-12
      export PATH=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/gcc-12.2.0-cx7pjxgmemcce4tohlmsekuo5qvgjqbl/bin:/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/ninja-1.11.1-plzpokehn3kdbcviteppqntkqun5752f/bin:/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/cmake-3.25.0-xlxorwhfz5jxpyx65ypsh2horyo7n3ef/bin:/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/binutils-2.38-ljwszfsihiioxxusqiwg76hsr3t7eazd/bin:/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/cuda-11.8.0-3fgvnzslqxogwa5djfpvthlj2ebs3w4f/bin:$PATH
      export LD_LIBRARY_PATH=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/gcc-12.2.0-cx7pjxgmemcce4tohlmsekuo5qvgjqbl/lib64:$LD_LIBRARY_PATH
      export CUDA_ROOT=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/cuda-11.8.0-3fgvnzslqxogwa5djfpvthlj2ebs3w4f
      export MAGMA_ROOT=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/magma-2.7.0-hzsyoya4k46b7ufbtccl7bh3yxgubcxd
      export HDF5_ROOT=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/hdf5-1.12.2-lisnrzzj7stud3m3okbsotk7eume7665
      export OPENBLAS_ROOT=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/openblas-0.3.21-hxlxmz56qgkexpvfku2dugwacjf2auvl
      export ADIOS2_ROOT=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/adios2-2.8.3-ln3ke2w6zishdslde526mvwjuurr4ayh
      export MPI_ROOT=/home/epd/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-12.2.0/openmpi-2.1.6-24emdbitdevs7ecwquhtc6433ausfff7
      # Make current environment variables available to subsequent steps
      echo "PATH=$PATH" >> $GITHUB_ENV
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $GITHUB_ENV
    fi
    
    case "${GH_JOBNAME}" in
      *"GCC12-MPI-CUDA-"*)
        echo 'Configure for debug mode to capture asserts with gcc'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
	      -DCMAKE_FORTRAN_COMPILER=gfortran -DCMAKE_EXE_LINKER_FLAGS_INIT="-lgfortran" -DCMAKE_SHARED_LINKER_FLAGS_INIT="-lgfortran" -DCMAKE_MODULE_LINKER_FLAGS_INIT="-lgfortran" \
              -DDCA_WITH_MPI=1 \
	      -DCMAKE_PREFIX_PATH=${MPI_ROOT}\;${CUDA_ROOT}\;${MAGMA_ROOT}\;${HDF5_ROOT}\;${OPENBLAS_ROOT}\;${ADIOS2_ROOT} \
              -DCMAKE_BUILD_TYPE=Release \
	      -DTEST_RUNNER="mpiexec" \
	      -DMPIEXEC_NUMPROC_FLAG="-n" -DMPIEXEC_PREFLAGS="-mca btl self,tcp" -DDCA_WITH_CUDA=1 -DDCA_WITH_ADIOS2=1 \
	      -DDCA_WITH_TESTS_FAST=1 \
	      -DLAPACK_LIBRARIES=${OPENBLAS_ROOT}/lib64/libopenblas.a \
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
    
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    
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
