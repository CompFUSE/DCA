module purge
module load cmake
module load ARM_Compiler_For_HPC/20.3_TX2 
module load openmpi/4.0.5_armclang
module load hdf5-1.10.6-arm-20.3-a7xz2pp
module load fftw-3.3.8-arm-20.3-qi3atrb
module load openblas-0.3.8-arm-20.3-jcthkv3
module load netlib-lapack-3.8.0-arm-20.3-bbtrjci

export CC=mpicc
export CXX=mpicxx

# cuda runtime requires when cuda = 10, gcc < 9 and clang < 9
# only wombat6 has temporarily installed cuda 10 
# cuda 11 modules on Wombat might be buggy
# gcc modules on Wombat all > 9 
# it leaves us only one combination, use cuda 10 and armclang on Wombat6
export nvcc=/usr/local/cuda-10.2/bin/nvcc

# module load cuda-11.2.1-arm-20.3-6ibfl33
#. /autofs/nccs-svm1_envoy_od/weile/dev/src/spack/share/spack/setup-env.sh
