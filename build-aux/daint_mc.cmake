# Initial cache list for Piz Daint multicore (Cray XC40)
#
# Usage: cmake /path/to/dca/source -C /path/to/this/file -D<option>=<value> ...

# Prevent CMake from searching for BLAS and LAPACK libraries.
# CC automatically links against them.
set(DCA_HAVE_LAPACK TRUE CACHE INTERNAL "")

# Use srun for executing the tests.
set(TEST_RUNNER "srun" CACHE STRING "Command for executing (MPI) programs.")
