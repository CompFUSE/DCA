# Initial cache list for Wombat
#
# Usage: cmake -C /path/to/this/file /path/to/DCA/source -D<option>=<value> ...

# Disable the GPU support.
option(DCA_WITH_CUDA "Enable GPU support." ON)
option(DCA_WITH_MPI "Enable MPI support." ON)

set(MAGMA_DIR /ccsopen/home/weile/dev/src/magma-2.5.4/ CACHE PATH
  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

