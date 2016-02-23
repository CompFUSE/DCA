################################################################################
# Build options for DCA++
#
# The options result in
# - preprocessor definitions (-D) and
# - direct replacements in the source code.
################################################################################

# Choose the Monte-Carlo algorithm.
set(DCA_QMC_TYPE "CT-AUX" CACHE STRING "Choose the Monte-Carlo algorithm, options are: CT-AUX|SS-CT-HYB|HTS.")
if (${DCA_QMC_TYPE} STREQUAL "CT-AUX")
  set(DCA_CLUSTER_SOLVER_TYPE "DCA::CT_AUX_CLUSTER_SOLVER")
elseif (${DCA_QMC_TYPE} STREQUAL "SS-CT-HYB")
  set(DCA_CLUSTER_SOLVER_TYPE "DCA::SS_CT_HYB")
elseif (${DCA_QMC_TYPE} STREQUAL "HTS")
  set(DCA_CLUSTER_SOLVER_TYPE "DCA::HIGH_TEMPERATURE_SERIES")
else()
  message(FATAL_ERROR "Please set QMC_TYPE to a valid option: CT-AUX|SS-CT-HYB|HTS.")
endif()

# Enable SSE acceleration.
option(DCA_SSE_ACCELERATION "Enable SSE acceleration." OFF)
mark_as_advanced(DCA_SSE_ACCELERATION)
if (DCA_SSE_ACCELERATION)
  add_definitions(-DUSE_SSE_ACCELERATION)
endif()

# Enable pthreads.
option(DCA_PTHREADS "Enable pthreads" ON)
if (DCA_PTHREADS)
  set(DCA_MC_INTEGRATOR_TYPE "posix_qmci_integrator<quantum_cluster_solver_type>")
else()
  set(DCA_MC_INTEGRATOR_TYPE "quantum_cluster_solver_type")
endif()

# Enable profiling and auto-tuning.
set(DCA_PROFILING_MODE "NONE" CACHE STRING "Choose the profiling mode, options are: NONE|COUNTING|PAPI.")
option(DCA_AUTO_TUNING "Enable auto-tuning." OFF)
if (${DCA_PROFILING_MODE} STREQUAL "COUNTING")
  add_definitions(-DCOUNTING_PROFILING)
  if(DCA_AUTO_TUNING)
    add_definitions(-DAUTOTUNING_ENABLED)
  endif()
elseif (${DCA_PROFILING_MODE} STREQUAL "PAPI")
  add_definitions(-DDCA_PAPI_PROFILING)
  if (DCA_AUTO_TUNING)
    add_definitions(-DAUTOTUNING_ENABLED)
  endif()
else()
  add_definitions(-DNO_PROFILING)
endif()

# Choose the random number generator.
set(DCA_RNG "RANQ2" CACHE STRING "Choose the random number generator, options are: RAN | RANQ2 | MT. The first two are from Numerical Recipes, the last from the standard Mersenne Twistor.")
if (${DCA_RNG} STREQUAL "RANQ2")
  set(RNG_TYPE "rng::ranq2")
elseif(${DCA_RNG} STREQUAL "RAN")
  set(RNG_TYPE "rng::ran")
elseif(${DCA_RNG} STREQUAL "MT")
  set(RNG_TYPE "rng::rng_mt")
else()
  message(FATAL_ERROR "Please set RNG to a valid value: NR|SPRNG.")
endif()

# Enable BIT.
option(DCA_QMC_INTEGRATOR_BIT "Enable BIT." OFF)
if (DCA_QMC_INTEGRATOR_BIT)
  set(DCA_QMC_INTEGRATOR_BIT_VALUE "true")
else()
  set(DCA_QMC_INTEGRATOR_BIT_VALUE "false")
endif()

# Compute the standard deviation.
option(DCA_STANDARD_DEVIATION "Compute the standard deviation." OFF)
if (DCA_STANDARD_DEVIATION)
  add_definitions(-DMEASURE_ERROR_BARS)
endif()

# Single precision measurement.
option(DCA_SINGLE_PRECISION_MEASUREMENTS "Measure in single precision." OFF)
mark_as_advanced(DCA_SINGLE_PRECISION_MEASUREMENTS)
if (DCA_SINGLE_PRECISION_MEASUREMENTS)
  add_definitions(-DSINGLE_PRECISION_MEASUREMENTS)
endif()

# Single precision coarsegraining.
option(DCA_SINGLE_PRECISION_COARSEGRAINING "Coarsegrain in single precision." OFF)
mark_as_advanced(DCA_SINGLE_PRECISION_COARSEGRAINING)
if (DCA_SINGLE_PRECISION_COARSEGRAINING)
  add_definitions(-DSINGLE_PRECISION_COARSEGRAINING)
endif()

# Use full vertex.
option(DCA_FULL_VERTEX "Use full vertex." OFF)
mark_as_advanced(DCA_FULL_VERTEX)
if (NOT DCA_FULL_VERTEX)
  add_definitions(-DUSE_REDUCED_VERTEX_FUNCTION)
endif()

# Allow Gnuplot.
option(DCA_GNUPLOT "Allow gnuplot." OFF)
if (DCA_GNUPLOT)
  find_package(Gnuplot)
  if (DCA_GNUPLOT_FOUND)
    message(STATUS "Found Gnuplot.")
    add_definitions(-DALLOW_GNUPLOT)
  else()
    message(WARNING "Could NOT find Gnuplot. Disable Gnuplot.")
    set(DCA_GNUPLOT OFF CACHE BOOL "Allow gnuplot." FORCE)
  endif()
endif()
