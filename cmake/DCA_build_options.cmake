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
  set(DCA_SOLVER_FILE "cluster_solver_mc_ctaux/ctaux_cluster_solver.h")
  add_definitions("-DUSE_CTAUX")
elseif (${DCA_QMC_TYPE} STREQUAL "SS-CT-HYB")
  set(DCA_CLUSTER_SOLVER_TYPE "DCA::SS_CT_HYB")
  set(DCA_SOLVER_FILE "cluster_solver_ss_hybridization/ss_hybridization_solver.h")
  add_definitions(-DUSE_SS_CT_HYB)
elseif (${DCA_QMC_TYPE} STREQUAL "HTS")
  set(DCA_CLUSTER_SOLVER_TYPE "DCA::HIGH_TEMPERATURE_SERIES")
  set(DCA_SOLVER_FILE "cluster_solver_series_expansion/high_temperature_series_expansion_solver.h")
  #TODO make self contained or remove.
  message(ERROR "High temperature solver is not fully supported")
  add_definitions(-DUSE_HTS)
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
  add_definitions(-DUSE_PTHREADS)
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
set(DCA_RNG "RAN" CACHE STRING  "Choose the random number generator, options are: SPRNG | RAN | RANQ2 | MT.")
if (${DCA_RNG} STREQUAL "RANQ2")
  add_definitions(-D CMAKE_RNG=rng::ranq2)
elseif(${DCA_RNG} STREQUAL "RAN")
  set(DCA_RNG_TYPE "rng::ran")
elseif(${DCA_RNG} STREQUAL "MT")
  set(DCA_RNG_TYPE "rng::rng_mt")
elseif(${DCA_RNG} STREQUAL "SPRNG")
  set(DCA_RNG_TYPE "rng::sprng_interface")
else()
  message(FATAL_ERROR "Please set RNG to a valid value: SPRNG | RAN | RANQ2 | MT.")
endif()

# Enable BIT.
option(DCA_QMC_INTEGRATOR_BIT "Enable BIT." OFF)
if (DCA_QMC_INTEGRATOR_BIT)
  add_definitions('-DDCA_QMC_INTEGRATOR_BIT')
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
