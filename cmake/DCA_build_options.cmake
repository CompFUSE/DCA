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
option(DCA_WITH_PTHREADS "Enable pthreads" ON)
if (DCA_WITH_PTHREADS)
  dca_add_config_define(DCA_HAVE_PTHREADS)
  set(DCA_MC_INTEGRATOR_TYPE "DCA::posix_qmci_integrator<ClusterSolverBaseType>")
  dca_add_config_define(DCA_MC_INTEGRATOR_TYPE "DCA::posix_qmci_integrator<ClusterSolverBaseType>")
  dca_add_config_define(DCA_PARALLELIZATION POSIX_LIBRARY)
else()
  set(DCA_MC_INTEGRATOR_TYPE "ClusterSolverBaseType")
endif()

# Enable profiling and auto-tuning.
set(DCA_PROFILING_MODE "NONE" CACHE STRING "Choose the profiling mode, options are: NONE|COUNTING|PAPI.")
option(DCA_AUTO_TUNING "Enable auto-tuning." OFF)
if (${DCA_PROFILING_MODE} STREQUAL "COUNTING")
  add_definitions(-DDCA_COUNTING_PROFILING)
  if(DCA_AUTO_TUNING)
    add_definitions(-DAUTOTUNING_ENABLED)
  endif()
elseif (${DCA_PROFILING_MODE} STREQUAL "PAPI")
  add_definitions(-DDCA_PAPI_PROFILING)
  if (DCA_AUTO_TUNING)
    add_definitions(-DAUTOTUNING_ENABLED)
  endif()
else()
  add_definitions(-DDCA_NO_PROFILING)
endif()

################################################################################
# Select the random number generator
set(DCA_RNG "std::ranlux48_base" CACHE STRING "Random number generator, options are: std::ranlux48_base | std::ranlux48 | std::mt19937_64 | SPRNG_LFG | SPRNG_MLFG | Ranq2.")
set_property(CACHE DCA_RNG
  PROPERTY STRINGS std::ranlux48_base std::ranlux48 std::mt19937_64 SPRNG_LFG SPRNG_MLFG Ranq2)

if(${DCA_RNG} STREQUAL "std::ranlux48_base")
  set(DCA_RNG dca::math::random::StdRandomWrapper<std::ranlux48_base>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")

elseif(${DCA_RNG} STREQUAL "std::ranlux48")
  set(DCA_RNG dca::math::random::StdRandomWrapper<std::ranlux48>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")

elseif(${DCA_RNG} STREQUAL "std::mt19937_64")
  set(DCA_RNG dca::math::random::StdRandomWrapper<std::mt19937_64>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")

elseif(${DCA_RNG} STREQUAL "SPRNG_LFG")
  set(DCA_RNG dca::math::random::SprngWrapper<dca::math::random::LFG>)
  set(DCA_RNG_INCLUDE "dca/math/random/sprng_wrapper.hpp")

elseif(${DCA_RNG} STREQUAL "SPRNG_MLFG")
  set(DCA_RNG dca::math::random::SprngWrapper<dca::math::random::MLFG>)
  set(DCA_RNG_INCLUDE "dca/math/random/sprng_wrapper.hpp")

elseif(${DCA_RNG} STREQUAL "Ranq2")
  set(DCA_RNG dca::math::random::Ranq2)
  set(DCA_RNG_INCLUDE "dca/math/random/ranq2.hpp")
  
else()
  message(FATAL_ERROR "Please set DCA_RNG to a valid option: std::ranlux48_base | std::ranlux48 |
                       std::mt19937_64 | SPRNG_LFG | SPRNG_MLFG | Ranq2.")
endif()

configure_file("${PROJECT_SOURCE_DIR}/cmake/templates/config_rng.hpp.in"
  "${PROJECT_BINARY_DIR}/include/dca/config/rng.hpp"
  @ONLY)

################################################################################
# Enable BIT.
option(DCA_QMC_INTEGRATOR_BIT "Enable BIT." OFF)
if (DCA_QMC_INTEGRATOR_BIT)
  dca_add_config_define(DCA_QMC_INTEGRATOR_BIT)
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

# gnuplot
option(DCA_WITH_GNUPLOT "Enable gnuplot." OFF)
if (DCA_WITH_GNUPLOT)
  find_package(Gnuplot REQUIRED)
  dca_add_config_define(DCA_HAVE_GNUPLOT)
endif()
