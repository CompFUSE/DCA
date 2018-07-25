################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Build options for DCA++.

################################################################################
# Enable MPI.
option(DCA_WITH_MPI "Enable MPI." ON)

if (DCA_WITH_MPI)
  include(dca_mpi)
  if (NOT DCA_HAVE_MPI)
    message(FATAL_ERROR "MPI not found but requested.")
  endif()

  set(DCA_CONCURRENCY_TYPE dca::parallel::MPIConcurrency)
  set(DCA_CONCURRENCY_INCLUDE "dca/parallel/mpi_concurrency/mpi_concurrency.hpp")
  set(DCA_CONCURRENCY_LIB parallel_mpi_concurrency)

else()
  set(DCA_CONCURRENCY_TYPE dca::parallel::NoConcurrency)
  set(DCA_CONCURRENCY_INCLUDE "dca/parallel/no_concurrency/no_concurrency.hpp")
  set(DCA_CONCURRENCY_LIB parallel_no_concurrency)
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/concurrency.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/concurrency.hpp" @ONLY)

################################################################################
# Enable CUDA.
option(DCA_WITH_CUDA "Enable GPU support." OFF)

if (DCA_WITH_CUDA)
  include(dca_cuda)
  if (NOT DCA_HAVE_CUDA)
    message(FATAL_ERROR "CUDA or MAGMA not found but requested.")
  endif()

  dca_add_config_define(DCA_WITH_CUDA)

  # Copy walker device config file for GPU.
  configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/walker_device_gpu.hpp"
    "${CMAKE_BINARY_DIR}/include/dca/config/walker_device.hpp")

else()
  # Copy walker device config file for CPU.
  configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/walker_device_cpu.hpp"
    "${CMAKE_BINARY_DIR}/include/dca/config/walker_device.hpp")
endif()

################################################################################
# Select the point group, the lattice type, and the model type.
# TODO: Add more point groups and lattices.

# Point group
set(DCA_POINT_GROUP "D4" CACHE STRING "Point group symmetry, options are: C6 | D4.")
set_property(CACHE DCA_POINT_GROUP PROPERTY STRINGS C6 D4)

if (DCA_POINT_GROUP STREQUAL "C6")
  set(DCA_POINT_GROUP_INCLUDE
    "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_hexagonal.hpp")

elseif (DCA_POINT_GROUP STREQUAL "D4")
  set(DCA_POINT_GROUP_INCLUDE "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp")

else()
  message(FATAL_ERROR "Please set DCA_POINT_GROUP to a valid option: C6 | D4.")
endif()

# Lattice type
set(DCA_LATTICE "square" CACHE STRING "Lattice type, options are: bilayer | square | triangular.")
set_property(CACHE DCA_LATTICE PROPERTY STRINGS bilayer square triangular)

if (DCA_LATTICE STREQUAL "bilayer")
  set(DCA_LATTICE_TYPE dca::phys::models::bilayer_lattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp")

elseif (DCA_LATTICE STREQUAL "square")
  set(DCA_LATTICE_TYPE dca::phys::models::square_lattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/square_lattice.hpp")

elseif (DCA_LATTICE STREQUAL "triangular")
  set(DCA_LATTICE_TYPE dca::phys::models::triangular_lattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/triangular_lattice.hpp")

else()
  message(FATAL_ERROR "Please set DCA_LATTICE to a valid option: bilayer | square | triangular.")
endif()

# Model type
set(DCA_MODEL "tight-binding" CACHE STRING "Model type, options are: tight-binding.")
set_property(CACHE DCA_MODEL PROPERTY STRINGS tight-binding)

if (DCA_MODEL STREQUAL "tight-binding")
  set(DCA_MODEL_TYPE dca::phys::models::TightBindingModel<Lattice>)
  set(DCA_MODEL_INCLUDE "dca/phys/models/tight_binding_model.hpp")

else()
  message(FATAL_ERROR "Please set DCA_MODEL to a valid option: tight-binding.")
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/lattice_model.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/lattice_model.hpp" @ONLY)

################################################################################
# Select the profiler type and enable auto-tuning.
set(DCA_PROFILER "None" CACHE STRING "Profiler type, options are: None | Counting | PAPI.")
set_property(CACHE DCA_PROFILER PROPERTY STRINGS None Counting PAPI)

if (DCA_PROFILER STREQUAL "Counting")
  set(DCA_PROFILING_EVENT_TYPE dca::profiling::time_event<std::size_t>)
  set(DCA_PROFILING_EVENT_INCLUDE "dca/profiling/events/time_event.hpp")
  set(DCA_PROFILER_TYPE dca::profiling::CountingProfiler<Event>)
  set(DCA_PROFILER_INCLUDE "dca/profiling/counting_profiler.hpp")

elseif (DCA_PROFILER STREQUAL "PAPI")
  # TODO: Replace long long with std::size_t?
  set(DCA_PROFILING_EVENT_TYPE "dca::profiling::papi_and_time_event<long long>")  # Need quotes because of space in 'long long'.
  set(DCA_PROFILING_EVENT_INCLUDE "dca/profiling/events/papi_and_time_event.hpp")
  set(DCA_PROFILER_TYPE dca::profiling::CountingProfiler<Event>)
  set(DCA_PROFILER_INCLUDE "dca/profiling/counting_profiler.hpp")

else()  # DCA_PROFILER = None
  # The NullProfiler doesn't have an event type.
  set(DCA_PROFILING_EVENT_TYPE void)
  set(DCA_PROFILING_EVENT_INCLUDE "dca/profiling/null_profiler.hpp")
  set(DCA_PROFILER_TYPE dca::profiling::NullProfiler)
  set(DCA_PROFILER_INCLUDE "dca/profiling/null_profiler.hpp")
endif()

option(DCA_WITH_AUTOTUNING "Enable auto-tuning. Needs a profiler type other than 'None'." OFF)
mark_as_advanced(DCA_WITH_AUTOTUNING)
if (DCA_WITH_AUTOTUNING AND NOT(DCA_PROFILER_TYPE STREQUAL "None"))
  dca_add_config_define(DCA_WITH_AUTOTUNING)
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/profiler.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/profiler.hpp" @ONLY)

################################################################################
# Select the random number generator.
set(DCA_RNG "std::mt19937_64" CACHE STRING
  "Random number generator, options are: std::mt19937_64 | std::ranlux48 | custom.")
set_property(CACHE DCA_RNG
  PROPERTY STRINGS std::mt19937_64 std::ranlux48 custom)

if (DCA_RNG STREQUAL "std::mt19937_64")
  set(DCA_RNG_TYPE dca::math::random::StdRandomWrapper<std::mt19937_64>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")
  set(DCA_RNG_LIBRARY random)

elseif (DCA_RNG STREQUAL "std::ranlux48")
  set(DCA_RNG_TYPE dca::math::random::StdRandomWrapper<std::ranlux48>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")
  set(DCA_RNG_LIBRARY random)

elseif (DCA_RNG STREQUAL "custom")
  if (NOT (DCA_RNG_CLASS AND EXISTS ${DCA_RNG_HEADER}))
    message(FATAL_ERROR
      "DCA_RNG_CLASS and DCA_RNG_HEADER must be set with the -D option, if 'custom' is chosen as RNG.")
  endif()
  set(DCA_RNG_TYPE ${DCA_RNG_CLASS})
  set(DCA_RNG_INCLUDE ${DCA_RNG_HEADER})

else()
  message(FATAL_ERROR "Please set DCA_RNG to a valid option: std::mt19937_64 | std::ranlux48 | custom.")
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/rng.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/rng.hpp" @ONLY)

################################################################################
# Select the cluster solver.
set(DCA_CLUSTER_SOLVER "CT-AUX" CACHE STRING
  "The cluster solver for the DCA(+) loop. Options are: CT-AUX | SS-CT-HYB.")
set_property(CACHE DCA_CLUSTER_SOLVER PROPERTY STRINGS CT-AUX SS-CT-HYB)

if (DCA_CLUSTER_SOLVER STREQUAL "CT-AUX")
  set(DCA_CLUSTER_SOLVER_NAME dca::phys::solver::CT_AUX)
  set(DCA_CLUSTER_SOLVER_TYPE "dca::phys::solver::CtauxClusterSolver<walker_device, ParametersType, DcaDataType>")
  set(DCA_CLUSTER_SOLVER_INCLUDE
    "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp")

elseif (DCA_CLUSTER_SOLVER STREQUAL "SS-CT-HYB")
  set(DCA_CLUSTER_SOLVER_NAME dca::phys::solver::SS_CT_HYB)
  set(DCA_CLUSTER_SOLVER_TYPE "dca::phys::solver::SsCtHybClusterSolver<walker_device, ParametersType, DcaDataType>")
  set(DCA_CLUSTER_SOLVER_INCLUDE
    "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_cluster_solver.hpp")

# elseif (DCA_CLUSTER_SOLVER STREQUAL "HTS")
#   set(DCA_CLUSTER_SOLVER_NAME dca::phys::solver::HIGH_TEMPERATURE_SERIES)
#   set(DCA_CLUSTER_SOLVER_INCLUDE
#     "dca/phys/dca_step/cluster_solver/high_temperature_series_expansion/high_temperature_series_expansion_solver.hpp")

else()
  message(FATAL_ERROR "Please set DCA_CLUSTER_SOLVER to a valid option: CT-AUX | SS-CT-HYB.")
endif()

################################################################################
# Threading options/settings
if (UNIX)
  set(DCA_THREADING_LIBS "pthread")
endif()

################################################################################
# Use threaded cluster solver.
option(DCA_WITH_THREADED_SOLVER "Use multiple walker and accumulator threads in the cluster solver." ON)

if (DCA_WITH_THREADED_SOLVER)
  dca_add_config_define(DCA_WITH_THREADED_SOLVER)
  set(DCA_THREADED_SOLVER_TYPE dca::phys::solver::StdThreadQmciClusterSolver<ClusterSolverBaseType>)
  set(DCA_THREADED_SOLVER_INCLUDE
      "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp")
endif()

################################################################################
# Enable the QMC solver built-in tests.
option(DCA_WITH_QMC_BIT "Enable QMC solver built-in tests." OFF)
mark_as_advanced(DCA_WITH_QMC_BIT)

if (DCA_WITH_QMC_BIT)
  dca_add_config_define(DCA_WITH_QMC_BIT)
endif()

################################################################################
# Single precision measurements
# TODO: change to ON by default after merging and testing the two particle accumulator.
option(DCA_WITH_SINGLE_PRECISION_MEASUREMENTS "Measure in single precision." OFF)
mark_as_advanced(DCA_WITH_SINGLE_PRECISION_MEASUREMENTS)

if (DCA_WITH_SINGLE_PRECISION_MEASUREMENTS)
  dca_add_config_define(DCA_WITH_SINGLE_PRECISION_MEASUREMENTS)
endif()

################################################################################
# Single precision coarsegraining
option(DCA_WITH_SINGLE_PRECISION_COARSEGRAINING "Coarsegrain in single precision." OFF)
mark_as_advanced(DCA_WITH_SINGLE_PRECISION_COARSEGRAINING)

if (DCA_WITH_SINGLE_PRECISION_COARSEGRAINING)
  dca_add_config_define(DCA_WITH_SINGLE_PRECISION_COARSEGRAINING)
endif()

################################################################################
# Gnuplot
option(DCA_WITH_GNUPLOT "Enable Gnuplot." OFF)

if (DCA_WITH_GNUPLOT)
  include(dca_gnuplot)
  if (NOT DCA_HAVE_GNUPLOT)
    message(FATAL_ERROR "Gnuplot not found but requested.")
  endif()

  dca_add_config_define(DCA_WITH_GNUPLOT)

  add_subdirectory(${PROJECT_SOURCE_DIR}/libs/gnuplot_i-2.11)
  set(GNUPLOT_INTERFACE_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/libs/gnuplot_i-2.11/src" CACHE INTERNAL "" FORCE)
  set(GNUPLOT_INTERFACE_LIBRARY "gnuplot_interface" CACHE INTERNAL "" FORCE)

else()
  set(GNUPLOT_INTERFACE_INCLUDE_DIR "" CACHE INTERNAL "" FORCE)
  set(GNUPLOT_INTERFACE_LIBRARY "" CACHE INTERNAL "" FORCE)
endif()

################################################################################
# Generate applications' config files.
configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/analysis.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/analysis.hpp" @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/cluster_solver_check.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/cluster_solver_check.hpp" @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/dca.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/dca.hpp" @ONLY)
