################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Build options for DCA++.

################################################################################
# Enable MPI.
option(DCA_WITH_MPI "Enable MPI." ON)

if (DCA_WITH_MPI)
  if (NOT DCA_HAVE_MPI)
    message(FATAL_ERROR "MPI not found but requested.")
  endif()

  set(DCA_CONCURRENCY_TYPE "dca::concurrency::MPI_LIBRARY")
  set(DCA_CONCURRENCY_INCLUDE "dca/concurrency/parallelization_mpi.h")

else()
  set(DCA_CONCURRENCY_TYPE "dca::concurrency::SERIAL_LIBRARY")
  set(DCA_CONCURRENCY_INCLUDE "dca/concurrency/parallelization_template.h")
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/concurrency.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/concurrency.hpp" @ONLY)

################################################################################
# Enable CUDA.
option(DCA_WITH_CUDA "Enable GPU support." OFF)
option(DCA_WITH_PINNED_HOST_MEMORY "Enable pinned host memory." OFF)
mark_as_advanced(DCA_WITH_PINNED_HOST_MEMORY)

if (DCA_WITH_CUDA)
  if (NOT DCA_HAVE_CUDA)
    message(FATAL_ERROR "CUDA or MAGMA not found but requested.")
  endif()

  if (DCA_WITH_PINNED_HOST_MEMORY)
    dca_add_config_define(ENABLE_PINNED_MEMORY_ALLOCATION)
  endif()

  add_subdirectory(src/DCA_GPU_routines)

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
    "phys_library/domains/cluster/symmetries/point_groups/2D/2D_hexagonal.h")

elseif (DCA_POINT_GROUP STREQUAL "D4")
  set(DCA_POINT_GROUP_INCLUDE "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h")

else()
  message(FATAL_ERROR "Please set DCA_POINT_GROUP to a valid option: C6 | D4.")
endif()

# Lattice type
set(DCA_LATTICE "square" CACHE STRING "Lattice type, options are: bilayer | square | triangular.")
set_property(CACHE DCA_LATTICE PROPERTY STRINGS bilayer square triangular)

if (DCA_LATTICE STREQUAL "bilayer")
  set(DCA_LATTICE_TYPE bilayer_lattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "phys_library/parameters/models/analytic_hamiltonians/lattices/2D_bilayer_lattice.h")

elseif (DCA_LATTICE STREQUAL "square")
  set(DCA_LATTICE_TYPE square_lattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "phys_library/parameters/models/analytic_hamiltonians/lattices/2D_square_lattice.h")

elseif (DCA_LATTICE STREQUAL "triangular")
  set(DCA_LATTICE_TYPE triangular_lattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "phys_library/parameters/models/analytic_Hamiltonians/lattices/2D_triangular_lattice.h")

else()
  message(FATAL_ERROR "Please set DCA_LATTICE to a valid option: bilayer | square | triangular.")
endif()

# Model type
set(DCA_MODEL "tight-binding" CACHE STRING "Model type, options are: tight-binding.")
set_property(CACHE DCA_MODEL PROPERTY STRINGS tight-binding)

if (DCA_MODEL STREQUAL "tight-binding")
  set(DCA_MODEL_TYPE tight_binding_model<Lattice>)
  set(DCA_MODEL_INCLUDE "phys_library/parameters/models/tight_binding_model.h")

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
  set(DCA_PROFILING_EVENT_TYPE PROFILER::time_event<std::size_t>)
  set(DCA_PROFILING_EVENT_INCLUDE "comp_library/profiler_library/events/time_events.h")
  set(DCA_PROFILER_TYPE PROFILER::CountingProfiler<Event>)
  set(DCA_PROFILER_INCLUDE "comp_library/profiler_library/profilers/counting_profiler.hpp")

elseif (DCA_PROFILER STREQUAL "PAPI")
  # TODO: Replace long long with std::size_t?
  set(DCA_PROFILING_EVENT_TYPE PROFILER::papi_and_time_event<long long>)
  set(DCA_PROFILING_EVENT_INCLUDE "comp_library/profiler_library/events/papi_events.h")
  set(DCA_PROFILER_TYPE PROFILER::CountingProfiler<Event>)
  set(DCA_PROFILER_INCLUDE "comp_library/profiler_library/profilers/counting_profiler.hpp")

else()  # DCA_PROFILER = None
  # The NullProfiler doesn't have an event type.
  set(DCA_PROFILING_EVENT_TYPE void)
  set(DCA_PROFILING_EVENT_INCLUDE "comp_library/profiler_library/profilers/null_profiler.hpp")
  set(DCA_PROFILER_TYPE PROFILER::NullProfiler)
  set(DCA_PROFILER_INCLUDE "comp_library/profiler_library/profilers/null_profiler.hpp")
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
set(DCA_RNG "std::ranlux48_base" CACHE STRING
  "Random number generator, options are: std::ranlux48_base | std::ranlux48 | std::mt19937_64 | SPRNG_LFG | SPRNG_MLFG | Ranq2.")
set_property(CACHE DCA_RNG
  PROPERTY STRINGS std::ranlux48_base std::ranlux48 std::mt19937_64 SPRNG_LFG SPRNG_MLFG Ranq2)

if (DCA_RNG STREQUAL "std::ranlux48_base")
  set(DCA_RNG_TYPE dca::math::random::StdRandomWrapper<std::ranlux48_base>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")

elseif (DCA_RNG STREQUAL "std::ranlux48")
  set(DCA_RNG_TYPE dca::math::random::StdRandomWrapper<std::ranlux48>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")

elseif (DCA_RNG STREQUAL "std::mt19937_64")
  set(DCA_RNG_TYPE dca::math::random::StdRandomWrapper<std::mt19937_64>)
  set(DCA_RNG_INCLUDE "dca/math/random/std_random_wrapper.hpp")

elseif (DCA_RNG STREQUAL "SPRNG_LFG")
  set(DCA_RNG_TYPE dca::math::random::SprngWrapper<dca::math::random::LFG>)
  set(DCA_RNG_INCLUDE "dca/math/random/sprng_wrapper.hpp")

elseif (DCA_RNG STREQUAL "SPRNG_MLFG")
  set(DCA_RNG_TYPE dca::math::random::SprngWrapper<dca::math::random::MLFG>)
  set(DCA_RNG_INCLUDE "dca/math/random/sprng_wrapper.hpp")

elseif (DCA_RNG STREQUAL "Ranq2")
  set(DCA_RNG_TYPE dca::math::random::Ranq2)
  set(DCA_RNG_INCLUDE "dca/math/random/ranq2.hpp")

else()
  message(FATAL_ERROR "Please set DCA_RNG to a valid option: std::ranlux48_base | std::ranlux48 |
                       std::mt19937_64 | SPRNG_LFG | SPRNG_MLFG | Ranq2.")
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/rng.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/rng.hpp" @ONLY)

################################################################################
# Select the cluster solver.
set(DCA_CLUSTER_SOLVER "CT-AUX" CACHE STRING
  "Choose the cluster solver, options are: CT-AUX | SS-CT-HYB | HTS.")
set_property(CACHE DCA_CLUSTER_SOLVER PROPERTY STRINGS CT-AUX SS-CT-HYB HTS)

if (DCA_CLUSTER_SOLVER STREQUAL "CT-AUX")
  set(DCA_CLUSTER_SOLVER_NAME DCA::CT_AUX_CLUSTER_SOLVER)
  set(DCA_CLUSTER_SOLVER_INCLUDE
    "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_cluster_solver.h")

elseif (DCA_CLUSTER_SOLVER STREQUAL "SS-CT-HYB")
  set(DCA_CLUSTER_SOLVER_NAME DCA::SS_CT_HYB)
  set(DCA_CLUSTER_SOLVER_INCLUDE
    "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_solver.h")

elseif (DCA_CLUSTER_SOLVER STREQUAL "HTS")
  # TODO: Remove this if HTS solver is fixed.
  message(FATAL_ERROR "High temperature series expansion solver is not yet supported.")
  set(DCA_CLUSTER_SOLVER_NAME DCA::HIGH_TEMPERATURE_SERIES)
  set(DCA_CLUSTER_SOLVER_INCLUDE
    "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/high_temperature_series_expansion_solver.h")

else()
  message(FATAL_ERROR "Please set DCA_CLUSTER_SOLVER to a valid option: CT-AUX | SS-CT-HYB | HTS.")
endif()

################################################################################
# Select the threading library.
# TODO: Implement HPX part including DCA_HPX.cmake.
set(DCA_THREADING_LIBRARY "POSIX" CACHE STRING "Threading library, options are: POSIX | HPX.")
set_property(CACHE DCA_THREADING_LIBRARY PROPERTY STRINGS POSIX HPX)

if (DCA_THREADING_LIBRARY STREQUAL POSIX)
  if (NOT DCA_HAVE_PTHREADS)
    message(FATAL_ERROR "PThreads not found but requested.")
  endif()

  dca_add_config_define(DCA_THREADING_LIBRARY POSIX_LIBRARY)
  dca_add_config_define(DCA_THREADING_INCLUDE \"dca/concurrency/parallelization_pthreads.h\")
  set(DCA_THREADING_FLAGS -pthread CACHE STRING "Flags needed for threading.")
  mark_as_advanced(DCA_THREADING_FLAGS)

elseif (DCA_THREADING_LIBRARY STREQUAL HPX)
  message(FATAL_ERROR "No HPX support yet.")

else()
  message(FATAL_ERROR "Please set DCA_THREADING_LIBRARY to a valid option: POSIX | HPX.")
endif()

################################################################################
# Use threaded cluster solver.
option(DCA_WITH_THREADED_SOLVER "Use threaded cluster solver." ON)

if (DCA_WITH_THREADED_SOLVER)
  if (DCA_THREADING_LIBRARY STREQUAL POSIX)
    set(DCA_THREADED_SOLVER_TYPE DCA::posix_qmci_integrator<ClusterSolverBaseType>)
    set(DCA_THREADED_SOLVER_INCLUDE
      "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_pthread_jacket/posix_qmci_cluster_solver.h")

  elseif (DCA_THREADING_LIBRARY STREQUAL HPX)
    set(DCA_THREADED_SOLVER_TYPE DCA::hpx_qmci_integrator<ClusterSolverBaseType>)
    set(DCA_THREADED_SOLVER_INCLUDE
      "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_hpx_jacket/hpx_qmci_cluster_solver.h")

  else()
    message(FATAL_ERROR "Need a threading library to use a threaded cluster solver.")
  endif()
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
# Use reduced vertex function.
option(DCA_WITH_REDUCED_VERTEX_FUNCTION "Use reduced vertex function." ON)
mark_as_advanced(DCA_WITH_REDUCED_VERTEX_FUNCTION)

if (DCA_WITH_REDUCED_VERTEX_FUNCTION)
  dca_add_config_define(DCA_WITH_REDUCED_VERTEX_FUNCTION)
endif()

################################################################################
# Gnuplot
option(DCA_WITH_GNUPLOT "Enable gnuplot." OFF)

if (DCA_WITH_GNUPLOT)
  if (NOT DCA_HAVE_GNUPLOT)
    message(FATAL_ERROR "Gnuplot not found but requested.")
  endif()
endif()

################################################################################
# Generate applications' config files.
configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/analysis.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/analysis.hpp" @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/cluster_solver_check.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/cluster_solver_check.hpp" @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/dca.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/dca.hpp" @ONLY)
