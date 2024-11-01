################################################################################
# Copyright (C) 2023 ETH Zurich
# Copyright (C) 2023 UT-Battelle, LLC
# All rights reserved.
#
# See LICENSE for terms of usage.
# See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
#
# Authors: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#          Giovanni Badlduzzi (gbalduzz@itp.phys.ethz.ch)
#          Peter Doak (doakpw@ornl.gov
#
# Build options for DCA++.

################################################################################
# Enable MPI.
option(DCA_WITH_MPI "Enable MPI." ON)
option(DCA_WITH_GPU_AWARE_MPI "Enable GPU aware MPI." OFF)

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

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/threading.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/threading.hpp" @ONLY)

################################################################################
# No GPU by default -- DCA_WITH_CUDA or DCA_WITH_HIP will set true if they
#                      are found
set(DCA_HAVE_GPU FALSE CACHE INTERNAL "")

################################################################################
# Enable CUDA.
option(DCA_WITH_CUDA "Enable Nvidia GPU support." OFF)

if (DCA_WITH_CUDA)
  include(dca_cuda)
  if (NOT DCA_HAVE_CUDA)
    message(FATAL_ERROR "CUDA and/or MAGMA not found but requested.")
  endif()
  dca_add_config_define(DCA_WITH_CUDA)
  dca_add_config_define(DCA_WITH_GPU)
endif()

################################################################################
# Enable HIP
option(DCA_WITH_HIP "Enable AMD GPU support." OFF)

if (DCA_WITH_HIP)
  include(dca_hip)
  if (NOT DCA_HAVE_HIP)
    message(FATAL_ERROR "HIP and/or MAGMA not found but requested.")
  endif()
  dca_add_config_define(DCA_WITH_HIP)
  dca_add_config_define(DCA_WITH_GPU)
endif()

################################################################################
# Extra parameters common to all GPU setups.
if (DCA_HAVE_GPU)
  if (NOT DCA_HAVE_MAGMA)
    message(FATAL_ERROR "At the moment the GPU code requires MAGMA.")
  endif()
  if (DCA_WITH_GPU_AWARE_MPI)
    # Don't have a good way to test for this at the moment.
    dca_add_haves_define(DCA_HAVE_GPU_AWARE_MPI)
  endif()
  # copy the appropriate walker_device header (defines device template parameter)
  configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/walker_device_gpu.hpp"
    "${CMAKE_BINARY_DIR}/include/dca/config/walker_device.hpp")
else()
  configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/walker_device_cpu.hpp"
    "${CMAKE_BINARY_DIR}/include/dca/config/walker_device.hpp")
endif()

################################################################################
# Enable ADIOS2.
option(DCA_WITH_ADIOS2 "Enable ADIOS2 support." OFF)
if (DCA_WITH_ADIOS2)
  dca_add_config_define(DCA_WITH_ADIOS2)
  dca_add_haves_define(DCA_HAVE_ADIOS2)
endif()
################################################################################
# Select the point group, the lattice type, and the model type.
# TODO: Add more point groups and lattices.

# Point group
set(DCA_POINT_GROUP "D4" CACHE STRING "Point group symmetry, options are: C6 | D4 | no_symmetry<2> | no_symmetry<3>.")
set_property(CACHE DCA_POINT_GROUP PROPERTY STRINGS C6 D4 no_symmetry<2> no_symmetry<3>)

if (DCA_POINT_GROUP STREQUAL "C6")
  set(DCA_POINT_GROUP_INCLUDE
    "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_hexagonal.hpp")

elseif (DCA_POINT_GROUP STREQUAL "D4")
  set(DCA_POINT_GROUP_INCLUDE "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp")

elseif (DCA_POINT_GROUP STREQUAL "no_symmetry<2>")
  set(DCA_POINT_GROUP_INCLUDE "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp")
elseif (DCA_POINT_GROUP STREQUAL "no_symmetry<3>")
  set(DCA_POINT_GROUP_INCLUDE "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp")

else()
  message(FATAL_ERROR "Please set DCA_POINT_GROUP to a valid option: C6 | D4 | no_symmetry<2> | no_symmetry<3>.")
endif()

# Lattice type
set(DCA_LATTICE "square" CACHE STRING "Lattice type, options are: bilayer | square | triangular |
    Kagome | hund | twoband_Cu | threeband | Rashba_Hubbard | Moire_Hubbard | FeAs | material_NiO | material_FeSn ")
set_property(CACHE DCA_LATTICE PROPERTY STRINGS bilayer square triangular Kagome hund twoband_Cu threeband
             Rashba_Hubbard Moire_Hubbard FeAs material_NiO material_FeSn)

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
elseif (DCA_LATTICE STREQUAL "Kagome")
  set(DCA_LATTICE_TYPE dca::phys::models::KagomeHubbard<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/Kagome_hubbard.hpp")
elseif (DCA_LATTICE STREQUAL "hund")
  set(DCA_LATTICE_TYPE dca::phys::models::HundLattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp")
elseif (DCA_LATTICE STREQUAL "threeband")
  set(DCA_LATTICE_TYPE dca::phys::models::ThreebandHubbard<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/threeband_hubbard.hpp")
elseif (DCA_LATTICE STREQUAL "Rashba_Hubbard")
  set(DCA_LATTICE_TYPE dca::phys::models::RashbaHubbard<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/rashba_hubbard.hpp")
elseif (DCA_LATTICE STREQUAL "Moire_Hubbard")
  set(DCA_LATTICE_TYPE dca::phys::models::moire_hubbard<PointGroup>)
  set(DCA_LATTICE_INCLUDE
    "dca/phys/models/analytic_hamiltonians/Moire_Hubbard.hpp")
elseif (DCA_LATTICE STREQUAL "twoband_chain")
  set(DCA_LATTICE_TYPE dca::phys::models::twoband_chain<dca::phys::domains::no_symmetry<1>>)
  set(DCA_LATTICE_INCLUDE
      "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp")
elseif (DCA_LATTICE STREQUAL "FeAs")
  set(DCA_LATTICE_TYPE dca::phys::models::FeAsLattice<PointGroup>)
  set(DCA_LATTICE_INCLUDE
      "dca/phys/models/analytic_hamiltonians/fe_as_lattice.hpp")
elseif (DCA_LATTICE STREQUAL "twoband_Cu")
  set(DCA_LATTICE_TYPE dca::phys::models::TwoBandCu<PointGroup>)
  set(DCA_LATTICE_INCLUDE
      "dca/phys/models/analytic_hamiltonians/twoband_Cu.hpp")
elseif (DCA_LATTICE STREQUAL "material_NiO")
  set(DCA_LATTICE_TYPE "dca::phys::models::material_lattice<dca::phys::models::Material::NiO_unsymmetric, dca::phys::domains::${DCA_POINT_GROUP}>")
  set(DCA_LATTICE_INCLUDE
      "dca/phys/models/material_hamiltonians/material_lattice.hpp")
  set(DCA_MODEL_IS_MATERIAL_LATTICE ON CACHE BOOL "is the model a material lattice")
elseif (DCA_LATTICE STREQUAL "material_FeSn")
  set(DCA_LATTICE_TYPE "dca::phys::models::material_lattice<dca::phys::models::Material::FeSn, dca::phys::domains::${DCA_POINT_GROUP}>")
  set(DCA_LATTICE_INCLUDE
      "dca/phys/models/material_hamiltonians/material_lattice.hpp")
  set(DCA_MODEL_IS_MATERIAL_LATTICE ON CACHE BOOL "is the model a material lattice")
else()
  message(FATAL_ERROR "Please set DCA_LATTICE to a valid option: bilayer | square | triangular | Kagome | hund | twoband_Cu | threeband | Rashba_Hubbard | Moire_Hubbard | FeAs | material_NiO | material_FeSn.")
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
set(DCA_PROFILER "None" CACHE STRING "Profiler type, options are: None | Counting | PAPI | Cuda.")
set_property(CACHE DCA_PROFILER PROPERTY STRINGS None Counting PAPI Cuda)

if (DCA_PROFILER STREQUAL "Counting")
  set(DCA_PROFILING_EVENT_TYPE dca::profiling::time_event<std::size_t>)
  set(DCA_PROFILING_EVENT_INCLUDE "dca/profiling/events/time_event.hpp")
  set(DCA_PROFILER_TYPE dca::profiling::CountingProfiler<Event>)
  set(DCA_PROFILER_INCLUDE "dca/profiling/counting_profiler.hpp")

elseif (DCA_PROFILER STREQUAL "PAPI")
  set(DCA_PROFILING_EVENT_TYPE "dca::profiling::PapiAndTimeEvent")
  set(DCA_PROFILING_EVENT_INCLUDE "dca/profiling/events/papi_and_time_event.hpp")
  set(DCA_PROFILER_TYPE dca::profiling::CountingProfiler<Event>)
  set(DCA_PROFILER_INCLUDE "dca/profiling/counting_profiler.hpp")

# Note: this profiler requires using the PTHREAD library and CUDA_TOOLS_EXT_LIBRARY
elseif (DCA_PROFILER STREQUAL "Cuda")
  set(DCA_PROFILING_EVENT_INCLUDE "dca/profiling/events/time.hpp")
  set(DCA_PROFILING_EVENT_TYPE "void")
  set(DCA_PROFILER_TYPE dca::profiling::CudaProfiler)
  set(DCA_PROFILER_INCLUDE "dca/profiling/cuda_profiler.hpp")
  link_libraries(${CUDA_nvToolsExt_LIBRARY})

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
  "The cluster solver for the DCA(+) loop. Options are: CT-AUX | CT-INT | SS-CT-HYB.")
set_property(CACHE DCA_CLUSTER_SOLVER PROPERTY STRINGS CT-AUX CT-INT SS-CT-HYB)

if (DCA_CLUSTER_SOLVER STREQUAL "CT-INT")
  set(DCA_CLUSTER_SOLVER_NAME dca::ClusterSolverId::CT_INT)
  set(DCA_CLUSTER_SOLVER_INCLUDE "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp")

  set(DCA_USE_CTINT_SUBMATRIX ON CACHE BOOL "Use submatrix updates if the CT-INT solver is selected.")
  if(DCA_USE_CTINT_SUBMATRIX)
    set(DCA_CLUSTER_SOLVER_TYPE
            "dca::phys::solver::CtintClusterSolver<walker_device, ParametersType, true, DIST>")
  else()
    set(DCA_CLUSTER_SOLVER_TYPE
            "dca::phys::solver::CtintClusterSolver<walker_device, ParametersType, false, DIST>")
  endif()

elseif (DCA_CLUSTER_SOLVER STREQUAL "CT-AUX")
  set(DCA_CLUSTER_SOLVER_NAME dca::ClusterSolverId::CT_AUX)
  set(DCA_CLUSTER_SOLVER_TYPE "dca::phys::solver::CtauxClusterSolver<walker_device, ParametersType, DcaDataType<DIST>, DIST>")
  set(DCA_CLUSTER_SOLVER_INCLUDE
      "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp")


elseif (DCA_CLUSTER_SOLVER STREQUAL "SS-CT-HYB")
  set(DCA_CLUSTER_SOLVER_NAME dca::ClusterSolverId::SS_CT_HYB)
  set(DCA_CLUSTER_SOLVER_TYPE "dca::phys::solver::SsCtHybClusterSolver<walker_device, ParametersType, DcaDataType<DIST>, DIST>")
  set(DCA_CLUSTER_SOLVER_INCLUDE
        "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_cluster_solver.hpp")

# elseif (DCA_CLUSTER_SOLVER STREQUAL "HTS")
#   set(DCA_CLUSTER_SOLVER_NAME dca::phys::solver::HIGH_TEMPERATURE_SERIES)
#   set(DCA_CLUSTER_SOLVER_INCLUDE
#     "dca/phys/dca_step/cluster_solver/high_temperature_series_expansion/high_temperature_series_expansion_solver.hpp")

else()
  message(FATAL_ERROR "Please set DCA_CLUSTER_SOLVER to a valid option: CT-AUX | CT_INT |
          SS-CT-HYB.")
endif()

################################################################################
# Threading options/settings with gtest this is nonoptional
if (UNIX)
  set(DCA_THREADING_LIBS pthread)
endif()

################################################################################
# Use threaded cluster solver.
option(DCA_WITH_THREADED_SOLVER "Use multiple walker and accumulator threads in the cluster solver." ON)

if (DCA_WITH_THREADED_SOLVER)
  dca_add_config_define(DCA_WITH_THREADED_SOLVER)
  set(DCA_THREADED_SOLVER_TYPE dca::phys::solver::StdThreadQmciClusterSolver<ClusterSolverBaseType<DIST>>)
  set(DCA_THREADED_SOLVER_INCLUDE
    "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp")
endif()

################################################################################
# Enable HPX threading support if desired
option(DCA_WITH_HPX "Enable HPX for multi-threading" OFF)
if (DCA_WITH_HPX)
  # if HPX is not found then DCA_HAVE_HPX will not be set
  include(dca_hpx)
  if (NOT DCA_HAVE_HPX)
    message(FATAL_ERROR "HPX library not found but requested.")
  endif()
  if (DCA_WITH_THREADED_SOLVER)
    set(DCA_THREADING_LIBS ${DCA_THREADING_LIBS} parallel_hpx)
  endif()
else()
  if (DCA_WITH_THREADED_SOLVER)
    set(DCA_THREADING_LIBS ${DCA_THREADING_LIBS} parallel_stdthread)
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
# MC options.
option(DCA_WITH_MEMORY_SAVINGS "Save memory in the two particle accumulation at a slight performance
       cost." OFF)
mark_as_advanced(DCA_WITH_MEMORY_SAVINGS)
if (DCA_WITH_MEMORY_SAVINGS)
  set(MEMORY_SAVINGS true)
else()
  set(MEMORY_SAVINGS false)
endif()

option(DCA_WITH_SINGLE_PRECISION_MC "Perform Monte Carlo and measurements in single precision." OFF)
option(DCA_WITH_SINGLE_PRECISION_TP_MEASUREMENTS "Measure two particle function in single precision." OFF)

if (DCA_WITH_SINGLE_PRECISION_MC)
  set(DCA_WITH_SINGLE_PRECISION_TP_MEASUREMENTS ON CACHE BOOL "Measure two particle function in single precision." FORCE)
  set(MC_REAL float)
else()
  set(MC_REAL double)
endif()

if (DCA_WITH_SINGLE_PRECISION_TP_MEASUREMENTS)
  set(TP_ACCUMULATION_PRECISION float)
else()
  set(TP_ACCUMULATION_PRECISION double)
endif()


option(DCA_WITH_MANAGED_MEMORY "Use managed memory allocator." OFF)
mark_as_advanced(DCA_WITH_MANAGED_MEMORY)
if (DCA_WITH_MANAGED_MEMORY)
  set(TWO_PARTICLE_ALLOCATOR "dca::linalg::util::ManagedAllocator<T>")
else()
  set(TWO_PARTICLE_ALLOCATOR "dca::linalg::util::DeviceAllocator<T>")
endif()

option(DCA_WITH_CTAUX_TRACING "special debug tracing of of delayed spin updates in ctaux" OFF)
mark_as_advanced(DCA_WITH_CTAUX_TRACING)
if(DCA_WITH_CTAUX_TRACING)
  add_compile_definitions(CTAUX_DEBUG_TRACING)
endif()

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/mc_options.hpp.in"
        "${CMAKE_BINARY_DIR}/include/dca/config/mc_options.hpp" @ONLY)

################################################################################
# Workarounds
option(DCA_FIX_BROKEN_MPICH "Re-define MPI_CXX_* datatypes as the corresponding MPI_C_* datatypes when mpich is the mpi provider."
       OFF)
if(DCA_FIX_BROKEN_MPICH)
  add_compile_definitions(DCA_FIX_BROKEN_MPICH)
endif()

if ((DCA_LATTICE STREQUAL "material") AND (NOT DCA_POINT_GROUP STREQUAL "no_symmetry<3>"))
  message( FATAL_ERROR "material lattice must be used with the no_symmetry<3> pointgroup")
endif()
################################################################################
# Generate applications' config files.
configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/analysis.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/analysis.hpp" @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/cluster_solver_check.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/cluster_solver_check.hpp" @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/dca.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/dca.hpp" @ONLY)
