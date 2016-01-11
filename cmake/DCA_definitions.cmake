################################################################################
# Common definitions for building DCA++ programs and tests.
################################################################################

# Libraries
# FIXME: Handle gitVersion and modules library nicer.
set(DCA_LIBRARIES "${DCA_EXTERNAL_LIBS};${CMAKE_BINARY_DIR}/gitVersion/libgitVersion.a;${CMAKE_BINARY_DIR}/modules/libmodules.a")
# message("DCA_LIBRARIES: ${DCA_LIBRARIES}")

# Includes
set(DCA_INCLUDES
  ${CMAKE_SOURCE_DIR}/src/LAPACK_PLANS
  ${CMAKE_SOURCE_DIR}/src/comp_library/LIN_ALG
  ${CMAKE_SOURCE_DIR}/src/comp_library/IO_library
  ${CMAKE_SOURCE_DIR}/src/comp_library/function_library
  ${CMAKE_SOURCE_DIR}/src/comp_library/function_library/domains
  ${CMAKE_SOURCE_DIR}/src/comp_library/function_plotting
  ${CMAKE_SOURCE_DIR}/src/comp_library/generic_methods_library
  ${CMAKE_SOURCE_DIR}/src/comp_library/parallelization_library
  ${CMAKE_SOURCE_DIR}/src/comp_library/profiler_library
  ${CMAKE_SOURCE_DIR}/src/comp_library/type_list
  ${CMAKE_SOURCE_DIR}/src/math_library
  ${CMAKE_SOURCE_DIR}/src/math_library/random_number_library
  ${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_algorithms/compute_band_structure
  ${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_analysis
  ${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_data
  ${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_loop
  ${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_step
  ${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_step/cluster_solver
  ${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_step/symmetrization
  ${CMAKE_SOURCE_DIR}/src/phys_library/domains
  ${CMAKE_SOURCE_DIR}/src/phys_library/domains/cluster/symmetries
  ${CMAKE_SOURCE_DIR}/src/phys_library/parameters
  ${CMAKE_SOURCE_DIR}/src/phys_library/parameters/models
  )
list(APPEND DCA_INCLUDES "${CMAKE_BINARY_DIR}/include")  # Directory of include_files.h and type_definitions.h.
list(APPEND DCA_INCLUDES "${CMAKE_SOURCE_DIR}/gitVersion")  # Directory of gitVersion.hpp.
list(APPEND DCA_INCLUDES "${CMAKE_SOURCE_DIR}/modules")  # Directory of modules.hpp.
list(APPEND DCA_INCLUDES "${DCA_EXTERNAL_INCLUDES}")
# message("DCA_INCLUDES: ${DCA_INCLUDES}")
