################################################################################
# Common definitions for building DCA++ programs and tests.
################################################################################

# Libraries
# FIXME: Handle gitVersion and modules library nicer.
set(DCA_LIBRARIES "${DCA_EXTERNAL_LIBS};gitVersion;modules")
# message("DCA_LIBRARIES: ${DCA_LIBRARIES}")

# Includes
set(DCA_INCLUDES
  ${PROJECT_SOURCE_DIR}/src/  # Directory of enumerations.hpp
  ${PROJECT_SOURCE_DIR}/src/LAPACK_PLANS
  ${PROJECT_SOURCE_DIR}/src/comp_library/LIN_ALG
  ${PROJECT_SOURCE_DIR}/src/comp_library/IO_library
  ${PROJECT_SOURCE_DIR}/src/comp_library/function_library
  ${PROJECT_SOURCE_DIR}/src/comp_library/function_library/domains
  ${PROJECT_SOURCE_DIR}/src/comp_library/function_plotting
  ${PROJECT_SOURCE_DIR}/src/comp_library/generic_methods_library
  ${PROJECT_SOURCE_DIR}/src/comp_library/parallelization_library
  ${PROJECT_SOURCE_DIR}/src/comp_library/profiler_library
  ${PROJECT_SOURCE_DIR}/src/comp_library/type_list
  ${PROJECT_SOURCE_DIR}/src/math_library
  ${PROJECT_SOURCE_DIR}/src/math_library/random_number_library
  ${PROJECT_SOURCE_DIR}/src/phys_library/DCA+_algorithms/compute_band_structure
  ${PROJECT_SOURCE_DIR}/src/phys_library/DCA+_analysis
  ${PROJECT_SOURCE_DIR}/src/phys_library/DCA+_data
  ${PROJECT_SOURCE_DIR}/src/phys_library/DCA+_loop
  ${PROJECT_SOURCE_DIR}/src/phys_library/DCA+_step
  ${PROJECT_SOURCE_DIR}/src/phys_library/DCA+_step/cluster_solver
  ${PROJECT_SOURCE_DIR}/src/phys_library/DCA+_step/symmetrization
  ${PROJECT_SOURCE_DIR}/src/phys_library/domains
  ${PROJECT_SOURCE_DIR}/src/phys_library/domains/cluster/symmetries
  ${PROJECT_SOURCE_DIR}/src/phys_library/parameters
  ${PROJECT_SOURCE_DIR}/src/phys_library/parameters/models
  )
list(APPEND DCA_INCLUDES "${PROJECT_SOURCE_DIR}/gitVersion")  # Directory of gitVersion.hpp.
list(APPEND DCA_INCLUDES "${PROJECT_SOURCE_DIR}/modules")  # Directory of modules.hpp.
list(APPEND DCA_INCLUDES "${DCA_EXTERNAL_INCLUDES}")
# message("DCA_INCLUDES: ${DCA_INCLUDES}")
