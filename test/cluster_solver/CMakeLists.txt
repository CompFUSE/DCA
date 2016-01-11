################################################################################
# test/qmc_solver/CMakeLists.txt
################################################################################

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input.qmc_solver_test.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (DCA_GPU_SUPPORT)
  cuda_add_executable(qmc_solver_test
    ${CMAKE_CURRENT_SOURCE_DIR}/qmc_solver_test.cpp
    ${CMAKE_BINARY_DIR}/version/version.cpp
    ${CMAKE_BINARY_DIR}/includes/include_files_GPU.cu
    )
  target_link_libraries(qmc_solver_test
    ${DCA_EXTERNAL_LIBRARIES}
    ${MAGMA_LIBRARY}
    ${CUDA_CUDA_LIBRARY}
    )
  cuda_add_cublas_to_target(qmc_solver_test)

else()
  add_executable(qmc_solver_test
    ${CMAKE_CURRENT_SOURCE_DIR}/qmc_solver_test.cpp
    ${CMAKE_BINARY_DIR}/version/version.cpp
    )
  target_link_libraries(qmc_solver_test ${DCA_EXTERNAL_LIBRARIES})
endif()

target_include_directories(qmc_solver_test PRIVATE
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
# Directory of include_files.h and type_definitions.h.
target_include_directories(qmc_solver_test PRIVATE
  ${CMAKE_BINARY_DIR}/includes)

# target_include_directories(analysis PRIVATE
#   ${NFFT_DIR}/build/include)
# target_include_directories(analysis PRIVATE
#   ${SPGLIB_DIR}/build/include/spglib)

target_include_directories(qmc_solver_test PRIVATE
  $ENV{EBROOTNFFT}/include)
target_include_directories(qmc_solver_test PRIVATE
  $ENV{EBROOTSPGLIB}/include)

if (HDF5_FOUND)
  target_include_directories(qmc_solver_test PRIVATE ${HDF5_INCLUDE_DIRS})
endif()
