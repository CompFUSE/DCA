################################################################################
# Description.
#
# References: - https://github.com/ALPSCore/ALPSCore
################################################################################

enable_testing()

#add_subdirectory(${CMAKE_SOURCE_DIR}/libs/gmock-1.7.0)
add_subdirectory($ENV{EBROOTGTEST} ${CMAKE_BINARY_DIR}/gtest)

#target_compile_options(gmock      PRIVATE "-w")
#target_compile_options(gmock_main PRIVATE "-w")
target_compile_options(gtest      PRIVATE "-w")
target_compile_options(gtest_main PRIVATE "-w")

include_directories(${gtest_SOURCE_DIR}/include)

function(add_gtest)
  set(options ADVANCED MPI GTEST_MAIN)
  set(oneValueArgs NAME MPI_NUMPROC)
  set(multiValueArgs INCLUDES SOURCES LIBS)
  cmake_parse_arguments(ADD_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (ADD_GTEST_ADVANCED AND NOT DCA_TESTS_INCLUDE_ADVANCED)
    # Do not build this advanced test.
    return()
  endif()

  if (ADD_GTEST_MPI AND NOT DCA_MPI_AVAILABLE)
    message(WARNING "Requested MPI for test ${ADD_GTEST_NAME} but MPI is not available. Test is skipped.")
    return()
  endif()

  add_executable(${ADD_GTEST_NAME} ${ADD_GTEST_NAME}.cpp "${ADD_GTEST_SOURCES}")
  
  if (ADD_GTEST_GTEST_MAIN)
    # Use gtest main.
    target_link_libraries(${ADD_GTEST_NAME} gtest_main "${ADD_GTEST_LIBS}")
  else()
    # Test has its own main.
    target_link_libraries(${ADD_GTEST_NAME} gtest "${ADD_GTEST_LIBS}")
  endif()

  target_include_directories(${ADD_GTEST_NAME} PRIVATE "${ADD_GTEST_INCLUDES}")

  if (ADD_GTEST_MPI)
    if (NOT DEFINED ADD_GTEST_MPI_NUMPROC)
      # message("MPI_NUMPROC not specified. Setting it to 4.")
      set(ADD_GTEST_MPI_NUMPROC 4)
    endif()
    
    add_test(
      NAME    ${ADD_GTEST_NAME}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ADD_GTEST_MPI_NUMPROC} ${MPIEXEC_PREFLAGS} "./${ADD_GTEST_NAME}" ${MPIEXEC_POSTFLAGS})

  else()
    add_test(
      NAME    ${ADD_GTEST_NAME}
      COMMAND ${ADD_GTEST_NAME})

  endif()
  
endfunction()
