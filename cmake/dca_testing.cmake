################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
#
# Enables testing.
# References: - https://github.com/ALPSCore/ALPSCore

include(CMakeParseArguments)

include(ProcessorCount)
ProcessorCount(CPUS)

if(DCA_HAVE_HPX)
  set(test_thread_option HPX)
endif()

# Adds a test written with Google Test.
#
# dca_add_gtest(name
#               [FAST | EXTENSIVE | STOCHASTIC | PERFORMANCE]
#               [GTEST_MAIN]
#               [THREADED]
#               [MPI [MPI_NUMPROC procs]]
#               [CUDA | CUDA_MPI]
#               [INCLUDE_DIRS dir1 [dir2 ...]]
#               [SOURCES src1 [src2 ...]]
#               [LIBS lib1 [lib2 ...]])
#
# Adds a test called 'name', the source is assumed to be 'name.cpp'.
# The type of the test can be FAST, EXTENSIVE, STOCHASTIC or PERFORMANCE (mutually exclusive
# options). If no option is specified, the default is FAST.
# MPI or CUDA may be given to indicate that the test requires these libraries. MPI_NUMPROC is the
# number of MPI processes to use for a test with MPI, the default value is 1.
function(dca_add_gtest name)
  set(options FAST EXTENSIVE STOCHASTIC PERFORMANCE GTEST_MAIN THREADED MPI CUDA CUDA_MPI)
  set(oneValueArgs MPI_NUMPROC)
  set(multiValueArgs INCLUDE_DIRS SOURCES LIBS)
  cmake_parse_arguments(DCA_ADD_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # FAST, EXTENSIVE and PERFORMANCE are mutually exclusive.
  if ((DCA_ADD_GTEST_FAST AND DCA_ADD_GTEST_EXTENSIVE) OR
      (DCA_ADD_GTEST_FAST AND DCA_ADD_GTEST_PERFORMANCE) OR
      (DCA_ADD_GTEST_EXTENSIVE AND DCA_ADD_GTEST_PERFORMANCE) OR
      (DCA_ADD_GTEST_STOCHASTIC AND DCA_ADD_GTEST_FAST) OR
      (DCA_ADD_GTEST_STOCHASTIC AND DCA_ADD_GTEST_EXTENSIVE) OR
      (DCA_ADD_GTEST_STOCHASTIC AND DCA_ADD_GTEST_PERFORMANCE))
    message(FATAL_ERROR "Incorrect use of dca_add_gtest.\n
                         dca_add_gtest(name\n
                                       [FAST | EXTENSIVE | STOCHASTIC | PERFORMANCE]\n
                                       [GTEST_MAIN]\n
                                       [THREADED]\n
                                       [MPI [MPI_NUMPROC procs]]\n
                                       [CUDA | CUDA_MPI]\n
                                       [INCLUDE_DIRS dir1 [dir2 ...]]\n
                                       [SOURCES src1 [src2 ...]]\n
                                       [LIBS lib1 [lib2 ...]])")
  endif()

  # Only build the test if the corresponding option is set.
  if (DCA_ADD_GTEST_PERFORMANCE)
    if (NOT DCA_WITH_TESTS_PERFORMANCE)
      return()
    endif()
    # Only build performance tests in Release mode.
    # if (NOT (CMAKE_BUILD_TYPE STREQUAL "Release"))
    #   return ()
    # endif()

  elseif (DCA_ADD_GTEST_EXTENSIVE)
    if (NOT DCA_WITH_TESTS_EXTENSIVE)
      return()
    endif()

  elseif (DCA_ADD_GTEST_STOCHASTIC)
    if (NOT DCA_WITH_TESTS_STOCHASTIC)
      return()
    endif()

  else()  # Default is FAST.
    if (NOT DCA_WITH_TESTS_FAST)
      return()
    endif()
  endif()

  # Only build the test if the required libraries are available.
  if (DCA_ADD_GTEST_MPI AND NOT DCA_HAVE_MPI)
    return()
  endif()

  # If we don't have GPU these tests aren't added.
  if (DCA_ADD_GTEST_CUDA AND NOT DCA_HAVE_GPU)
    return()
  endif()

  # Right now we're only testing GPU distributed code on one node so its pointless
  # without more than one GPU per node.
  if (DCA_ADD_GTEST_CUDA_MPI AND DCA_HAVE_GPU_AWARE_MPI AND (DCA_TEST_GPU_COUNT LESS 2) )
    return()
  endif()

  if (DCA_ADD_GTEST_GTEST_MAIN)
    set(DCA_ADD_GTEST_SOURCES ${PROJECT_SOURCE_DIR}/test/dca_gtest_main.cpp ${DCA_ADD_GTEST_SOURCES})
  endif()

  add_executable(${name} ${name}.cpp ${DCA_ADD_GTEST_SOURCES})

  target_compile_definitions(${name} PRIVATE DCA_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")

  # this is hacky but allows continued use of DCA_THREADING_LIBS
  # if (DCA_ADD_GTEST_LIBS MATCHES "parallel_hpx")
  #   set(oldname ${name})
  #   set(name ${name}_hpx)
  #   add_executable(${name} ${oldname}.cpp ${DCA_ADD_GTEST_SOURCES})
  #   hpx_setup_target(${name})
  #   set(DCA_TESTING_ARGS_HPX "--hpx:threads=5")
  #   if (DEFINED DCA_ADD_GTEST_MPI_NUMPROC)
  #     if(DCA_ADD_GTEST_MPI_NUMPROC EQUAL 0)
  #       set(DCA_ADD_GTEST_MPI_NUMPROC 1)
  #     endif()
  #     math(EXPR NUMTHREADS "${CPUS} / ${DCA_ADD_GTEST_MPI_NUMPROC}")
  #     if ("${NUMTHREADS}" LESS "1")
  #       set(NUMTHREADS "1")
  #     endif()
  #     set(DCA_TESTING_ARGS_HPX "--hpx:threads=5")
  #   endif()
  # endif()

  # Create a macro with the project source dir. We use this as the root path for reading files in
  # tests.

  target_link_libraries(${name} PUBLIC ${DCA_ADD_GTEST_LIBS} gtest)

  if (DCA_HAVE_CUDA)
    message("test has cuda!")
    set_property(TARGET ${name} PROPERTY CMAKE_CUDA_ARCHITECTURES 70)
    set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_target_properties(${name} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    set_target_properties(${name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${name} PRIVATE CUDA::cudart)
  endif()

  if (DCA_HAVE_HIP)
    target_link_libraries(${name} PUBLIC hip::host)
  endif()

  if (DCA_ADD_GTEST_CUDA OR DCA_ADD_GTEST_CUDA_MPI)
    IF(DCA_HAVE_CUDA)
      target_include_directories(${name} PRIVATE ${CUDATookit_INCLUDE_DIRS})
      target_compile_definitions(${name} PRIVATE DCA_HAVE_CUDA)
      target_link_libraries(${name} PRIVATE ${DCA_GPU_LIBS})
      target_link_libraries(${name} PRIVATE CUDA::cublas)
    ENDIF()
    if(DCA_HAVE_MAGMA)
      target_include_directories(${name} PRIVATE ${MAGMA_INCLUDE_DIR})
      target_compile_definitions(${name} PRIVATE DCA_HAVE_MAGMA)
    endif()
    if(DCA_WITH_GPU_AWARE_MPI)
      target_compile_definitions(${name} PRIVATE DCA_HAVE_GPU_AWARE_MPI)
    endif()
    if (DCA_ADD_GTEST_CUDA_MPI)
      #We need to document which cuda aware openmpi requires this and which doesn't.
      set(THIS_CVD_LAUNCHER "${CMAKE_SOURCE_DIR}/${CVD_LAUNCHER}")
    endif()
  endif()

  if (DCA_ADD_GTEST_THREADED AND DCA_HAVE_HPX)
    message("adding HPX link")
    set(DCA_TESTING_FLAGS "--hpx:threads=5")
    target_compile_definitions(${name} PUBLIC DCA_HAVE_HPX)
    target_link_libraries(${name} PUBLIC HPX::hpx HPX::wrap_main)
  endif()
  
  target_include_directories(${name} PRIVATE
    ${gtest_SOURCE_DIR}/include
    ${DCA_ADD_GTEST_INCLUDE_DIRS})

  if (DCA_ADD_GTEST_MPI)
    if (NOT DEFINED DCA_ADD_GTEST_MPI_NUMPROC)
      set(DCA_ADD_GTEST_MPI_NUMPROC 1)
    endif()

    add_test(NAME ${name}
             COMMAND ${TEST_RUNNER} ${MPIEXEC_NUMPROC_FLAG} ${DCA_ADD_GTEST_MPI_NUMPROC}
                     ${MPIEXEC_PREFLAGS} ${SMPIARGS_FLAG_MPI} ${CVD_LAUNCHER} "$<TARGET_FILE:${name}>"
                     ${DCA_TESTING_FLAGS})
                 target_link_libraries(${name} PRIVATE ${MPI_C_LIBRARIES})
  else()
    if (TEST_RUNNER)
      add_test(NAME ${name}
               COMMAND ${TEST_RUNNER} ${MPIEXEC_NUMPROC_FLAG} 1
                   ${MPIEXEC_PREFLAGS} ${SMPIARGS_FLAG_NOMPI} ${CVD_LAUNCHER} "$<TARGET_FILE:${name}>"
                   ${DCA_TESTING_FLAGS})
    else (TEST_RUNNER)
      add_test(NAME ${name}
               COMMAND "$<TARGET_FILE:${name}>"
               ${DCA_TESTING_FLAGS})
    endif (TEST_RUNNER)
  endif()

endfunction()

# Adds a performance test written with the Google benchmark.
#
# dca_add_gtest(name
#               [INCLUDE_DIRS dir1 [dir2 ...]]
#               [LIBS lib1 [lib2 ...]])
#
# If DCA_WITH_TESTS_PERFORMANCE is defined, adds a test called 'name', the source is assumed to be
# 'name.cpp'.

function(dca_add_perftest name)
  set(multiValueArgs INCLUDE_DIRS LIBS)
  cmake_parse_arguments(DCA_ADD_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (NOT DCA_WITH_TESTS_PERFORMANCE)
    return()
  endif()

  add_executable(${name} ${name}.cpp)

  target_include_directories(${name} PRIVATE ${benchmark_dir}/include ${PROJECT_SOURCE_DIR}/include
                             ${DCA_ADD_PERFTEST_INCLUDE_DIRS})


  target_link_libraries(${name} PRIVATE ${DCA_ADD_PERFTEST_LIBS};gtest benchmark_main benchmark)

  target_compile_definitions(${name} PRIVATE DCA_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
endfunction()


