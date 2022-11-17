if(DCA_WITH_TESTS_STOCHASTIC)

  #  This application was used to build the reference data for the NiO test.
  #  leave it to break the build if its not maintained.  I found it completely inoperable after 4 years of not being built.
  add_executable(NiO_fullDCA NiO_fullDCA.cpp)
  target_link_libraries(NiO_fullDCA ${DCA_LIBS} ${DCA_THREADING_LIBS})
  target_include_directories(NiO_fullDCA PRIVATE ${DCA_INCLUDE_DIRS})
  target_compile_definitions(NiO_fullDCA PRIVATE TEST_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/")
  target_compile_definitions(NiO_fullDCA PRIVATE DCA_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

  include(FetchContent)
  # you need to have a ssh key setup with code.ornl.gov or this is just not
  # going to work
  FetchContent_Declare(
    nio_ed_test_data
    GIT_REPOSITORY "git@code.ornl.gov:dcaplusplus/dca_binary_test_data.git"
    GIT_TAG "582a07e53f081506a70e5f0cfad060e4d22e4701"
  )

  FetchContent_MakeAvailable(nio_ed_test_data)


  dca_add_gtest(NiO_stattest
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS  ${DCA_LIBS};statistical_testing
  )

  dca_add_gtest(NiO_stattest_ctint
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS  ${DCA_LIBS};statistical_testing
  )
  
  dca_add_gtest(NiO_stattest_gs
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS  ${DCA_LIBS};statistical_testing
  )

  configure_file(${nio_ed_test_data_SOURCE_DIR}/NiO.data.ed.hdf5 ${CMAKE_CURRENT_SOURCE_DIR}/NiO.data.ed.hdf5 COPYONLY)
  
  configure_file(t_ij_NiO.txt ${CMAKE_CURRENT_BINARY_DIR}/t_ij_NiO.txt COPYONLY)
  configure_file(U_ij_NiO_8_lit.txt ${CMAKE_CURRENT_BINARY_DIR}/U_ij_NiO_8_lit.txt COPYONLY)

endif()
