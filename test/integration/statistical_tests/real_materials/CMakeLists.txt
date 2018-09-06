#  This application was used to build the reference data for the NiO test.
#add_executable(NiO_fullDCA NiO_fullDCA.cpp)
#target_link_libraries(NiO_fullDCA ${TEST_LIBRARIES})
#target_include_directories(NiO_fullDCA PRIVATE ${TEST_INCLUDES})
#target_compile_definitions(NiO_fullDCA PRIVATE TEST_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/")

dca_add_gtest(NiO_stattest
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS  ${DCA_LIBS};statistical_testing
    )
