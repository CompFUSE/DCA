################################################################################
# applications/analysis/test/CMakeLists.txt
################################################################################

# TODO: Avoid copying this file.
file(COPY input_data.analysis_DCA+_mpi_test.hdf5
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (DCA_WITH_MPI)
add_gtest(NAME analysis_DCA+_mpi_test
          ADVANCED
          MPI MPI_NUMPROC 8
          INCLUDES ${DCA_INCLUDES}
          LIBS ${DCA_LIBRARIES})
endif()
