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

include(FetchContent)
  
# FetchContent_Declare(
#     niotestdata
#     URL https://www.dropbox.com/s/0jbzz6e12vgu3go/NiO.data.ed.hdf5
#     URL_HASH SHA256=e9645af0b0e1e98af9ccfa808ca7895a100a372d
#     )

# FetchContent_MakeAvailable(niotestdata)
  
configure_file(t_ij_NiO.txt ${CMAKE_CURRENT_BINARY_DIR}/t_ij_NiO.txt COPYONLY)
configure_file(U_ij_NiO_8_lit.txt ${CMAKE_CURRENT_BINARY_DIR}/U_ij_NiO_8_lit.txt COPYONLY)

