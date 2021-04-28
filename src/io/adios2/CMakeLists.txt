# ADIOS2

add_library(dca_adios2 STATIC adios2_writer.cpp adios2_reader.cpp)
target_link_libraries(dca_adios2 PUBLIC adios2::adios2)
target_include_directories(dca_adios2 PUBLIC adios2::adios2)
