# HDF5

add_library(dca_hdf5 STATIC hdf5_reader.cpp hdf5_writer.cpp)
target_link_libraries(dca_hdf5 PUBLIC ${HDF5_LIBRARIES})
target_include_directories(dca_hdf5 PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(dca_hdf5 PUBLIC ${HDF5_LIBRARIES})

message("CMAKE_CXX_COMPILER_EXTERNAL_TOOLCHAIN : ${CMAKE_CXX_STANDARD_LIBRARIES}")

# Only GCC < 9 needs to link to a separate filesystem library
target_link_libraries(dca_hdf5 PRIVATE
  std::filesystem
)
#$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>
