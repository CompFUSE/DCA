# HDF5

add_library(dca_hdf5 STATIC hdf5_reader.cpp hdf5_writer.cpp)
target_include_directories(dca_hdf5 PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(dca_hdf5 PUBLIC ${HDF5_LIBRARIES})
