# HDF5

add_library(hdf5 STATIC hdf5_reader.cpp hdf5_writer.cpp)
target_include_directories(hdf5 PUBLIC ${HDF5_INCLUDE_DIRS})
