# HDF5
# TODO: Create one hdf5 library as soon as the multiple definition problem is fixed.

add_library(hdf5_reader STATIC hdf5_reader.cpp)
add_library(hdf5_writer STATIC hdf5_writer.cpp)

target_include_directories(hdf5_reader PUBLIC ${HDF5_INCLUDE_DIRS})
target_include_directories(hdf5_writer PUBLIC ${HDF5_INCLUDE_DIRS})
