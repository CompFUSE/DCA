# HDF5
# TODO: Create one hdf5 library as soon as the multiple definition problem is fixed.

if (DCA_HAVE_MULDEFS)
  add_library(hdf5_reader STATIC hdf5_reader.cpp)
  add_library(hdf5_writer STATIC hdf5_writer.cpp)
else()
  add_library(hdf5_reader SHARED hdf5_reader.cpp)
  add_library(hdf5_writer SHARED hdf5_writer.cpp)
  target_link_libraries(hdf5_reader ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
  target_link_libraries(hdf5_writer ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

target_include_directories(hdf5_reader PUBLIC ${HDF5_INCLUDE_DIRS})
target_include_directories(hdf5_writer PUBLIC ${HDF5_INCLUDE_DIRS})
