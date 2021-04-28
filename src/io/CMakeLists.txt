add_subdirectory(hdf5)
add_subdirectory(json)

if (DCA_WITH_ADIOS2 AND DCA_HAVE_ADIOS2)
add_subdirectory(adios2)
endif()
