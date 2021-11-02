# Tests

# Deal with the gpu allocators needing gpu utils.  This causes many cpu tests to need them
# when DCA_HAVE_GPU=true
set(UTIL_LIBS)
if(DCA_HAVE_GPU)
set(UTIL_LIBS gpu_utils)
endif()

# Unit tests
add_subdirectory(unit)

# Integration tests
add_subdirectory(integration)

# System-level tests
add_subdirectory(system-level)

# Performance tests
add_subdirectory(performance)
