# Function transform

add_library(function_transform STATIC
  basis_expansions.cpp boundary_conditions.cpp domain_representations.cpp element_spacings.cpp)

if (DCA_HAVE_GPU)
  add_library(special_transform_kernels special_transforms_kernels.cu)
  if (DCA_HAVE_HIP)
    set_source_files_properties(special_transform_kernels special_transforms_kernels.cu PROPERTIES LANGUAGE HIP)
  endif()
  dca_gpu_device_link(special_transform_kernels)
  target_link_libraries(function_transform PUBLIC special_transform_kernels)
endif()
