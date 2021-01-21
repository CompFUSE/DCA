# Function transform

add_library(function_transform STATIC
  basis_expansions.cpp boundary_conditions.cpp domain_representations.cpp element_spacings.cpp)

if (DCA_HAVE_CUDA)
  add_library(special_transform_kernels special_transforms_kernels.cu)
  set_target_properties(special_transform_kernels PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
