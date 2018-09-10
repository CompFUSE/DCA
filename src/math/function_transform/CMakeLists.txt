# Function transform

add_library(function_transform STATIC
  basis_expansions.cpp boundary_conditions.cpp domain_representations.cpp element_spacings.cpp)

if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(special_transform_kernels special_transforms_kernels.cu)
endif()
