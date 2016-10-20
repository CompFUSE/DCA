# Function transform

set(FUNCTION_TRANSFORM_SRC basis_expansions.cpp boundary_conditions.cpp domain_representations.cpp element_spacings.cpp)

if (DCA_HAVE_MULDEFS)
  add_library(function_transform STATIC ${FUNCTION_TRANSFORM_SRC})
else()
    add_library(function_transform SHARED ${FUNCTION_TRANSFORM_SRC})
endif()
