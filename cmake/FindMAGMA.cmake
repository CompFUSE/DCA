if(TARGET magma::magma)
    # This module has already been processed. Don't do it again.
    return()
endif()

find_path(MAGMA_INCLUDE_DIR
  NAMES "magma_v2.h"
  DOC "Magma include directory")

find_library(MAGMA_LIBRARY
  NAMES magma libmagma)
message("MAGMA_LIBRARY: ${MAGMA_LIBRARY}")

find_library(MAGMA_sparse_LIBRARY
  NAMES magma_sparse libmagma_sparse)
message("MAGMA_sparse_LIBRARY: ${MAGMA_sparse_LIBRARY}")

mark_as_advanced(MAGMA_LIBRARY MAGMA_INCLUDE_DIR)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MAGMA
                                  REQUIRED_VARS MAGMA_LIBRARY MAGMA_sparse_LIBRARY MAGMA_INCLUDE_DIR)

if(MAGMA_FOUND)
  set(MAGMA_LIBRARIES ${MAGMA_LIBRARY} ${MAGMA_sparse_LIBRARY})
  set(MAGMA_INCLUDE_DIRS ${MAGMA_INCLUDE_DIR})
  if(NOT TARGET magma::magma)
    add_library(magma::magma UNKNOWN IMPORTED)
    set_target_properties(magma::magma PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${MAGMA_INCLUDE_DIRS}"
      INTERFACE_COMPILE_DEFINITION "${MAGMA_DEFINITIONS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C;CXX"
      IMPORTED_LOCATION "${MAGMA_LIBRARY}")
    target_link_libraries(magma::magma INTERFACE BLAS::BLAS LAPACK::LAPACK ${DCA_GPU_LIBS})
    message("Added magma::magma target")
  endif()
  if(NOT TARGET magma::sparse)
    add_library(magma::sparse UNKNOWN IMPORTED)
    set_target_properties(magma::sparse PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${MAGMA_INCLUDE_DIRS}"
      INTERFACE_COMPILE_DEFINITION "${MAGMA_DEFINITIONS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C;CXX"
      IMPORTED_LOCATION "${MAGMA_sparse_LIBRARY}")
    target_link_libraries(magma::sparse INTERFACE ${DCA_GPU_LIBS})
    message("Added magma::sparse target")
  endif()
endif()
