################################################################################
# Author: Urs R. Haehner
#
# Checks for SPRNG and accordingly sets DCA_HAVE_SPRNG.

set(DCA_HAVE_SPRNG FALSE CACHE INTERNAL "")

find_library(SPRNG_LIBRARY sprng HINTS ${SPRNG_DIR}/lib)
find_path(SPRNG_INCLUDE_DIR sprng_cpp.h HINTS ${SPRNG_DIR}/include)

mark_as_advanced(SPRNG_LIBRARY SPRNG_INCLUDE_DIR)

if (SPRNG_LIBRARY AND SPRNG_INCLUDE_DIR)
  set(DCA_HAVE_SPRNG TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_SPRNG)
endif()
