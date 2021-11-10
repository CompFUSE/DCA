################################################################################
# Find HPX
################################################################################

set(DCA_HAVE_HPX FALSE CACHE INTERNAL "")

find_package(HPX REQUIRED)
if (HPX_FOUND)
  set(DCA_HAVE_HPX TRUE)
  dca_add_haves_define(DCA_HAVE_HPX)
  include_directories(${HPX_INCLUDE_DIRS})
endif()
