################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for Pthreads and accordingly sets DCA_HAVE_PTHREADS.

set(DCA_HAVE_PTHREADS FALSE CACHE INTERNAL "")

set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)

if (CMAKE_USE_PTHREADS_INIT)
  set(DCA_HAVE_PTHREADS TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_PTHREADS)
endif()
