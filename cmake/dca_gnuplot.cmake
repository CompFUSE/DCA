################################################################################
# Author: Urs R. Haehner
#
# Checks for Gnuplot and accordingly sets DCA_HAVE_GNUPLOT.

set(DCA_HAVE_GNUPLOT FALSE CACHE INTERNAL "")

find_package(Gnuplot REQUIRED)

if (GNUPLOT_FOUND)
  set(DCA_HAVE_GNUPLOT TRUE CACHE INTERNAL "")
endif()
