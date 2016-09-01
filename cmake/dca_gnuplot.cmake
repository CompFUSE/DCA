################################################################################
# Author: Urs R. Haehner
#
# Checks for Gnuplot and accordingly sets DCA_HAVE_GNUPLOT.

set(DCA_HAVE_GNUPLOT FALSE CACHE INTERNAL "")

find_package(Gnuplot)

if (GNUPLOT_FOUND)
  set(DCA_HAVE_GNUPLOT TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_GNUPLOT)
endif()
