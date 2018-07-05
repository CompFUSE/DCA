// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the C++ wrapper for the Gnuplot C interface.

// Defines DCA_WITH_GNUPLOT if Gnuplot is available and requested.
// Needs to be included before plot.hpp.
#include "dca/config/config_defines.hpp"

#include "dca/util/plot.hpp"

#include <complex>
#include <numeric>
#include <vector>

int main() {
  std::vector<double> x(5);
  std::iota(x.begin(), x.end(), .2);
  std::vector<double> y(5);
  std::iota(y.begin(), y.end(), -.2);
  std::vector<std::complex<double>> yc(5);
  for (size_t i = 0; i < yc.size(); ++i)
    yc[i] = std::complex<double>(.25 * i * i, -.5 * i);

  dca::util::Plot::plotPoints(x, y, "Real");
  dca::util::Plot::plotPoints(x, yc, "Complex");

  dca::util::Plot plot("lines");
  plot.plot(x, y, "Real");
  plot.plot(x, yc, "Complex");

  dca::util::Plot plot2("lines");

  plot2.setStyle("points");
  plot2.plot(x, yc, "Points");

  plot2.setStyle("lines");
  plot2.setXLabel("X");
  plot2.setYLabel("Y label");
  plot2.plot(x, y, "Lines");
}
