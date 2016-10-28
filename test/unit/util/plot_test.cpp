#include <complex>
#include <numeric>
#include <vector>
#include "dca/util/plot.hpp"

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
