// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides a C++ wrapper for the Gnuplot C interface.
//
// TODO: const correctness for dca::func::functions.

#ifndef DCA_UTIL_PLOT_HPP
#define DCA_UTIL_PLOT_HPP

#include <cassert>
#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

#ifdef DCA_WITH_GNUPLOT
extern "C" {
#include "gnuplot_i.h"
}
#endif  // DCA_WITH_GNUPLOT

namespace dca {
namespace util {
// dca::util::

#ifdef DCA_WITH_GNUPLOT
class Plot {
public:
  Plot(std::string style = "points") : gp_handle_(gnuplot_init()) {
    gnuplot_setstyle(gp_handle_, &style[0]);
  }

  Plot(const Plot& rhs) = delete;

  Plot(Plot&& rhs) : gp_handle_(rhs.gp_handle_) {
    rhs.gp_handle_ = nullptr;
  }

  Plot& operator=(const Plot& rhs) = delete;

  Plot& operator=(Plot&& rhs) {
    if (gp_handle_ != nullptr)
      gnuplot_close(gp_handle_);
    gp_handle_ = rhs.gp_handle_;
    rhs.gp_handle_ = nullptr;
    return *this;
  }

  ~Plot() {
    if (gp_handle_ != nullptr)
      gnuplot_close(gp_handle_);
  }

  template <typename ScalarTypeX, typename ScalarTypeY>
  void plot(const std::vector<ScalarTypeX>& x, const std::vector<ScalarTypeY>& y,
            std::string label = "") {
    std::vector<double> x_copy = copyToDouble(x);
    auto y_copy = copyToDouble(y);
    plotReference(x_copy, y_copy, label);
  }

  template <typename ScalarType, typename DmnType>
  void plot(/*const*/ func::function<ScalarType, func::dmn_0<DmnType>>& f) {
    std::vector<double> x = copyToDouble(DmnType::get_elements());
    auto y = copyToDouble(f.size(), &f(0), 1);

    setXLabel(DmnType::get_name());
    setYLabel(f.get_name());
    plotReference(x, y, f.get_name());
  }

  template <typename ScalarType, typename DmnType0, typename DmnType1, typename DmnType2>
  void plot(
      /*const*/ func::function<ScalarType, func::dmn_variadic<DmnType0, DmnType0, DmnType1, DmnType2>>& f) {
    std::vector<double> x = copyToDouble(DmnType2::get_elements());
    int inc = &f(0, 0, 0, 0) - &f(0, 0, 0, 1);

    setXLabel(DmnType2::parameter_type::get_name());
    setYLabel(f.get_name());

    for (int j = 0; j < DmnType0::dmn_size(); j++) {
      auto y = copyToDouble(DmnType2::dmn_size(), &f(0, 0, j, 0), inc);
      plotReference(x, y, "");
    }
  }

  template <typename ScalarType, typename DmnType0, typename DmnType1, typename DmnType2>
  void plotBands(
      /*const*/ func::function<ScalarType, func::dmn_variadic<DmnType0, DmnType0, DmnType1, DmnType2>>& f) {
    int start = 0;
    int size = (DmnType2::dmn_size() + 1) / 2;

    std::vector<double> x = copyToDouble(size, &DmnType2::get_elements()[0], 1);
    int inc = &f(0, 0, 0, 0) - &f(0, 0, 0, 1);

    setXLabel(DmnType2::parameter_type::get_name());
    setYLabel(f.get_name());

    for (int j = 0; j < DmnType0::dmn_size(); j++) {
      const ScalarType* ptr = &f(j, j, 0, start);
      auto y = copyToDouble(size, ptr, inc);
      plotReference(x, y, "");
    }
  }

  template <typename ScalarType>
  void plotLine2D(std::vector<ScalarType> p1, std::vector<ScalarType> p2) {
    assert(p1.size() == 2);
    assert(p2.size() == 2);

    std::vector<double> x(2);
    std::vector<double> y(2);

    x[0] = p1[0];
    y[0] = p1[1];
    x[1] = p2[0];
    y[1] = p2[1];

    plotReference(x, y);
  }

  template <typename ScalarType>
  void plotLine3D(std::vector<ScalarType> p1, std::vector<ScalarType> p2) {
    assert(p1.size() == 3);
    assert(p2.size() == 3);
    // TODO: gnuplot C interface doesn't have plot_xyz.
  }

  void setStyle(std::string style) {
    gnuplot_setstyle(gp_handle_, &style[0]);
  }

  void setXLabel(std::string label) {
    gnuplot_set_xlabel(gp_handle_, &label[0]);
  }

  void setYLabel(std::string label) {
    gnuplot_set_ylabel(gp_handle_, &label[0]);
  }

private:
  template <typename ScalarType>
  std::vector<double> copyToDouble(const std::vector<ScalarType>& x) {
    return std::vector<double>(x.begin(), x.end());
  }
  template <typename ScalarType>
  std::pair<std::vector<double>, std::vector<double>> copyToDouble(
      const std::vector<std::complex<ScalarType>>& x) {
    return copyToDouble(x.size(), &x[0], 1);
  }
  template <typename ScalarType>
  std::vector<double> copyToDouble(size_t size, const ScalarType* x, int incx) {
    std::vector<double> cp(size);
    for (size_t i = 0; i < size; ++i) {
      cp[i] = x[i * incx];
    }
    return cp;
  }
  template <typename ScalarType>
  std::pair<std::vector<double>, std::vector<double>> copyToDouble(size_t size,
                                                                   const std::complex<ScalarType>* x,
                                                                   int incx) {
    std::vector<double> re(size);
    std::vector<double> im(size);
    for (size_t i = 0; i < size; ++i) {
      re[i] = x[i * incx].real();
      im[i] = x[i * incx].imag();
    }
    return std::make_pair(std::move(re), std::move(im));
  }

  void plotReference(std::vector<double>& x, std::vector<double>& y, std::string label = "") {
    assert(x.size() == y.size());
    gnuplot_plot_xy(gp_handle_, &x[0], &y[0], x.size(), &label[0]);
  }
  void plotReference(std::vector<double>& x, std::pair<std::vector<double>, std::vector<double>>& y,
                     std::string label = "") {
    assert(x.size() == y.first.size());
    assert(x.size() == y.second.size());
    plotReference(x, y.first, "Re[" + label + "]");
    plotReference(x, y.second, "Im[" + label + "]");
  }

public:
  template <typename ScalarType>
  static void heatMap(const std::vector<ScalarType>& x, const std::vector<ScalarType>& y,
                      const std::vector<ScalarType>& z, std::string /*label*/ = "") {
    assert(x.size() == y.size());
    assert(x.size() == z.size());
    // TODO: gnuplot C interface doesn't have plot_xyz.
    // Commands to generate the heatmap:
    // set dgrid3d  100,100,16
    // set pm3d at b
    // "gnuplot_plot_xyz"(x, y, z, label)
    // unset surface
    // set pm3d map
  }

  template <typename ScalarType, typename DmnType>
  static void plotErrorBars(func::function<ScalarType, func::dmn_0<DmnType>>& /*f*/,
                            func::function<ScalarType, func::dmn_0<DmnType>>& /*g*/) {
    // TODO: gnuplot C interface doesn't have plot_xy_err.
  }

  template <typename ScalarTypeX, typename ScalarTypeY>
  static void plotPoints(const std::vector<ScalarTypeX>& x, const std::vector<ScalarTypeY>& y,
                         std::string label = "") {
    Plot plot("points");
    plot.plot(x, y, label);
  }

  template <typename ScalarType, typename DmnType>
  static void plotLine(/*const*/ func::function<ScalarType, func::dmn_0<DmnType>>& f) {
    Plot plot("lines");
    plot.plot(f);
  }

  template <typename ScalarType, typename DmnType0, typename DmnType1, typename DmnType2>
  static void plotLinesPoints(
      /*const*/ func::function<ScalarType, func::dmn_variadic<DmnType0, DmnType0, DmnType1, DmnType2>>& f) {
    Plot plot("linespoints");
    plot.plot(f);
  }

  template <typename ScalarType, typename DmnType0, typename DmnType1, typename DmnType2>
  static void plotBandsLinesPoints(
      /*const*/ func::function<ScalarType, func::dmn_variadic<DmnType0, DmnType0, DmnType1, DmnType2>>& f) {
    Plot plot("linespoints");
    plot.plotBands(f);
  }

private:
  gnuplot_ctrl* gp_handle_;
};

#else
class Plot {
public:
  Plot(std::string = "") {}

  Plot(const Plot& rhs) = delete;

  Plot(Plot&&) {}

  Plot& operator=(const Plot& rhs) = delete;

  Plot& operator=(Plot&&) {
    return *this;
  }

  template <typename... Args>
  void plot(Args...) {}
  template <typename... Args>
  void plotBands(Args...) {}
  template <typename... Args>
  void plotLine2D(Args...) {}
  template <typename... Args>
  void plotLine3D(Args...) {}
  template <typename... Args>
  void setStyle(Args...) {}
  template <typename... Args>
  void setXLabel(Args...) {}
  template <typename... Args>
  void setYLabel(Args...) {}

public:
  template <typename... Args>
  static void heatMap(Args...) {}
  template <typename... Args>
  static void plotErrorBars(Args...) {}
  template <typename... Args>
  static void plotPoints(Args...) {}
  template <typename... Args>
  static void plotLine(Args...) {}
  template <typename... Args>
  static void plotLinesPoints(Args...) {}
  template <typename... Args>
  static void plotBandsLinesPoints(Args...) {}
};
#endif  // DCA_WITH_GNUPLOT

}  // util
}  // dca

#endif  // DCA_UTIL_PLOT_HPP
