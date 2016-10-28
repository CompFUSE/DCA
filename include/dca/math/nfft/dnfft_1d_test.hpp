// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a class for testing dnfft_1d.hpp.
//
// TODO: Move to test dir and create proper (unit) test.

#ifndef DCA_MATH_NFFT_DNFFT_1D_TEST_HPP
#define DCA_MATH_NFFT_DNFFT_1D_TEST_HPP

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/nfft/dnfft_1d.hpp"
#include "dca/util/plot.hpp"

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <typename scalartype, typename t_dmn_t, typename w_dmn_t>
class dnfft_1D_test {
  typedef func::dmn_0<func::dmn<1, int>> p_dmn_t;

public:
  static void execute() {
    srand(55423543);

    double begin = t_dmn_t::get_elements().front();
    double delta = (t_dmn_t::get_elements().back() - t_dmn_t::get_elements().front());

    std::vector<scalartype> t(0);
    std::vector<scalartype> f(0);
    for (int l = 0; l < 1.e6; l++) {
      t.push_back(t_dmn_t::get_elements()[0] + double(rand()) / double(RAND_MAX) * delta);
      f.push_back(exp(-2. * M_PI / delta * (t[l] - begin)));
    }

    func::function<std::complex<scalartype>, w_dmn_t> f_w_1("f_w_1");
    func::function<std::complex<scalartype>, w_dmn_t> f_w_2("f_w_2");
    func::function<std::complex<scalartype>, w_dmn_t> error("error");

    compute_f_w(t, f, f_w_1);
    compute_f_w_dnfft(t, f, f_w_2);

    for (int i = 0; i < w_dmn_t::dmn_size(); i++)
      error(i) = log10(abs(f_w_1(i) - f_w_2(i)));

    util::Plot::plotLinesPoints(error);

    for (int i = 0; i < w_dmn_t::dmn_size(); i++)
      error(i) = log10(abs(f_w_1(i) - f_w_2(i)) / (abs(f_w_1(i)) + 1.e-6));

    util::Plot::plotLinesPoints(error);

    std::cout << "\n\n\n\t STOP !!!!  \n\n\n" << std::endl;

    throw std::logic_error(__FUNCTION__);
  }

private:
  static void compute_f_w(std::vector<scalartype>& t, std::vector<scalartype>& f,
                          func::function<std::complex<scalartype>, w_dmn_t>& f_w) {
    std::cout << __FUNCTION__ << std::endl;

    std::complex<scalartype> I(0, 1);

    {
      clock_t t0 = clock();

      for (int j = 0; j < t.size(); j++)
        for (int i = 0; i < w_dmn_t::dmn_size(); i++) {
          scalartype t_val = t[j];
          scalartype w_val = w_dmn_t::get_elements()[i];

          f_w(i) += f[j] * std::exp(I * t_val * w_val);
        }
      clock_t t1 = clock();

      std::cout << "\n\n\t time : " << double(t1 - t0) / double(CLOCKS_PER_SEC) << "\n";
    }

    util::Plot::plotLinesPoints(f_w);
  }

  static void compute_f_w_dnfft(std::vector<scalartype>& t, std::vector<scalartype>& f,
                                func::function<std::complex<scalartype>, w_dmn_t>& f_w) {
    std::cout << __FUNCTION__ << std::endl;

    func::function<std::complex<scalartype>, func::dmn_variadic<w_dmn_t, p_dmn_t>> f_w_tmp;

    dnfft_1D<scalartype, w_dmn_t, p_dmn_t> nfft_obj;

    std::vector<scalartype> t_vals(0);
    for (int j = 0; j < t.size(); j++)
      t_vals.push_back((t[j] - t_dmn_t::get_elements().front()) /
                           (t_dmn_t::get_elements().back() - t_dmn_t::get_elements().front()) -
                       0.5);

    {
      clock_t t0 = clock();
      for (int j = 0; j < t.size(); j++) {
        nfft_obj.accumulate_at(0, t_vals[j], f[j]);
      }

      clock_t t1 = clock();

      std::cout << "\n\n\t time : " << double(t1 - t0) / double(CLOCKS_PER_SEC) << "\n";
    }

    nfft_obj.finalize(f_w_tmp);

    for (int i = 0; i < w_dmn_t::dmn_size(); i++)
      f_w(i) = f_w_tmp(i);

    util::Plot::plotLinesPoints(f_w);
  }
};

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DNFFT_1D_TEST_HPP
