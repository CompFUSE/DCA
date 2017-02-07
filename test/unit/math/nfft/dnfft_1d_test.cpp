// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file tests the delayed non uniform fast fourier transform implemented in dnfft_1d.hpp.

#include "gtest/gtest.h"
#include "dca/math/nfft/dnfft_1d.hpp"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/random/random.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
//#include "dca/util/plot.hpp"

using dca::func::function;
using Tdmn = dca::func::dmn_0<dca::phys::domains::time_domain>;
using Wdmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
using ParsDmn = dca::func::dmn_0<dca::func::dmn<1, int>>;

void computeF_w_dnfft(const std::vector<double>& t, const std::vector<double>& f,
                      function<std::complex<double>, Wdmn>& f_w);
void computeF_w_direct(const std::vector<double>& t, const std::vector<double>& f,
                       function<std::complex<double>, Wdmn>& f_w);

TEST(Dnfft, TransformTest) {
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);

  const double beta = 10;
  const int n_t = 100;
  const int n_w = n_t - 2;
  const int n_samples = 1.e4;

  dca::phys::domains::time_domain::initialize(beta, n_t, 1e-10);
  dca::phys::domains::frequency_domain::initialize(beta, n_w);

  const double begin = Tdmn::get_elements().front();
  const double delta = (Tdmn::get_elements().back() - Tdmn::get_elements().front());

  std::vector<double> t(0);
  std::vector<double> f(0);
  for (int l = 0; l < n_samples; l++) {
    t.push_back(begin + rng() * delta);
    f.push_back(exp(-2. * M_PI / delta * (t[l] - begin)));
  }

  function<std::complex<double>, Wdmn> f_w_1("f_w_1");
  function<std::complex<double>, Wdmn> f_w_2("f_w_2");

  computeF_w_direct(t, f, f_w_1);
  computeF_w_dnfft(t, f, f_w_2);

  double max_error = 0;
  double avg_error = 0;

  for (int i = 0; i < Wdmn::dmn_size(); i++) {
    const double err = abs(f_w_1(i) - f_w_2(i)) / (abs(f_w_1(i)) + 1.e-6);
    max_error = std::max(max_error, err);
    avg_error += err;
  }
  avg_error /= Wdmn::dmn_size();

  EXPECT_GE(1e-8, avg_error);
  EXPECT_GE(1e-6, max_error);

//util::Plot::plotLinesPoints(error);
}

void computeF_w_dnfft(const std::vector<double>& t, const std::vector<double>& f,
                      function<std::complex<double>, Wdmn>& f_w) {

  function<std::complex<double>, dca::func::dmn_variadic<Wdmn, ParsDmn>> f_w_tmp;

  dca::math::nfft::Dnfft1D<double, Wdmn, ParsDmn> nfft_obj;

  const double begin = Tdmn::get_elements().front();
  const double delta = (Tdmn::get_elements().back() - Tdmn::get_elements().front());

  for (int j = 0; j < t.size(); j++){
    const double scaled_t = (t[j] - begin) /  delta - 0.5;
    nfft_obj.accumulate(0, scaled_t, f[j]);
    }

  nfft_obj.finalize(f_w_tmp);

  for (int i = 0; i < Wdmn::dmn_size(); i++)
    f_w(i) = f_w_tmp(i);

  // util::Plot::plotLinesPoints(f_w);
}

void computeF_w_direct(const std::vector<double>& t, const std::vector<double>& f,
                       function<std::complex<double>, Wdmn>& f_w) {
  std::complex<double> I(0, 1);

  for (int j = 0; j < t.size(); j++)
    for (int i = 0; i < Wdmn::dmn_size(); i++) {
      const double t_val = t[j];
      const double w_val = Wdmn::get_elements()[i];

      f_w(i) += f[j] * std::exp(I * t_val * w_val);
    }

 // util::Plot::plotLinesPoints(f_w);
}
