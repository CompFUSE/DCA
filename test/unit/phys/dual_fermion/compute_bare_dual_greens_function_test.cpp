// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests compute_bare_dual_greens_function.hpp.

#include "dca/phys/dual_fermion/compute_bare_dual_greens_function.hpp"

#include <complex>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"

using namespace dca;

class ComputeBareDualGreensFunctionTest : public ::testing::Test {
protected:
  using SpinDmn = func::dmn_0<func::dmn<2, int>>;
  using OrbitalDmn = func::dmn_0<func::dmn<1, int>>;
  using OrbitalSpinDmn = func::dmn_variadic<OrbitalDmn, SpinDmn>;
  using KClusterDmn = func::dmn_0<func::dmn<2, std::vector<double>>>;
  using KSuperlatticeDmn = func::dmn_0<func::dmn<4, std::vector<double>>>;
  using MatsubaraFreqDmn = func::dmn_0<func::dmn<2, double>>;

  ComputeBareDualGreensFunctionTest() : concurrency_(0, nullptr) {}

  const parallel::NoConcurrency concurrency_;
};

TEST_F(ComputeBareDualGreensFunctionTest, TwoClusterSites) {
  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KClusterDmn>>>
      eps_cg;

  func::function<std::complex<double>, func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn>>
      eps_tilde;

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KClusterDmn>,
                                    func::dmn_0<MatsubaraFreqDmn>>>
      Delta;

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KClusterDmn>,
                                    func::dmn_0<MatsubaraFreqDmn>>>
      G;

  func::function<std::complex<double>, func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn,
                                                          func::dmn_0<MatsubaraFreqDmn>>>

      G0_tilde;

  phys::df::computeBareDualGreensFunction(eps_cg, eps_tilde, Delta, G, concurrency_, G0_tilde);
}
