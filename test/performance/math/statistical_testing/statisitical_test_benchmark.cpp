// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the speed of the statistical test for a large number of d.o.f. Compare the speed
// of the computation performed in the eigenbases, were the normalized samples are computed, with
// the LU decomposition.

#include "dca/math/statistical_testing/statistical_testing.hpp"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/linalg/blas/blas3.hpp"
#include "dca/math/random/random.hpp"
#include "dca/profiling/events/time.hpp"

constexpr int k = 1000;
using Domain = dca::func::dmn_0<dca::func::dmn<k>>;
using Function = dca::func::function<double, Domain>;
using Covariance = dca::func::function<double, dca::func::dmn_variadic<Domain, Domain>>;
// Compute cov <- cov * cov^t
void makePositive(Covariance& cov);

int main() {
  dca::func::function<double, Domain> f(""), f0("");
  dca::func::function<double, dca::func::dmn_variadic<Domain, Domain>> cov("");

  // Create random input.
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1);
  for (std::size_t i = 0; i < f.size(); ++i)
    f(i) = rng();
  for (std::size_t i = 0; i < cov.size(); ++i)
    cov(i) = rng();
  makePositive(cov);

  // 1st timed section. Eigenbase computation.
  {
    dca::profiling::WallTime start;
    dca::math::StatisticalTesting test(f, f0, cov);
    double pval = test.computePValue(true, 1, /*allow fast*/ false);
    dca::profiling::WallTime end;
    dca::profiling::Duration time(end, start);
    std::cout << "allow fast = false" << std::endl;
    std::cout << "Final dof: " << test.get_dof() << " pvalue: " << pval << std::endl;
    std::cout << "Time taken: " << time.sec + 1e-6 * time.usec << std::endl;
  }
  // 2nd timed section. LU decomposition.
  {
    dca::profiling::WallTime start;
    dca::math::StatisticalTesting test(f, f0, cov);
    double pval = test.computePValue(true, 1, /*allow fast*/ true);
    dca::profiling::WallTime end;
    dca::profiling::Duration time(end, start);
    std::cout << "allow fast = true" << std::endl;
    std::cout << "Final dof: " << test.get_dof() << " pvalue: " << pval << std::endl;
    std::cout << "Time taken: " << time.sec + 1e-6 * time.usec << std::endl;
  }
}

void makePositive(Covariance& cov) {
  Covariance result("");
  int n = cov[0];

  dca::linalg::blas::gemm("N", "T", n, n, n, 1., &cov(0), n, &cov(0), n, 0., &result(0), n);
  cov = result;
}
