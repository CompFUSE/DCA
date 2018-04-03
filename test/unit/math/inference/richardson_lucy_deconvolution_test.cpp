// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests richardson_lucy_deconvolution.hpp

#include "dca/math/inference/richardson_lucy_deconvolution.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"

namespace testing {

class DeconvolutionParameters {
  struct Concurrency {
    int get_id() {
      return 0;
    }
    int first() {
      return 0;
    }
  };

public:
  using concurrency_type = Concurrency;

  DeconvolutionParameters(int iterations = 16, double tolerance = 1.e-2)
      : iterations_(iterations), tolerance_(tolerance) {}

  int get_deconvolution_iterations() const {
    return iterations_;
  }
  double get_deconvolution_tolerance() const {
    return tolerance_;
  }
  concurrency_type& get_concurrency() {
    return concurrency_;
  }

private:
  int iterations_;
  double tolerance_;

  concurrency_type concurrency_;
};

}  // testing

TEST(RichardsonLucyDeconvolutionTest, Constructor) {
  using DeconvolutionDmn = dca::func::dmn_0<dca::func::dmn<4, int>>;
  using OtherDmn = dca::func::dmn_0<dca::func::dmn<1, int>>;

  testing::DeconvolutionParameters params(1, 1.e-3);

  dca::math::inference::RichardsonLucyDeconvolution<testing::DeconvolutionParameters,
                                                    DeconvolutionDmn, OtherDmn>
      deconvolution(params);
}
