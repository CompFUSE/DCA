// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// This file tests phase

#include "dca/math/util/phase.hpp"
#include <complex>
#include <utility>
#include "gtest/gtest.h"

TEST(Phase, metaprogramming) {
  using PhaseTrue = dca::math::Phase<std::complex<double>>;
  if constexpr (!dca::math::IsPhase<PhaseTrue>::value)
    FAIL() << "a Phase is a phase";
}
