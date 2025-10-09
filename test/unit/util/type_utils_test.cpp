// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// This file tests type_help.hpp
// It's primary use is to allow refactoring or reimplementing as metaprogramming gets easier
// with advancing C++ standards.
#include "dca/util/type_help.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"


#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu_complex.h"
#endif

TEST(TypeListTest, Sublist) {
  EXPECT_EQ(dca::util::IsComplex_t<std::complex<double>>::value, true);
  EXPECT_EQ(dca::util::IsComplex_t<double>::value, false);
#ifdef DCA_HAVE_GPU
  EXPECT_EQ(dca::util::IsComplex_t<cuComplex>::value, true);
  EXPECT_EQ(dca::util::IsComplex_t<cuDoubleComplex>::value, true);
#endif
}  
