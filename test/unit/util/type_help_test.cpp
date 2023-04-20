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
#include "dca/platform/dca_gpu.h"
#include "dca/util/type_help.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"


#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu_complex.h"
#endif

TEST(TypeHelpTest, IsComplex) {
  EXPECT_EQ(dca::util::IsComplex_t<std::complex<double>>::value, true);
  EXPECT_EQ(dca::util::IsComplex_t<double>::value, false);
#ifdef DCA_HAVE_GPU
  EXPECT_EQ(dca::util::IsComplex_t<cuComplex>::value, false);
  EXPECT_EQ(dca::util::IsComplex_t<cuDoubleComplex>::value, false);
#endif
}  

#ifdef DCA_HAVE_GPU
TEST(TypeHelpTest, IsCudaComplex) {
  EXPECT_TRUE(dca::util::IsCudaComplex_t<cuComplex>::value);
  EXPECT_TRUE(dca::util::IsCudaComplex_t<float2>::value);
  EXPECT_TRUE(dca::util::IsCudaComplex_t<cuDoubleComplex>::value);
  EXPECT_TRUE(dca::util::IsCudaComplex_t<double2>::value);
  EXPECT_FALSE(dca::util::IsCudaComplex_t<std::complex<double>>::value);
}
#endif

TEST(TypeHelpTest, ComplexAlias) {
  bool is_same = std::is_same_v<dca::util::ComplexAlias<float>, std::complex<float>>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<std::complex<float>>, std::complex<float>>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<double>, std::complex<double>>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<std::complex<double>>, std::complex<double>>;
  EXPECT_TRUE(is_same);
#ifdef DCA_HAVE_GPU
  is_same = std::is_same_v<dca::util::ComplexAlias<float2>, float2>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<float2*>, float2*>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<float2**>, float2**>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<double2>, double2>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<double2*>, double2*>;
  EXPECT_TRUE(is_same);
  is_same = std::is_same_v<dca::util::ComplexAlias<double2**>, double2**>;
  EXPECT_TRUE(is_same);
#endif
}  

TEST(TypeHelpTest, castHostType) {
  std::complex<double> a_complex{1.0, 2.0};
  std::complex<double>* p_a_complex = &a_complex;
  std::complex<double>** pp_a_complex = &p_a_complex;

  const std::complex<double>* c_p_a_complex = &a_complex;
  const std::complex<double>** c_pp_a_complex = &c_p_a_complex;
#ifdef DCA_HAVE_GPU
  auto cast_pp_a_complex = dca::util::castHostType(pp_a_complex);
  auto cast_c_pp_a_complex = dca::util::castHostType(c_pp_a_complex);
  EXPECT_TRUE(**cast_pp_a_complex == **cast_c_pp_a_complex);
#endif
}  


TEST(TypeHelpTest, LinalgConstants) {
  std::complex<double> one{1.0,0.0};
  std::complex<double> made_one;
  dca::util::makeOne(made_one);
  EXPECT_EQ(made_one, one);
  double d_one{1.0};
  double d_made_one;
  dca::util::makeOne(d_made_one);
  EXPECT_EQ(d_made_one, d_one);
#ifdef DCA_HAVE_GPU
  double2 d2_one{1.0,0.0};
  double2 made_d2_one;
  dca::util::makeOne(made_d2_one);
  EXPECT_EQ(made_d2_one.x, d2_one.x);
  EXPECT_EQ(made_d2_one.y, d2_one.y);
  double2 d2_zero{0.0,0.0};
  double2 made_d2_zero;
  dca::util::makeZero(made_d2_zero);
  EXPECT_EQ(made_d2_zero.x, d2_zero.x);
  EXPECT_EQ(made_d2_zero.y, d2_zero.y);
  #endif
}
