// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests reduced_domain.hpp.

#include "comp_library/function_library/domains/special_domains/reduced_domain.hpp"
#include <vector>
#include "gtest/gtest.h"
#include "comp_library/function_library/domains/special_domains/dmn.h"

TEST(ReducedDomainTest, InitializeWithIndices) {
  // Default name
  using BaseDmn1 = dmn<2, int>;
  BaseDmn1::set_elements(std::vector<int>{0, 10});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  const std::vector<std::size_t> indices1{1};
  ReducedDmn1::initialize(indices1);

  const std::vector<int> check1{10};
  EXPECT_EQ(indices1.size(), ReducedDmn1::get_size());
  EXPECT_EQ(check1, ReducedDmn1::get_elements());
  EXPECT_EQ("dmn<2, int>_reduced", ReducedDmn1::get_name());

  // Specify name
  using BaseDmn2 = dmn<3, int>;
  BaseDmn2::set_elements(std::vector<int>{-1, -2, -3});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  const std::vector<std::size_t> indices2{1, 2};
  ReducedDmn2::initialize(indices2, "ReducedDomain-with-custom-name");

  const std::vector<int> check2{-2, -3};
  EXPECT_EQ(indices2.size(), ReducedDmn2::get_size());
  EXPECT_EQ(check2, ReducedDmn2::get_elements());
  EXPECT_EQ("ReducedDomain-with-custom-name", ReducedDmn2::get_name());

  // Empty indices vector
  using BaseDmn3 = dmn<4, int>;
  BaseDmn3::set_elements(std::vector<int>{42, 314, -11, 101});

  using ReducedDmn3 = dca::func::ReducedDomain<BaseDmn3>;
  const std::vector<std::size_t> indices3;
  ReducedDmn3::initialize(indices3);

  EXPECT_EQ(indices3.size(), ReducedDmn3::get_size());
  EXPECT_EQ(std::vector<int>{}, ReducedDmn3::get_elements());

  // All base domain elements
  using BaseDmn4 = dmn<5, int>;
  const std::vector<int> elements4{1, 3, 5, 7, 9};
  BaseDmn4::set_elements(elements4);

  using ReducedDmn4 = dca::func::ReducedDomain<BaseDmn4>;
  const std::vector<std::size_t> indices4{0, 1, 2, 3, 4};
  ReducedDmn4::initialize(indices4);

  EXPECT_EQ(indices4.size(), ReducedDmn4::get_size());
  EXPECT_EQ(elements4, ReducedDmn4::get_elements());
}

TEST(ReducedDomainTest, InitializeWithIndicesExceptions) {
  // More indicies than elements in the base domain
  using BaseDmn1 = dmn<2, double>;
  BaseDmn1::set_elements(std::vector<double>{0., 10.});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  EXPECT_THROW(ReducedDmn1::initialize(std::vector<std::size_t>{0, 1, 2}), std::out_of_range);

  // Index >= base-domain-size
  using BaseDmn2 = dmn<3, double>;
  BaseDmn2::set_elements(std::vector<double>{0., 10., 20.});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  EXPECT_THROW(ReducedDmn2::initialize(std::vector<std::size_t>{1, 3}), std::out_of_range);
}

TEST(ReducedDomainTest, ReinitializationException) {
  using BaseDmn1 = dmn<4, double>;
  BaseDmn1::set_elements(std::vector<double>{0., 10., 20., 30.});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  ReducedDmn1::initialize(std::vector<std::size_t>{0, 1, 2});
  EXPECT_THROW(ReducedDmn1::initialize(std::vector<std::size_t>{0, 1}), std::logic_error);
  EXPECT_THROW(ReducedDmn1::initialize(2, false), std::logic_error);
}

TEST(ReducedDomainTest, InitializeWithRangeFromFirst) {
  using BaseDmn1 = dmn<4, long>;
  BaseDmn1::set_elements(std::vector<long>{0, 10, 20, 30});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  ReducedDmn1::initialize(3, false, "ReducedDomain-intialized-with-range");

  const std::vector<long> check1{0, 10, 20};
  EXPECT_EQ(3, ReducedDmn1::get_size());
  EXPECT_EQ(check1, ReducedDmn1::get_elements());
  EXPECT_EQ("ReducedDomain-intialized-with-range", ReducedDmn1::get_name());

  // Zero
  using BaseDmn2 = dmn<5, long>;
  BaseDmn2::set_elements(std::vector<long>{0, 10, 20, 30, 40});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  ReducedDmn2::initialize(0, false);

  EXPECT_EQ(0, ReducedDmn2::get_size());
  EXPECT_EQ(std::vector<long>{}, ReducedDmn2::get_elements());
  EXPECT_EQ("dmn<5, long>_reduced", ReducedDmn2::get_name());
}

TEST(ReducedDomainTest, InitializeWithRangeFromFirstExeceptions) {
  // Range to copy goes beyond end of base domain
  using BaseDmn1 = dmn<4, float>;
  BaseDmn1::set_elements(std::vector<float>{-1., -2., -3., -4.});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  EXPECT_THROW(ReducedDmn1::initialize(5, false), std::out_of_range);
}

TEST(ReducedDomainTest, InitializeWithRangeFromMiddle) {
  // Even
  using BaseDmn1 = dmn<4, unsigned int>;
  BaseDmn1::set_elements(std::vector<unsigned int>{0, 10, 20, 30});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  ReducedDmn1::initialize(2, true);

  const std::vector<unsigned int> check1{20, 30};
  EXPECT_EQ(2, ReducedDmn1::get_size());
  EXPECT_EQ(check1, ReducedDmn1::get_elements());
  EXPECT_EQ("dmn<4, unsigned int>_reduced", ReducedDmn1::get_name());

  // Odd
  using BaseDmn2 = dmn<5, unsigned int>;
  BaseDmn2::set_elements(std::vector<unsigned int>{0, 10, 20, 30, 40});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  ReducedDmn2::initialize(3, true);

  const std::vector<unsigned int> check2{20, 30, 40};
  EXPECT_EQ(3, ReducedDmn2::get_size());
  EXPECT_EQ(check2, ReducedDmn2::get_elements());

  // Zero
  using BaseDmn3 = dmn<6, unsigned int>;
  BaseDmn3::set_elements(std::vector<unsigned int>{0, 10, 20, 30, 40, 50});

  using ReducedDmn3 = dca::func::ReducedDomain<BaseDmn3>;
  ReducedDmn3::initialize(0, true);

  EXPECT_EQ(0, ReducedDmn3::get_size());
  EXPECT_EQ(std::vector<unsigned int>{}, ReducedDmn3::get_elements());
}

TEST(ReducedDomainTest, InitializeWithRangeFromMiddleExceptions) {
  // Range to copy goes beyond end of base domain
  // Even
  using BaseDmn1 = dmn<6, float>;
  BaseDmn1::set_elements(std::vector<float>{-1., -2., -3., -4., -5., -6.});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  EXPECT_THROW(ReducedDmn1::initialize(4, true), std::out_of_range);

  // Odd
  using BaseDmn2 = dmn<7, float>;
  BaseDmn2::set_elements(std::vector<float>{-1., -2., -3., -4., -5., -6., -7.});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  EXPECT_THROW(ReducedDmn2::initialize(5, true), std::out_of_range);
}
