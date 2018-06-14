// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests reduced_domain.hpp.

#include "dca/function/domains/reduced_domain.hpp"
#include <vector>
#include "gtest/gtest.h"
#include "dca/function/domains/dmn.hpp"

TEST(ReducedDomainTest, InitializeWithIndices) {
  // Default name
  using BaseDmn1 = dca::func::dmn<2, int>;
  BaseDmn1::set_elements(std::vector<int>{0, 10});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  const std::vector<std::size_t> indices1{1};
  ReducedDmn1::initialize(indices1);

  const std::vector<int> check1{10};
  EXPECT_EQ(indices1.size(), ReducedDmn1::get_size());
  EXPECT_EQ(check1, ReducedDmn1::get_elements());
  EXPECT_EQ("dca::func::dmn<2, int>_reduced", ReducedDmn1::get_name());

  // Specify name
  using BaseDmn2 = dca::func::dmn<3, int>;
  BaseDmn2::set_elements(std::vector<int>{-1, -2, -3});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  const std::vector<std::size_t> indices2{1, 2};
  ReducedDmn2::initialize(indices2, "ReducedDomain-with-custom-name");

  const std::vector<int> check2{-2, -3};
  EXPECT_EQ(indices2.size(), ReducedDmn2::get_size());
  EXPECT_EQ(check2, ReducedDmn2::get_elements());
  EXPECT_EQ("ReducedDomain-with-custom-name", ReducedDmn2::get_name());

  // Empty indices vector
  using BaseDmn3 = dca::func::dmn<4, int>;
  BaseDmn3::set_elements(std::vector<int>{42, 314, -11, 101});

  using ReducedDmn3 = dca::func::ReducedDomain<BaseDmn3>;
  const std::vector<std::size_t> indices3;
  ReducedDmn3::initialize(indices3);

  EXPECT_EQ(indices3.size(), ReducedDmn3::get_size());
  EXPECT_EQ(std::vector<int>{}, ReducedDmn3::get_elements());

  // All base domain elements
  using BaseDmn4 = dca::func::dmn<5, int>;
  const std::vector<int> elements4{1, 3, 5, 7, 9};
  BaseDmn4::set_elements(elements4);

  using ReducedDmn4 = dca::func::ReducedDomain<BaseDmn4>;
  const std::vector<std::size_t> indices4{0, 1, 2, 3, 4};
  ReducedDmn4::initialize(indices4);

  EXPECT_EQ(indices4.size(), ReducedDmn4::get_size());
  EXPECT_EQ(elements4, ReducedDmn4::get_elements());
}

TEST(ReducedDomainTest, InitializeWithIndicesExceptions) {
  // Index >= base-domain-size
  using BaseDmn1 = dca::func::dmn<2, double>;
  BaseDmn1::set_elements(std::vector<double>{0., 10.});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  EXPECT_THROW(ReducedDmn1::initialize(std::vector<std::size_t>{0, 2}), std::out_of_range);

  // Vector of indices isn't sorted.
  using BaseDmn2 = dca::func::dmn<3, double>;
  BaseDmn2::set_elements(std::vector<double>{0., 10., 20.});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  EXPECT_THROW(ReducedDmn2::initialize(std::vector<std::size_t>{2, 0}), std::logic_error);

  // Vector of indices contains duplicates.
  using BaseDmn3 = dca::func::dmn<4, double>;
  BaseDmn3::set_elements(std::vector<double>{0., 10., 20., 40.});

  using ReducedDmn3 = dca::func::ReducedDomain<BaseDmn3>;
  EXPECT_THROW(ReducedDmn3::initialize(std::vector<std::size_t>{1, 1, 2}), std::logic_error);
}

TEST(ReducedDomainTest, ReinitializationException) {
  using BaseDmn1 = dca::func::dmn<5, double>;
  BaseDmn1::set_elements(std::vector<double>{0., 10., 20., 30., 40.});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  ReducedDmn1::initialize(std::vector<std::size_t>{0, 1, 2});
  EXPECT_THROW(ReducedDmn1::initialize(std::vector<std::size_t>{0, 1}), std::logic_error);
  EXPECT_THROW(ReducedDmn1::initialize(2, 4), std::logic_error);
}

TEST(ReducedDomainTest, InitializeWithRange) {
  using BaseDmn1 = dca::func::dmn<5, long>;
  BaseDmn1::set_elements(std::vector<long>{0, 10, 20, 30, 40});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  ReducedDmn1::initialize(1, 4, "ReducedDomain-intialized-with-range");

  const std::vector<long> check1{10, 20, 30};
  EXPECT_EQ(3, ReducedDmn1::get_size());
  EXPECT_EQ(check1, ReducedDmn1::get_elements());
  EXPECT_EQ("ReducedDomain-intialized-with-range", ReducedDmn1::get_name());

  // first = 0
  using BaseDmn2 = dca::func::dmn<6, long>;
  BaseDmn2::set_elements(std::vector<long>{0, 10, 20, 30, 40, 50});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  ReducedDmn2::initialize(0, 4);

  const std::vector<long> check2{0, 10, 20, 30};
  EXPECT_EQ(4, ReducedDmn2::get_size());
  EXPECT_EQ(check2, ReducedDmn2::get_elements());
  EXPECT_EQ("dca::func::dmn<6, long>_reduced", ReducedDmn2::get_name());

  // last = base-domain-size
  using BaseDmn3 = dca::func::dmn<7, long>;
  BaseDmn3::set_elements(std::vector<long>{0, 10, 20, 30, 40, 50, 60});

  using ReducedDmn3 = dca::func::ReducedDomain<BaseDmn3>;
  ReducedDmn3::initialize(5, 7);

  const std::vector<long> check3{50, 60};
  EXPECT_EQ(2, ReducedDmn3::get_size());
  EXPECT_EQ(check3, ReducedDmn3::get_elements());
  EXPECT_EQ("dca::func::dmn<7, long>_reduced", ReducedDmn3::get_name());

  // first = last
  using BaseDmn4 = dca::func::dmn<8, long>;
  BaseDmn4::set_elements(std::vector<long>{0, 10, 20, 30, 40, 50, 60, 70});

  using ReducedDmn4 = dca::func::ReducedDomain<BaseDmn4>;
  ReducedDmn4::initialize(3, 3);

  EXPECT_EQ(0, ReducedDmn4::get_size());
  EXPECT_EQ(std::vector<long>{}, ReducedDmn4::get_elements());
  EXPECT_EQ("dca::func::dmn<8, long>_reduced", ReducedDmn4::get_name());
}

TEST(ReducedDomainTest, InitializeWithRangeExceptions) {
  // first >= base-domain-size
  using BaseDmn1 = dca::func::dmn<6, double>;
  BaseDmn1::set_elements(std::vector<double>{0., 10., 20., 30., 40., 50.});

  using ReducedDmn1 = dca::func::ReducedDomain<BaseDmn1>;
  EXPECT_THROW(ReducedDmn1::initialize(6, 6), std::out_of_range);

  // last > base-domain-size
  using BaseDmn2 = dca::func::dmn<7, double>;
  BaseDmn2::set_elements(std::vector<double>{0., 10., 20., 30., 40., 50., 60.});

  using ReducedDmn2 = dca::func::ReducedDomain<BaseDmn2>;
  EXPECT_THROW(ReducedDmn2::initialize(1, 8), std::out_of_range);

  // first > last
  using BaseDmn3 = dca::func::dmn<8, double>;
  BaseDmn3::set_elements(std::vector<double>{0., 10., 20., 30., 40., 50., 60., 70.});

  using ReducedDmn3 = dca::func::ReducedDomain<BaseDmn3>;
  EXPECT_THROW(ReducedDmn3::initialize(2, 1), std::logic_error);
}
