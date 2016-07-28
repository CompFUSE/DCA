// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests print_type.hpp.

#include "dca/util/print_type.hpp"
#include "gtest/gtest.h"

TEST(PrintTypeTest, FundamentalTypes) {
  EXPECT_EQ("void", dca::util::Type<void>::print());
  EXPECT_EQ("bool", dca::util::Type<bool>::print());
  EXPECT_EQ("wchar_t", dca::util::Type<wchar_t>::print());
  EXPECT_EQ("int", dca::util::Type<int>::print());
  EXPECT_EQ("unsigned long", dca::util::Type<unsigned long>::print());
  EXPECT_EQ("float", dca::util::Type<float>::print());
}

class MySimpleClass {};

template <typename T>
class MyClassTemplate {};

namespace ns1 {
class MyClassInNamespace1 {};
}  // ns1

namespace ns2 {
class MyClassInNamespace2 {};
}  // ns2

TEST(PrintTypeTest, CustomTypes) {
  EXPECT_EQ("MySimpleClass", dca::util::Type<MySimpleClass>::print());

  EXPECT_EQ("MyClassTemplate<double>", dca::util::Type<MyClassTemplate<double>>::print());

  EXPECT_EQ("ns1::MyClassInNamespace1", dca::util::Type<ns1::MyClassInNamespace1>::print());

  using namespace ns2;
  EXPECT_EQ("ns2::MyClassInNamespace2", dca::util::Type<MyClassInNamespace2>::print());
}
