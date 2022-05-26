// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests integration between type_list and function

#include "dca/util/type_list.hpp"
#include "function_testing.hpp"
#include "gtest/gtest.h"
#include <sstream>

TEST(TypeUtils, FunctionTypes) {
  using Domain16v = dca::testing::Domain16v;
  std::ostringstream domain_16v_print_out;
  dca::util::print_type<Domain16v>::print(domain_16v_print_out);
  std::cout << "Printing typelist test_domain_16v \n" << domain_16v_print_out.str() << '\n';
  // Be careful with formatting this file, there are tabs in the literal reference output below.
  constexpr std::string_view expected{R"(	dca::func::dmn_0<dca::func::dmn<1, double> >
	dca::func::dmn_0<dca::func::dmn<2, double> >
	dca::func::dmn_0<dca::func::dmn<4, double> >
	dca::func::dmn_0<dca::func::dmn<8, double> >
	dca::func::dmn_0<dca::func::dmn<1, double> >
	dca::func::dmn_0<dca::func::dmn<2, double> >
	dca::func::dmn_0<dca::func::dmn<4, double> >
	dca::func::dmn_0<dca::func::dmn<8, double> >
	dca::func::dmn_0<dca::func::dmn<4, double> >
	dca::func::dmn_0<dca::func::dmn<8, double> >
)"};
  EXPECT_EQ(expected, domain_16v_print_out.str());
}
