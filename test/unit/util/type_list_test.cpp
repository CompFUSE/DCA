// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the new methods in type_lists.hpp

#include "dca/util/type_utils.hpp"

#include "gtest/gtest.h"

TEST(TypeListTest, Sublist) {
  using List = dca::util::Typelist<int, float, double, char>;
  using Sublist = dca::util::Sublist<2, List>;

  EXPECT_EQ(4, dca::util::Length<List>::value);
  EXPECT_EQ(2, dca::util::Length<Sublist>::value);

  constexpr bool elements_eq = std::is_same<dca::util::Typelist<int, float>, Sublist>::value;
  EXPECT_TRUE(elements_eq);
}

TEST(TypeListTest, lists) {
    
    using test_list = dca::util::Typelist<int*, double, float&, int, int*, int, int, const float*,
                                          const int, long long, int&>;
    using test_list2 = dca::util::Typelist<float, int, int*, double*, int, int, int*, double, int,
                                           float, const int, long long, int&>;
    std::cout << "Testing Typelist Length " << dca::util::mp_size<test_list>::value << std::endl;
    std::cout << "Testing Typelist Length " << dca::util::Length<test_list>::value << std::endl;

    std::cout << "Testing Typelist NumberOf " << dca::util::mp_count<test_list, int>::value
              << std::endl;
    std::cout << "Testing Typelist NumberOf " << dca::util::NumberOf<test_list, int>::value
              << std::endl;

    std::cout << "Testing Typelist IndexOf " << dca::util::mp_index_of<long long, test_list>::value
              << std::endl;
    std::cout << "Testing Typelist IndexOf " << dca::util::IndexOf<long long, test_list>::value
              << std::endl;

    std::cout << "Testing Typelist TypeAt "
              << dca::util::type_name<dca::util::mp_element<9, test_list>::type>().c_str()
              << std::endl;
    std::cout << "Testing Typelist TypeAt "
              << dca::util::type_name<dca::util::TypeAt<9, test_list>::type>().c_str() << std::endl;

    std::cout << "Testing Typelist Append "
              << dca::util::mp_size<dca::util::mp_append<test_list, test_list2>::type>::value
              << std::endl;
    std::cout << "Testing Typelist Append "
              << dca::util::Length<dca::util::Append<test_list, test_list2>::type>::value
              << std::endl;
    std::cout << "Testing Typelist Append/Index "
              << dca::util::mp_index_of<const float*, dca::util::mp_append<test_list, test_list2>::type>::value
              << std::endl;
    std::cout << "Testing Typelist Append/Index "
              << dca::util::IndexOf<const float*, dca::util::Append<test_list, test_list2>::type>::value
              << std::endl;

    std::cout << "Testing Typelist Prepend "
              << dca::util::mp_size<dca::util::mp_prepend<test_list, test_list2>::type>::value
              << std::endl;
    std::cout << "Testing Typelist Prepend "
              << dca::util::Length<dca::util::Prepend<test_list, test_list2>::type>::value
              << std::endl;
    std::cout << "Testing Typelist Prepend/Index "
              << dca::util::mp_index_of<const float*,
                                        dca::util::mp_prepend<test_list, test_list2>::type>::value
              << std::endl;
    std::cout << "Testing Typelist Prepend/Index "
              << dca::util::IndexOf<const float*, dca::util::Prepend<test_list, test_list2>::type>::value
              << std::endl;
    std::cout << std::endl;

    std::cout << "Printing typelist test_list \n";
    dca::util::print_type<test_list>::print(std::cout);
    std::cout << std::endl;
}
