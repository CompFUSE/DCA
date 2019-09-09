// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This file tests the function library.
//
// TODO: Move domains-only tests to separate files and rename this file to function_test.cpp and the
//       test cases to FunctionTest.

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

#include <cassert>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "dca/util/type_list.hpp"
#include "dca/util/type_utils.hpp"

using dca::func::dmn;
using dca::func::dmn_0;
using dca::func::dmn_variadic;

namespace dca {
namespace testing {
// dca::testing::

bool compare_to_file(const std::string& filename, const std::string& check) {
  // Open the file.
  std::ifstream known_result(filename);

  if (known_result.good()) {
    std::string contents;

    // Get the right size to reserve it.
    known_result.seekg(0, std::ios::end);
    contents.reserve(known_result.tellg());

    known_result.seekg(0, std::ios::beg);

    // Read contents into string directly.
    contents.assign((std::istreambuf_iterator<char>(known_result)), std::istreambuf_iterator<char>());

    return (contents == check);
  }

  else {
    std::ofstream new_result(filename);
    new_result << check.c_str();
    std::cout << "No baseline file exists, writing new one " << filename.c_str() << std::endl;
  }

  return false;
}

// A selection of domain types we can use for testing.
typedef dmn_0<dmn<1, double>> test_domain_0a;
typedef dmn_0<dmn<2, double>> test_domain_0b;
typedef dmn_0<dmn<4, double>> test_domain_0c;
typedef dmn_0<dmn<8, double>> test_domain_0d;
typedef dmn_variadic<test_domain_0d> test_domain_1d;
typedef dmn_variadic<test_domain_0a, test_domain_0b> test_domain_2a;
typedef dmn_variadic<test_domain_0c, test_domain_0d> test_domain_2c;
typedef dmn_variadic<test_domain_2a, test_domain_2c> test_domain_4a;
typedef dmn_variadic<test_domain_4a, test_domain_4a, test_domain_2c> test_domain_16;

typedef dmn_variadic<test_domain_0d> test_domain_0v;
typedef dmn_variadic<test_domain_4a, test_domain_4a, test_domain_2c> test_domain_16v;

dca::func::function<double, test_domain_0a> function_0a("test_domain_0a");
dca::func::function<double, test_domain_0b> function_0b("test_domain_0b");
dca::func::function<double, test_domain_0c> function_0c("test_domain_0c");
dca::func::function<double, test_domain_0d> function_0d("test_domain_0d");
dca::func::function<double, test_domain_1d> function_1d("test_domain_1d");
dca::func::function<double, test_domain_2a> function_2a("test_domain_2a");
dca::func::function<double, test_domain_4a> function_4a("test_domain_4a");
dca::func::function<double, test_domain_16> function_16("test_domain_16");

test_domain_0v dummy;
test_domain_16v dummy2;

template <typename T1>
struct function_test {};

template <int N, typename Dmn>
struct function_test<dca::func::function<double, dmn<N, Dmn>>> {
  typedef dca::func::function<double, dmn<N, Dmn>> fType;

  function_test(fType& func) : f(func) {}

  template <typename Arg>
  bool check_1(Arg /*arg*/) {
    /*
    std::cout << "Sub branch size " << std::endl;
    for (int i = 0; i < f.get_Nb_branch_domains(); i++) {
      std::cout << f.get_branch_size(i) << "\t";
    }
    std::cout << std::endl;

    std::cout << "Sub branch steps " << std::endl;
    for (int i = 0; i < f.get_Nb_branch_domains(); i++) {
      std::cout << f.get_branch_domain_steps()[i] << "\t";
    }
    std::cout << std::endl;
    */
    return true;
  }

  fType& f;
};

template <typename Domain>
struct function_test<dca::func::function<double, Domain>> {
  typedef dca::func::function<double, Domain> fType;
  typedef typename fType::this_scalar_type scalartype;
  typedef Domain domainType;
  // typedef typename Domain::this_type sub_type;

  // const int Ntypes = dca::util::Length<sub_type>::value;

  function_test(fType& func) : f(func) {}

  int signature() {
    return f.signature();
  }
  int size() {
    return f.size();
  }

  void fill_sequence() {
    int N = f.size();
    for (int i = 0; i < N; ++i) {
      f(i) = i;
      // if (i<1024) std::cout << i << ",";
    }
  }

  void check_sequence() {
    int N = f.size();
    for (int i = 0; i < N; ++i) {
      if (f(i) != i)
        throw(std::runtime_error("fault"));
    }
    // std::cout << "Ntypes " << Ntypes << " signature " << signature() << " size " << size()
    //           << std::endl;
  }

  template <typename Arg>
  bool check_1(Arg /*arg*/) {
    std::cout << "Sub branch size " << std::endl;
    for (int i = 0; i < f.get_domain().get_Nb_branch_domains(); i++) {
      std::cout << f.get_domain().get_branch_size(i) << "\t";
    }
    std::cout << std::endl;

    std::cout << "Sub branch steps " << std::endl;
    for (int i = 0; i < f.get_domain().get_Nb_branch_domains(); i++) {
      std::cout << f.get_domain().get_branch_domain_steps()[i] << "\t";
    }
    std::cout << std::endl;

    // IsSame<mp_append<test_domain_0a::this_type, test_domain_0b::this_type>::type, double>();
    // IsSame<test_domain_0a::this_type, test_domain_0b::this_type>();
    // IsSame<test_domain_0a, test_domain_0b>();
    // IsSame<test_domain_2a, test_domain_2c>();

    std::cout << std::endl;
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

    std::cout << "Printing typelist test_domain_16v \n";
    dca::util::print_type<test_domain_16v>::print(std::cout);

    std::cout << "Printing typelist test_list \n";
    dca::util::print_type<test_list>::print(std::cout);
    std::cout << std::endl;

    // std::cout << "\nTesting Typelist count "
    //           << dca::util::type_name<dca::util::TypeAt<2, test_list>>().c_str() << std::endl;
    // std::cout << "\nTesting Typelist count "
    //           << dca::util::type_name<dca::util::TypeAt<2, test_list>>().c_str() << std::endl;

    // typedef typename TypeAt<typename Domain::domain_typelist_0, 0>::Result dom_0;
    // std::cout << "Getting first subdomain "
    //           << "Type Id is " << typeid(dom_0).name() << std::endl;
    // dca::func::function<double, dom_0> sub_function;
    // function_test<decltype(sub_function)> sub_domain(sub_function);
    // sub_domain.check_1(1);

    return true;
  }

  template <typename... Args>
  scalartype expand(Args... /*args*/) {
    return scalartype(0);
  }

  template <typename... Args>
  bool check_value(Args... args) {
    // if (f(args...) == arg1 * offset<f, 1> + arg2 * offset<f, 2> +) {
    // }
    return f.operator()(args...) == f(args...);
    // return check_value(args...);
  }
  fType& f;
};

}  // namespace testing
}  // namespace dca

TEST(Function, TestDomain4a) {
  try {
    std::cout << "Leaf indexing \n";
    int index = 0;
    for (int i0 = 0; i0 < 1; ++i0) {
      for (int i1 = 0; i1 < 2; ++i1) {
        for (int i2 = 0; i2 < 4; ++i2) {
          for (int i3 = 0; i3 < 8; ++i3) {
            std::cout << i0 << "," << i1 << "," << i2 << "," << i3 << "\n";
            dca::testing::function_4a.operator()(i0, i1, i2, i3) = index++;
            // dca::testing::function_4a.operator()(i3,i2,i1,i0) = index; bad ordering
          }
        }
      }
    }
    std::cout << "Branch indexing \n";
    index = 0;
    for (int i0 = 0; i0 < 2; ++i0) {
      for (int i1 = 0; i1 < 32; ++i1) {
        std::cout << i0 << "," << i1 << "\n";
        dca::testing::function_4a.operator()(i0, i1) = index++;
      }
    }
  }
  catch (std::string& err) {
    // Check exception.
    std::cout << "Caught " << err.c_str() << std::endl;
    FAIL();
  }
}

TEST(Function, FillDomain) {
  ::testing::internal::TimeInMillis elapsed1(::testing::UnitTest::GetInstance()->elapsed_time());

  dca::testing::function_test<decltype(dca::testing::function_16)> ftest(dca::testing::function_16);

  std::cout << "Print fingerprint of function 16 \n";
  dca::testing::function_16.print_fingerprint(std::cout);

  std::cout << "ftest.fill_sequence \n";
  ftest.fill_sequence();

  std::cout << "ftest.check_sequence \n";
  ftest.check_sequence();

  std::cout << "ftest.check_1 \n";
  ftest.check_1(1);

  ::testing::internal::TimeInMillis elapsed2(::testing::UnitTest::GetInstance()->elapsed_time());

  std::cout << "Elapsed time is " << (elapsed2 - elapsed1) << std::endl;
}

TEST(Function, FingerPrint) {
  std::stringstream result;

  dca::testing::function_0a.print_fingerprint(result);
  dca::testing::function_0b.print_fingerprint(result);
  dca::testing::function_0c.print_fingerprint(result);
  dca::testing::function_0d.print_fingerprint(result);
  dca::testing::function_1d.print_fingerprint(result);
  dca::testing::function_2a.print_fingerprint(result);
  dca::testing::function_4a.print_fingerprint(result);
  dca::testing::function_16.print_fingerprint(result);

  EXPECT_TRUE(dca::testing::compare_to_file(DCA_SOURCE_DIR "/test/unit/function/fingerprint.txt",
                                            result.str()));
}

TEST(Function, PrintElements) {
  std::stringstream result;

  dca::testing::function_4a.print_elements(result);

  EXPECT_TRUE(dca::testing::compare_to_file(DCA_SOURCE_DIR "/test/unit/function/print_elements.txt",
                                            result.str()));
}

TEST(Function, to_JSON) {
  std::stringstream result;

  dca::util::print_type<dca::testing::test_domain_16::this_type>::to_JSON(std::cout);
  dca::util::print_type<dca::testing::test_domain_0a::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_0b::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_0c::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_0d::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_1d::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_2a::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_2c::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_4a::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::test_domain_16::this_type>::to_JSON(result);
  result << "\n";

  EXPECT_TRUE(
      dca::testing::compare_to_file(DCA_SOURCE_DIR "/test/unit/function/json.txt", result.str()));
}

TEST(FunctionTest, DefaultConstructor) {
  // Default name
  dca::func::function<double, dca::testing::test_domain_2a> f1;

  EXPECT_EQ(2, f1.signature());
  EXPECT_EQ(2, f1.size());

  for (int linind = 0; linind < f1.size(); ++linind)
    EXPECT_EQ(0., f1(linind));

  // Custom name
  const std::string name = "my-function";
  dca::func::function<double, dca::testing::test_domain_4a> f2(name);

  EXPECT_EQ(name, f2.get_name());
  EXPECT_EQ(4, f2.signature());
  EXPECT_EQ(64, f2.size());

  for (int linind = 0; linind < f2.size(); ++linind)
    EXPECT_EQ(0., f2(linind));
}

TEST(FunctionTest, CopyConstructor) {
  using FunctionType = dca::func::function<double, dca::testing::test_domain_2a>;

  FunctionType f1("original");
  f1(0) = 3.14;
  f1(1) = 2.72;

  // Default name
  FunctionType f2(f1);

  EXPECT_EQ(f1.signature(), f2.signature());
  EXPECT_EQ(f1.size(), f2.size());

  for (int linind = 0; linind < f2.size(); ++linind)
    EXPECT_EQ(f1(linind), f2(linind));

  // Same name
  FunctionType f3(f1);

  EXPECT_EQ(f1.get_name(), f3.get_name());
  EXPECT_EQ(f1.signature(), f3.signature());
  EXPECT_EQ(f1.size(), f3.size());

  for (int linind = 0; linind < f3.size(); ++linind)
    EXPECT_EQ(f1(linind), f3(linind));

  // Different name
  FunctionType f4(f3, "another name");
  EXPECT_EQ("another name", f4.get_name());
}

TEST(FunctionTest, MoveConstructor) {
  using FunctionType = dca::func::function<double, dca::testing::test_domain_2a>;

  FunctionType f1("original");
  f1(0) = 3.14;
  f1(1) = 2.72;

  // Default name
  FunctionType f1_copy(f1);
  FunctionType f2(std::move(f1_copy));

  EXPECT_EQ(f1.signature(), f2.signature());
  EXPECT_EQ(f1.size(), f2.size());

  for (int linind = 0; linind < f2.size(); ++linind)
    EXPECT_EQ(f1(linind), f2(linind));

  // Copy name
  FunctionType f1_copy_2(f1);
  FunctionType f3(std::move(f1_copy_2));

  EXPECT_EQ(f1.get_name(), f3.get_name());
  EXPECT_EQ(f1.signature(), f3.signature());
  EXPECT_EQ(f1.size(), f3.size());

  for (int linind = 0; linind < f3.size(); ++linind)
    EXPECT_EQ(f1(linind), f3(linind));

  // Different name
  FunctionType f4(std::move(f3), "another name");
  EXPECT_EQ("another name", f4.get_name());
}

TEST(FunctionTest, CopyAssignment) {
  using FunctionType = dca::func::function<double, dca::testing::test_domain_2a>;

  FunctionType f1("original");
  f1(0) = 3.14;
  f1(1) = 2.72;

  // Self-assignment
  f1 = f1;
  EXPECT_EQ("original", f1.get_name());
  EXPECT_EQ(2, f1.signature());
  EXPECT_EQ(2, f1.size());

  EXPECT_EQ(3.14, f1(0));
  EXPECT_EQ(2.72, f1(1));

  // Other assignment
  FunctionType f2("f2-assigned");
  f2 = f1;

  EXPECT_EQ("f2-assigned", f2.get_name());
  EXPECT_EQ(f1.signature(), f2.signature());
  EXPECT_EQ(f1.size(), f2.size());

  for (int linind = 0; linind < f2.size(); ++linind)
    EXPECT_EQ(f1(linind), f2(linind));

  // Check the return type of the copy assignment operator.
  FunctionType f3("f3-assigned");
  const FunctionType* const ptr_f3 = &f3;
  EXPECT_EQ(ptr_f3, &(f3 = f1));
}

TEST(FunctionTest, MoveAssignment) {
  using FunctionType = dca::func::function<double, dca::testing::test_domain_2a>;

  FunctionType f1("original");
  f1(0) = 3.14;
  f1(1) = 2.72;

  // Self-assignment
  f1 = std::move(f1);
  EXPECT_EQ("original", f1.get_name());
  EXPECT_EQ(2, f1.signature());
  EXPECT_EQ(2, f1.size());

  EXPECT_EQ(3.14, f1(0));
  EXPECT_EQ(2.72, f1(1));

  // Other assignment
  FunctionType f1_copy(f1);
  FunctionType f2("f2-assigned");
  f2 = std::move(f1_copy);

  EXPECT_EQ("f2-assigned", f2.get_name());
  EXPECT_EQ(f1.signature(), f2.signature());
  EXPECT_EQ(f1.size(), f2.size());

  for (int linind = 0; linind < f2.size(); ++linind)
    EXPECT_EQ(f1(linind), f2(linind));

  // Check the return type of the move assignment operator.
  FunctionType f1_copy_2(f1);
  FunctionType f3("f3-assigned");
  const FunctionType* const ptr_f3 = &f3;
  EXPECT_EQ(ptr_f3, &(f3 = std::move(f1)));
}

namespace dca {
namespace testing {
// dca::testing::

// Domain class with changeable size.
// Only for testing of the function class's reset method.
class VariableDmn {
public:
  using element_type = void;
  static void initialize(const int size) {
    size_ = size;
  }
  static int get_size() {
    return size_;
  }

private:
  static int size_;
};

int VariableDmn::size_ = 0;

}  // namespace testing
}  // namespace dca

TEST(FunctionTest, Reset) {
  // Size of test_domain_0b = 2.
  using Domain =
      dca::func::dmn_variadic<dca::testing::test_domain_0b, dca::func::dmn_0<dca::testing::VariableDmn>>;

  dca::testing::VariableDmn::initialize(3);

  dca::func::function<double, Domain> f;
  for (int i = 0; i < f.size(); ++i)
    f(i) = i;

  EXPECT_EQ(6, f.size());

  // Reinitialize with larger size.
  dca::testing::VariableDmn::initialize(10);

  EXPECT_EQ(6, f.size());

  f.reset();

  EXPECT_EQ(20, f.size());
  for (int i = 0; i < f.size(); ++i)
    EXPECT_EQ(0., f(i));

  // Set elements again to non-default (non-zero) values.
  for (int i = 0; i < f.size(); ++i)
    f(i) = i + 3.14;

  // Reinitialize with smaller size.
  dca::testing::VariableDmn::initialize(4);

  EXPECT_EQ(20, f.size());

  f.reset();

  EXPECT_EQ(8, f.size());
  for (int i = 0; i < f.size(); ++i)
    EXPECT_EQ(0., f(i));
}

TEST(FunctionTest, ComparisonOperatorEqual) {
  using FunctionType = dca::func::function<double, dca::testing::test_domain_2a>;

  FunctionType f1("f1");
  f1(0) = 3.14;
  f1(1) = 2.72;

  FunctionType f2("f2");
  f2(0) = 3.14;
  f2(1) = 2.72;

  FunctionType f3("f3");
  f3(0) = 3.14;
  f3(1) = 9.9;

  EXPECT_TRUE(f1 == f1);
  EXPECT_TRUE(f1 == f2);
  EXPECT_TRUE(f2 == f1);
  EXPECT_FALSE(f1 == f3);
}

TEST(FunctionTest, MemoryLayout) {
  // The function's data need to be stored in a column major order, as the code relies on this for
  // interactions with the Vector and Matrix classes.

  using Dmn1 = dmn_0<dmn<2>>;
  using Dmn2 = dmn_0<dmn<5>>;
  using Dmn3 = dmn_0<dmn<3>>;

  dca::func::function<int, dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>> f;

  std::vector<int> raw_data(f.size());
  int count = 0;

  for (auto& x : raw_data)
    x = ++count;

  // Copy from raw_data to f.
  std::copy_n(raw_data.data(), f.size(), f.values());

  // Check that the first index is the fast one.
  count = 0;
  for (int k = 0; k < Dmn3::dmn_size(); ++k)
    for (int j = 0; j < Dmn2::dmn_size(); ++j)
      for (int i = 0; i < Dmn1::dmn_size(); ++i) {
        EXPECT_EQ(raw_data[count], f(i, j, k));
        EXPECT_EQ(count, &f(i, j, k) - f.values());
        ++count;
      }
}
