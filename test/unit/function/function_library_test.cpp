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

#include "dca/testing/gtest_h_w_warning_blocking.h"
#include <array>
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

#include "dca/util/type_list.hpp"
#include "dca/util/type_help.hpp"
#include "function_testing.hpp"
#include <ModernStringUtils.hpp>

namespace dca {
namespace testing {
// dca::testing::

/** basically this ignores whitespace for comparing output of the type ouput.
 *  this is necessary because libc++ inserts spaces between trailing template type brackets which is
 * legal and was required by some older c++ stds.
 */
bool checkTypeLine(const std::string_view& expected, const std::string_view& test) {
  auto expected_tokens = modernstrutil::split(expected, " \t");
  auto test_tokens = modernstrutil::split(test, " \t");
  if (expected_tokens.size() != test_tokens.size())
    return false;
  for (int i = 0; i < expected_tokens.size(); ++i)
    if (expected_tokens[i] != test_tokens[i])
      return false;
  return true;
}

bool compare_to_file(const std::string& filename, const std::string& check) {
  // Open the file.
  std::ifstream known_result(filename);
  auto test_lines = modernstrutil::split(check, "\n");
  if (known_result.good()) {
    std::string expected_line;
    auto test_line = test_lines.begin();
    while (std::getline(known_result, expected_line)) {
      if (test_line == test_lines.end())
        return false;
      else if (!checkTypeLine(expected_line, *test_line))
        return false;
      ++test_line;
      return true;
    }
  }
  else {
    std::ofstream new_result(filename);
    auto test_line = test_lines.begin();
    new_result << check.c_str();
    std::cout << "No baseline file exists, writing new one " << filename.c_str() << std::endl;
  }
  return false;
}

dca::func::function<double, Domain0a> function_0a("Domain0a");
dca::func::function<double, Domain0b> function_0b("Domain0b");
dca::func::function<double, Domain0c> function_0c("Domain0c");
dca::func::function<double, Domain0d> function_0d("Domain0d");
dca::func::function<double, Domain1d> function_1d("Domain1d");
dca::func::function<double, Domain2a> function_2a("Domain2a");
dca::func::function<double, Domain4a> function_4a("Domain4a");
dca::func::function<double, Domain16> function_16("Domain16");

Domain0v dummy;
Domain16v dummy2;

}  // namespace testing
}  // namespace dca

TEST(Function, TestDomain4a) {
  try {
    // std::cout << "Leaf indexing \n";
    int index = 0;
    for (int i0 = 0; i0 < 1; ++i0) {
      for (int i1 = 0; i1 < 2; ++i1) {
        for (int i2 = 0; i2 < 4; ++i2) {
          for (int i3 = 0; i3 < 8; ++i3) {
            // std::cout << i0 << "," << i1 << "," << i2 << "," << i3 << "\n";
            dca::testing::function_4a.operator()(i0, i1, i2, i3) = index++;
            // dca::testing::function_4a.operator()(i3,i2,i1,i0) = index; bad ordering
          }
        }
      }
    }
    // std::cout << "Branch indexing \n";
    index = 0;
    for (int i0 = 0; i0 < 2; ++i0) {
      for (int i1 = 0; i1 < 32; ++i1) {
        // std::cout << i0 << "," << i1 << "\n";
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

  std::cout << result.str();
  EXPECT_TRUE(dca::testing::compare_to_file(DCA_SOURCE_DIR "/test/unit/function/fingerprint.txt",
                                            result.str()))
      << result.str();
}

TEST(Function, PrintElements) {
  std::stringstream result;

  dca::testing::function_4a.print_elements(result);

  EXPECT_TRUE(dca::testing::compare_to_file(DCA_SOURCE_DIR "/test/unit/function/print_elements.txt",
                                            result.str()));
}

TEST(Function, to_JSON) {
  std::stringstream result;

  dca::util::print_type<dca::testing::Domain16::this_type>::to_JSON(std::cout);
  dca::util::print_type<dca::testing::Domain0a::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain0b::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain0c::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain0d::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain1d::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain2a::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain2c::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain4a::this_type>::to_JSON(result);
  result << "\n";
  dca::util::print_type<dca::testing::Domain16::this_type>::to_JSON(result);
  result << "\n";

  EXPECT_TRUE(
      dca::testing::compare_to_file(DCA_SOURCE_DIR "/test/unit/function/json.txt", result.str()));
}

TEST(FunctionTest, DefaultConstructor) {
  // Default name
  dca::func::function<double, dca::testing::Domain2a> f1;

  EXPECT_EQ(2, f1.signature());
  EXPECT_EQ(2, f1.size());

  for (int linind = 0; linind < f1.size(); ++linind)
    EXPECT_EQ(0., f1(linind));

  // Custom name
  const std::string name = "my-function";
  dca::func::function<double, dca::testing::Domain4a> f2(name);

  EXPECT_EQ(name, f2.get_name());
  EXPECT_EQ(4, f2.signature());
  EXPECT_EQ(64, f2.size());

  for (int linind = 0; linind < f2.size(); ++linind)
    EXPECT_EQ(0., f2(linind));
}

TEST(FunctionTest, CopyConstructor) {
  using FunctionType = dca::func::function<double, dca::testing::Domain2a>;

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
  using FunctionType = dca::func::function<double, dca::testing::Domain2a>;

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
  using FunctionType = dca::func::function<double, dca::testing::Domain2a>;

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
  using FunctionType = dca::func::function<double, dca::testing::Domain2a>;

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
      dca::func::dmn_variadic<dca::testing::Domain0b, dca::func::dmn_0<dca::testing::VariableDmn>>;

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
  using FunctionType = dca::func::function<double, dca::testing::Domain2a>;

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

TEST(FunctionTest, RangeBasedLoop) {
  using Dmn = dmn_variadic<dmn_0<dmn<5>>, dmn_0<dmn<7>>>;

  dca::func::function<float, Dmn> f;

  int i = 0;
  for (auto& x : f)
    x = i++;

  for (int i = 0; i < f.size(); ++i)
    EXPECT_EQ(i, f(i));
}

TEST(FunctionTest, BranchIndexingAndAssignment) {
  using namespace dca::testing;
  using dca::func::dmn_variadic;
  using dca::vectorToString;

  dca::func::function<double, dmn_variadic<dmn_variadic<Domain2c, Domain2c>, Domain0c>> branched_function;
  dca::func::function<double, dmn_variadic<dmn_variadic<Domain2c, Domain2c>, Domain0c>> bf_test;

  function_test<decltype(branched_function)> branched_function_test(bf_test);
  branched_function_test.fill_sequence();
  int size_domain2c = Domain2c::dmn_size();
  EXPECT_EQ(size_domain2c, 32);
  auto& branched_steps = branched_function.get_domain().get_branch_domain_steps();
  // std::cout << "branched_steps: " << vectorToString(branched_steps) << '\n';
  auto& branch_sizes = branched_function.get_domain().get_branch_domain_sizes();
  // std::cout << "branched_sizes: " << vectorToString(branch_sizes) << '\n';
  auto get_steps = [](const std::vector<unsigned long>& sbdm_sizes) {
    std::vector<unsigned long> sbdm_steps(sbdm_sizes.size(), 1);
    for (int i = 1; i < sbdm_steps.size(); ++i)
      sbdm_steps[i] = sbdm_steps[i - 1] * sbdm_sizes[i - 1];
    return sbdm_steps;
  };

  auto other_steps = get_steps(branch_sizes);
  EXPECT_EQ(branched_steps, other_steps);
  ASSERT_EQ(bf_test(0, 0, 0, 0, 1), 1024);
  ASSERT_EQ(bf_test(0, 1, 0, 0, 0), 4);
  ASSERT_EQ(bf_test(0, 1), 1024);
  for (int c0_ind = 0; c0_ind < Domain0c::dmn_size(); ++c0_ind) {
    for (int j = 0; j < size_domain2c; ++j)
      for (int i = 0; i < size_domain2c; ++i) {
        branched_function(i + j * size_domain2c, c0_ind) =
            i + j * size_domain2c + c0_ind * size_domain2c * size_domain2c;
      }
  }

  EXPECT_EQ(branched_function.size(), 4096);
  EXPECT_EQ(branched_function.getValues().size(), 4096);
  // std::cout << vectorToString(branched_function.getValues()) << '\n';
  // std::cout << vectorToString(branched_function_test.f.getValues()) << '\n';
}

TEST(FunctionTest, ArrayBasedIndexing) {
  using namespace dca::testing;
  dca::func::function<double, Domain2c0c0c> f2c0c0c;
  function_test<decltype(f2c0c0c)> f2c0c0c_test(f2c0c0c);
  f2c0c0c_test.fill_sequence();
  std::array<int, 3> index{11, 0, 0};
  EXPECT_EQ(f2c0c0c(index), 11);
  index = {0, 1, 0};
  EXPECT_EQ(f2c0c0c(index), 32);

  // using it to to do remapping
  using DomainNot2c = dmn_variadic<Domain0d, Domain0c>;
  using DomainNot2c0c0c = dmn_variadic<DomainNot2c, Domain0c, Domain0c>;
  dca::func::function<double, DomainNot2c0c0c> fNot2c0c0c;
  std::array<int, 4> subind;
  std::array<int, 4> subind_transpose;
  for (int c1 = 0; c1 < Domain0c::dmn_size(); c1++)
    for (int c2 = 0; c2 < Domain0c::dmn_size(); c2++)
      for (int i2c = 0; i2c < Domain2c::dmn_size(); i2c++) {
        f2c0c0c.linind_2_subind(
            i2c + c2 * Domain2c::dmn_size() + c1 * Domain0c::dmn_size() * Domain2c::dmn_size(),
            subind);
        subind_transpose[1] = subind[0];
        subind_transpose[0] = subind[1];
        subind_transpose[2] = subind[2];
        subind_transpose[3] = subind[3];
        fNot2c0c0c(subind_transpose) = f2c0c0c(subind);
      }
  EXPECT_EQ(fNot2c0c0c(0, 0, 0, 0), f2c0c0c(0, 0, 0, 0));
  EXPECT_EQ(fNot2c0c0c(1, 0, 0, 0), f2c0c0c(0, 1, 0, 0));
  EXPECT_EQ(fNot2c0c0c(7, 3, 0, 0), f2c0c0c(3, 7, 0, 0));
  EXPECT_EQ(fNot2c0c0c(1, 0, 1, 0), f2c0c0c(0, 1, 1, 0));
  EXPECT_EQ(fNot2c0c0c(1, 0, 1, 0), f2c0c0c(0, 1, 1, 0));
  EXPECT_EQ(fNot2c0c0c(1, 0, 1, 1), f2c0c0c(0, 1, 1, 1));
  EXPECT_EQ(fNot2c0c0c(1, 0, 1, 1), f2c0c0c(0, 1, 1, 1));

  using Domain2c0a0c = dmn_variadic<Domain2c, Domain0a, Domain0c>;
  using DomainNot2c0a0c = dmn_variadic<DomainNot2c, Domain0a, Domain0c>;
  dca::func::function<double, Domain2c0a0c> f2c0a0c;
  function_test<decltype(f2c0a0c)> f2c0a0c_test(f2c0a0c);
  f2c0a0c_test.fill_sequence();
  dca::func::function<double, DomainNot2c0a0c> fNot2c0a0c;
  std::size_t linind = 0;
  for (int c1 = 0; c1 < Domain0c::dmn_size(); c1++)
    for (int c2 = 0; c2 < Domain0a::dmn_size(); c2++)
      for (int i2c = 0; i2c < Domain2c::dmn_size(); i2c++) {
        f2c0a0c.linind_2_subind(linind++, subind);
        subind_transpose[1] = subind[0];
        subind_transpose[0] = subind[1];
        subind_transpose[2] = subind[2];
        subind_transpose[3] = subind[3];
        fNot2c0a0c(subind_transpose) = f2c0a0c(subind);
      }
  EXPECT_EQ(fNot2c0a0c(0, 0, 0, 0), f2c0a0c(0, 0, 0, 0));
  EXPECT_EQ(fNot2c0a0c(1, 0, 0, 0), f2c0a0c(0, 1, 0, 0));
  EXPECT_EQ(fNot2c0a0c(7, 3, 0, 0), f2c0a0c(3, 7, 0, 0));
  EXPECT_EQ(fNot2c0a0c(1, 0, 0, 1), f2c0a0c(0, 1, 0, 1));
  EXPECT_EQ(fNot2c0a0c(7, 3, 0, 1), f2c0a0c(3, 7, 0, 1));
}
