// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// Check the integration between the parallel distribution of large function objects
// provided by parallel util and various function methods
//

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/function/function.hpp"

#include <iostream>
#include <numeric>
#include <functional>

using dca::func::dmn;
using dca::func::dmn_0;
using dca::func::dmn_variadic;

namespace dca {
namespace testing {

template < class T >
std::ostream& operator<< (std::ostream& os, const std::vector<T>& v) 
{
    os << "[ ";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << "]";
    return os;
}

using Index = uint64_t;

// independent way to count number of elements in column major subindex span
int countSubIndexSpan(std::vector<int> sub_start, std::vector<int> sub_end, std::vector<size_t> dims)
{
  assert(sub_start.size() == sub_end.size());
  if (sub_start.size() > 1)
  {
    size_t mult = std::accumulate(dims.begin(), dims.end() - (dims.size() - sub_start.size() + 1), 1, std::multiplies<size_t>());
    return (sub_end.back() - sub_start.back()) *  mult + countSubIndexSpan(std::vector<int>(sub_start.begin(), sub_start.end() - 1),
                                                                         std::vector<int>(sub_end.begin(), sub_end.end() - 1),
                                                                             dims);
  }
  else if(sub_start.size() == 1)
    return sub_end.front() - sub_start.front() + 1;
  return 0;
}


TEST(FunctionTest, SubindicesOfDistributionOverLinearIndexBalanced) {
  // The function's data need to be stored in a column major order, as the code relies on this for
  // interactions with the Vector and Matrix classes.

  constexpr int dim1_size{15};
  constexpr int dim2_size{5};
  constexpr int dim3_size{3};
  using Dmn1 = dmn_0<dmn<dim1_size>>;
  using Dmn2 = dmn_0<dmn<dim2_size>>;
  using Dmn3 = dmn_0<dmn<dim3_size>>;

  dca::func::function<int, dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>> f;

  std::vector<Index> raw_data(f.size());
  Index count = 0;

  for (auto& x : raw_data)
    x = ++count;

  // Copy from raw_data to f.
  std::copy_n(raw_data.data(), f.size(), f.values());

  Index start = 0;
  Index end = 0;

  int num_ranks = 3;

  auto& this_domain = f.get_domain();
  auto& f_dims = this_domain.get_branch_domain_sizes();

  // balanced case
  for (int i = 0; i < num_ranks; ++i) {
    parallel::util::getComputeRange(i, num_ranks, f.size(), start, end);
    std::vector<int> subind_start = f.linind_2_subind(start);
    std::vector<int> subind_end = f.linind_2_subind(end);
    std::cerr << "rank: " << i << '\n'
              << "start lin indices " << start 
              << "start subindicies " << subind_start << '\n';
    std::cerr               << "end lin indices " << end 
              << "end subindicies " << subind_end << '\n';
    EXPECT_EQ(start, static_cast<uint64_t>(i) * ( raw_data.size() / num_ranks ));
    EXPECT_EQ(end - start + 1, raw_data.size() / num_ranks );
    EXPECT_EQ(countSubIndexSpan(subind_start, subind_end, f_dims), raw_data.size() / num_ranks);
  }
}

TEST(FunctionTest, SubindicesOfDistributionOverLinearIndexUnbalanced) {
  // The function's data need to be stored in a column major order, as the code relies on this for
  // interactions with the Vector and Matrix classes.

  constexpr int dim1_size{16};
  constexpr int dim2_size{5};
  constexpr int dim3_size{4};
  using Dmn1 = dmn_0<dmn<dim1_size>>;
  using Dmn2 = dmn_0<dmn<dim2_size>>;
  using Dmn3 = dmn_0<dmn<dim3_size>>;

  dca::func::function<int, dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>> f;

  std::vector<Index> raw_data(f.size());
  Index count = 0;

  for (auto& x : raw_data)
    x = ++count;

  // Copy from raw_data to f.
  std::copy_n(raw_data.data(), f.size(), f.values());

  Index start = 0;
  Index end = 0;

  int num_ranks = 3;

  int more_work_ranks = f.size() % static_cast<uint64_t>(num_ranks);

  auto& this_domain = f.get_domain();
  auto& f_dims = this_domain.get_branch_domain_sizes();
  // unbalanced case
  for (int i = 0; i < num_ranks; ++i) {
    parallel::util::getComputeRange(i, num_ranks, f.size(), start, end);
    std::vector<int> subind_start = f.linind_2_subind(start);
    std::vector<int> subind_end = f.linind_2_subind(end);
    std::cerr << "rank: " << i << '\n'
              << "start lin indices " << start 
              << "start subindicies " << subind_start << '\n';
    std::cerr               << "end lin indices " << end 
              << "end subindicies " << subind_end << '\n';

    if (i < more_work_ranks)
    {
      EXPECT_EQ(start, static_cast<uint64_t>(i) * ( raw_data.size() / num_ranks + 1));
      EXPECT_EQ(end - start + 1, raw_data.size() / num_ranks + 1 );
      EXPECT_EQ(countSubIndexSpan(subind_start, subind_end, f_dims), raw_data.size() / num_ranks + 1);
    }
    else
    {
      EXPECT_EQ(start, ( raw_data.size() / num_ranks + 1) * more_work_ranks + static_cast<uint64_t>(i - more_work_ranks) * ( raw_data.size() / num_ranks ));
      EXPECT_EQ(end - start + 1, raw_data.size() / num_ranks );
      EXPECT_EQ(countSubIndexSpan(subind_start, subind_end, f_dims), raw_data.size() / num_ranks);
    }
  }
}

}  // namespace testing
}  // namespace dca
