// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests dmn_variadic.hpp.

#include "dca/function/domains/dmn_variadic.hpp"
#include "gtest/gtest.h"
#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

#ifndef NDEBUG

TEST(DmnVariadicTest, CheckIndices) {
  using dca::func::dmn_0;
  using dca::func::dmn;

  using SubDmn1 = dmn<2, int>;
  using SubDmn2 = dmn<4, double>;
  using SubDmn3 = dmn<3, unsigned int>;

  using ProductDmn = dca::func::dmn_variadic<dmn_0<SubDmn1>, dmn_0<SubDmn2>, dmn_0<SubDmn3>>;

  ProductDmn test_dmn;

  EXPECT_NO_THROW(test_dmn(0, 0, 0));
  EXPECT_NO_THROW(test_dmn(SubDmn1::get_size() - 1, 0, 0));
  EXPECT_NO_THROW(test_dmn(1, 3, 2));

  // Index too big.
  EXPECT_THROW(test_dmn(SubDmn1::get_size(), 0, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(0, 4, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(1, 0, 3), std::runtime_error);

  // Index too small.
  EXPECT_THROW(test_dmn(-1, 0, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(1, -1, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(0, 2, -3), std::runtime_error);
}

#endif  // NDEBUG

TEST(DmnVariadicTest, OperatorParentheses) {
  using dca::func::dmn_0;
  using dca::func::dmn;
  using SubDmn1 = dmn_0<dmn<2, short>>;
  using SubDmn2 = dmn_0<dmn<5, ulong>>;
  using SubDmn3 = dmn_0<dmn<7, float>>;
  using SubDmn4 = dmn_0<dmn<6, float>>;

  auto get_steps = [](const std::vector<int>& sbdm_sizes) {
    std::vector<int> sbdm_steps(sbdm_sizes.size(), 1);
    for (int i = 1; i < sbdm_steps.size(); ++i)
      sbdm_steps[i] = sbdm_steps[i - 1] * sbdm_sizes[i - 1];
    return sbdm_steps;
  };
  auto compute_linindex = [](const auto& indices, const auto& steps) {
    int linidex = 0;
    for (int i = 0; i < indices.size(); ++i)
      linidex += indices[i] * steps[i];
    return linidex;
  };

  using dca::func::dmn_variadic;
  {
    // Test a domain variadic given by a composition of simple subdomains.
    dmn_variadic<SubDmn1, SubDmn2, SubDmn3> domain;
    const std::vector<int> size{2, 5, 7};
    const auto steps = get_steps(size);

    std::array<int, 3> idx;
    for (idx[0] = 0; idx[0] < size[0]; ++idx[0])
      for (idx[1] = 0; idx[1] < size[1]; ++idx[1])
        for (idx[2] = 0; idx[2] < size[2]; ++idx[2])
          EXPECT_EQ(compute_linindex(idx, steps), domain(idx[0], idx[1], idx[2]));
  }

  {
    // Test a domain variadic given by a composition of domain variadics.
    dmn_variadic<dmn_variadic<dmn_variadic<SubDmn1, SubDmn2>, SubDmn3>, SubDmn4> domain;
    {
      // Access the domain with the innermost indices.
      const std::vector<int> size_leaves{2, 5, 7, 6};
      const auto steps = get_steps(size_leaves);

      std::array<int, 4> idx;
      for (idx[0] = 0; idx[0] < size_leaves[0]; ++idx[0])
        for (idx[1] = 0; idx[1] < size_leaves[1]; ++idx[1])
          for (idx[2] = 0; idx[2] < size_leaves[2]; ++idx[2])
            for (idx[3] = 0; idx[3] < size_leaves[3]; ++idx[3])
              EXPECT_EQ(compute_linindex(idx, steps), domain(idx[0], idx[1], idx[2], idx[3]));
    }
    {
      // Access the domain with the composite subdomain indices.
      const std::vector<int> size_branches{2 * 5 * 7, 6};
      const auto steps = get_steps(size_branches);

      std::array<int, 2> idx;
      for (idx[0] = 0; idx[0] < size_branches[0]; ++idx[0])
        for (idx[1] = 0; idx[1] < size_branches[1]; ++idx[1])
          EXPECT_EQ(compute_linindex(idx, steps), domain(idx[0], idx[1]));
    }
  }
}
