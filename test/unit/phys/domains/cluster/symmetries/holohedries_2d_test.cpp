// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Tyler Sax
//
// Completeness oracle for the 2D holohedry candidate pool. Proves the pool is a complete
// crystallographic op set: filtering the pool by each 2D Bravais lattice reproduces that
// lattice's full holohedry (square -> D4 order 8, hexagonal -> D6 order 12, rectangular -> D2
// order 4, oblique -> C2 order 2). Because every one of the ten 2D crystallographic point groups
// is a subgroup of one of these holohedries, realizing all four maximal groups proves the pool
// misses no op.

#include <array>
#include <cmath>
#include <vector>

#include "dca/phys/domains/cluster/symmetries/point_groups/2d/holohedries_2d.hpp"

#include "gtest/gtest.h"

#include "dca/util/type_list.hpp"

namespace {

// A 2x2 op stored row-major: m[i][j].
using Mat2 = std::array<std::array<double, 2>, 2>;

constexpr double kTol = 1e-10;

// The op classes store their matrix column-major (m[i + 2*j] = row i, col j); convert to row-major.
Mat2 fromColumnMajor(const double* m) {
  return {{{m[0], m[2]}, {m[1], m[3]}}};
}

// Materialize a point_group_type_list into its op matrices. Recurses over the typelist by index so
// the same helper serves any pool/group.
template <class List, int I, int N>
struct Collect {
  static void run(std::vector<Mat2>& out) {
    using Op = typename dca::util::TypeAt<I, List>::type;
    out.push_back(fromColumnMajor(Op::matrix()));
    Collect<List, I + 1, N>::run(out);
  }
};
template <class List, int N>
struct Collect<List, N, N> {
  static void run(std::vector<Mat2>&) {}
};

template <class PointGroup>
std::vector<Mat2> collectOps() {
  using List = typename PointGroup::point_group_type_list;
  std::vector<Mat2> out;
  Collect<List, 0, dca::util::Length<List>::value>::run(out);
  return out;
}

// Note that the algebra below is implemented directly in this test, as opposed
// to being imported via the existing linalg path. Given this is a logic-only
// test operating on 2x2 matrices, hand rolling makes more sense than adding
// a dependency on LAPACK/BLAS, which is overkill here.
bool equal(const Mat2& a, const Mat2& b) {
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      if (std::abs(a[i][j] - b[i][j]) > kTol)
        return false;
  return true;
}

std::vector<Mat2> dedup(const std::vector<Mat2>& ops) {
  std::vector<Mat2> out;
  for (const auto& o : ops) {
    bool seen = false;
    for (const auto& k : out)
      seen = seen || equal(o, k);
    if (!seen)
      out.push_back(o);
  }
  return out;
}

Mat2 matmul(const Mat2& a, const Mat2& b) {
  Mat2 c{};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j];
  return c;
}

double det(const Mat2& a) {
  return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

Mat2 inverse(const Mat2& a) {
  const double d = det(a);
  return {{{a[1][1] / d, -a[0][1] / d}, {-a[1][0] / d, a[0][0] / d}}};
}

bool isInteger(const Mat2& a) {
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      if (std::abs(a[i][j] - std::round(a[i][j])) > 1e-8)
        return false;
  return true;
}

// Order of the holohedry of the lattice whose primitive vectors are the columns of
// `basis`: the number of pool ops that map the lattice onto itself (verified by applying
// the operation to the basis vectors, multiplying by the inverse of the basis, and
// checking whether the resulting matrix has integer values.)
int holohedryOrder(const std::vector<Mat2>& ops, const Mat2& basis) {
  const Mat2 binv = inverse(basis);
  int count = 0;
  for (const auto& O : ops)
    if (isInteger(matmul(binv, matmul(O, basis))))
      ++count;
  return count;
}

const double kHalfRoot3 = std::sqrt(3.0) / 2.0;

}  // namespace

// Completeness: filtering the pool by each 2D Bravais lattice yields that lattice's full holohedry.
// Covering all four maximal groups (D4, D6, D2, C2) covers all ten 2D point groups as subgroups, so
// this is the load-bearing check that the pool misses no op.
//
// The pool lists D4 (8) + D6 (12) = 20 ops; the shared {identity, C2, mirror@0, mirror@90} collapse
// under dedup to 16 distinct ops, which is what the per-lattice filter runs over.
TEST(HolohedryPool2DTest, RealizesEvery2DHolohedry) {
  const auto ops = dedup(collectOps<dca::phys::domains::holohedry_pool_2D>());
  ASSERT_EQ(ops.size(), 16u);  // 20 listed, 16 distinct after the D4/D6 overlap collapses

  const Mat2 square{{{1, 0}, {0, 1}}};
  const Mat2 hexagonal{{{1, -0.5}, {0, kHalfRoot3}}};
  const Mat2 rectangular{{{1, 0}, {0, 2}}};
  const Mat2 oblique{{{1, 0}, {0.3, 1.7}}};

  EXPECT_EQ(holohedryOrder(ops, square), 8);       // D4
  EXPECT_EQ(holohedryOrder(ops, hexagonal), 12);   // D6
  EXPECT_EQ(holohedryOrder(ops, rectangular), 4);  // D2
  EXPECT_EQ(holohedryOrder(ops, oblique), 2);      // C2
}
