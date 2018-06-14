// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests vector_operations.hpp.

#include "dca/math/util/vector_operations.hpp"
#include <limits>
#include <vector>
#include "gtest/gtest.h"

template <typename ScalarType>
class VectorOperationsRealTest : public ::testing::Test {
protected:
  VectorOperationsRealTest() : x({1.23, 4.56, 7.89}), y({3.14, 2.72, 1.0}) {}

  static const ScalarType epsilon;

  const std::vector<ScalarType> x;
  const std::vector<ScalarType> y;
};
template <typename ScalarType>
const ScalarType VectorOperationsRealTest<ScalarType>::epsilon =
    std::numeric_limits<ScalarType>::epsilon();

using FloatingPointTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(VectorOperationsRealTest, FloatingPointTypes);

template <typename ComplexType>
class VectorOperationsComplexTest : public ::testing::Test {
protected:
  VectorOperationsComplexTest()
      : v({ComplexType(1., 2.), ComplexType(3., 3.14)}),
        w({ComplexType(3.1, -1.), ComplexType(2.7, 1.)}) {}

  static const typename ComplexType::value_type epsilon;

  const std::vector<ComplexType> v;
  const std::vector<ComplexType> w;
};
template <typename ComplexType>
const typename ComplexType::value_type VectorOperationsComplexTest<ComplexType>::epsilon =
    std::numeric_limits<typename ComplexType::value_type>::epsilon();

using ComplexFloatingPointTypes = ::testing::Types<std::complex<float>, std::complex<double>>;
TYPED_TEST_CASE(VectorOperationsComplexTest, ComplexFloatingPointTypes);

//
// Print
//
TYPED_TEST(VectorOperationsRealTest, Print) {
  dca::math::util::print(this->x);
}

TYPED_TEST(VectorOperationsComplexTest, Print) {
  dca::math::util::print(this->v);
}

//
// Scale
//
TYPED_TEST(VectorOperationsRealTest, Scale) {
  using ScalarType = TypeParam;
  const ScalarType a = 3.14;

  const std::vector<ScalarType> res = dca::math::util::scale(a, this->x);

  EXPECT_NEAR(ScalarType(3.8622), res[0], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(14.3184), res[1], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(24.7746), res[2], 500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, Scale) {
  using ComplexType = TypeParam;
  const ComplexType a(1.5, 2.0);

  const std::vector<ComplexType> res = dca::math::util::scale(a, this->v);

  EXPECT_NEAR(typename ComplexType::value_type(-2.5), res[0].real(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(5.), res[0].imag(), 500 * this->epsilon);
}

//
// Add
//
TYPED_TEST(VectorOperationsRealTest, Add) {
  using ScalarType = TypeParam;
  const std::vector<ScalarType> res = dca::math::util::add(this->x, this->y);

  EXPECT_NEAR(ScalarType(4.37), res[0], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(7.28), res[1], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(8.89), res[2], 500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, Add) {
  using ComplexType = TypeParam;
  const std::vector<ComplexType> res = dca::math::util::add(this->v, this->w);

  EXPECT_NEAR(typename ComplexType::value_type(4.1), res[0].real(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(1.), res[0].imag(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(5.7), res[1].real(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(4.14), res[1].imag(), 500 * this->epsilon);
}

//
// Subtract
//
TYPED_TEST(VectorOperationsRealTest, Subtract) {
  using ScalarType = TypeParam;
  const std::vector<ScalarType> res = dca::math::util::subtract(this->x, this->y);

  EXPECT_NEAR(ScalarType(1.91), res[0], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(-1.84), res[1], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(-6.89), res[2], 500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, Subtract) {
  using ComplexType = TypeParam;
  const std::vector<ComplexType> res = dca::math::util::subtract(this->v, this->w);

  EXPECT_NEAR(typename ComplexType::value_type(2.1), res[0].real(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(-3.), res[0].imag(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(-0.3), res[1].real(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(-2.14), res[1].imag(), 500 * this->epsilon);
}

//
// Inner product
//
TYPED_TEST(VectorOperationsRealTest, InnerProduct) {
  using ScalarType = TypeParam;
  EXPECT_NEAR(ScalarType(24.1554), dca::math::util::innerProduct(this->x, this->y),
              500 * this->epsilon);
}

TYPED_TEST(VectorOperationsRealTest, ScalarInnerProduct) {
  using ScalarType = TypeParam;
  ScalarType x(3.14), y(42);
  EXPECT_EQ(x * y, dca::math::util::innerProduct(x, y));
}

TYPED_TEST(VectorOperationsComplexTest, InnerProduct) {
  using ComplexType = TypeParam;
  const ComplexType res = dca::math::util::innerProduct(this->v, this->w);

  EXPECT_NEAR(typename ComplexType::value_type(12.34), res.real(), 500 * this->epsilon);
  EXPECT_NEAR(typename ComplexType::value_type(12.678), res.imag(), 500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, ScalarInnerProduct) {
  using ComplexType = TypeParam;
  ComplexType x(3.14, -1), y(42, 2.71);
  EXPECT_EQ(x * std::conj(y), dca::math::util::innerProduct(x, y));
}

//
// L^2 norm squared
//
TYPED_TEST(VectorOperationsRealTest, L2Norm2) {
  using ScalarType = TypeParam;
  EXPECT_NEAR(ScalarType(84.5586), dca::math::util::l2Norm2(this->x), 500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, L2Norm2) {
  using ComplexType = TypeParam;
  EXPECT_NEAR(typename ComplexType::value_type(23.8596), dca::math::util::l2Norm2(this->v),
              500 * this->epsilon);
}

//
// L^2 norm squared
//
TYPED_TEST(VectorOperationsRealTest, L2Norm) {
  using ScalarType = TypeParam;
  EXPECT_NEAR(ScalarType(9.19557502280309), dca::math::util::l2Norm(this->x), 500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, L2Norm) {
  using ComplexType = TypeParam;
  EXPECT_NEAR(typename ComplexType::value_type(4.88462895213137), dca::math::util::l2Norm(this->v),
              500 * this->epsilon);
}

//
// Distance squared
//
TYPED_TEST(VectorOperationsRealTest, Distance2) {
  using ScalarType = TypeParam;
  EXPECT_NEAR(ScalarType(54.5058), dca::math::util::distance2(this->x, this->y), 500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, Distance2) {
  using ComplexType = TypeParam;
  EXPECT_NEAR(typename ComplexType::value_type(18.0796),
              dca::math::util::distance2(this->v, this->w), 500 * this->epsilon);
}

//
// Distance
//
TYPED_TEST(VectorOperationsRealTest, Distance) {
  using ScalarType = TypeParam;
  EXPECT_NEAR(ScalarType(7.38280434523359), dca::math::util::distance(this->x, this->y),
              500 * this->epsilon);
}

TYPED_TEST(VectorOperationsComplexTest, Distance) {
  using ComplexType = TypeParam;
  EXPECT_NEAR(typename ComplexType::value_type(4.25201128879028),
              dca::math::util::distance(this->v, this->w), 500 * this->epsilon);
}

//
// Area
//
TEST(VectorOperationsTest, Area) {
  const std::vector<int> a{2, 1};
  const std::vector<int> b{-4, 7};

  EXPECT_EQ(18, dca::math::util::area(a, b));
}

//
// Volume
//
TEST(VectorOperationsTest, Volume) {
  const std::vector<int> a{3, 2, 1};
  const std::vector<int> b{-1, 3, 0};
  const std::vector<int> c{2, 2, 5};

  EXPECT_EQ(47, dca::math::util::volume(a, b, c));
}

//
// Coordinates
//
TYPED_TEST(VectorOperationsRealTest, Coordinates) {
  using ScalarType = TypeParam;
  std::vector<ScalarType> r{3, 2, 1};
  std::vector<std::vector<ScalarType>> basis{std::vector<ScalarType>{1, 0, 0},
                                             std::vector<ScalarType>{1, 1, 0},
                                             std::vector<ScalarType>{1, 1, 1}};
  std::vector<ScalarType> res = dca::math::util::coordinates(r, basis);

  EXPECT_NEAR(ScalarType(1), res[0], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(1), res[1], 500 * this->epsilon);
  EXPECT_NEAR(ScalarType(1), res[2], 500 * this->epsilon);
}

//
// Is less vector
//
TYPED_TEST(VectorOperationsRealTest, IsLessVector) {
  using ScalarType = TypeParam;
  const std::vector<ScalarType> v1{0., 1.};
  const std::vector<ScalarType> v2{0., 2.};
  const std::vector<ScalarType> v3{1., 0.};

  EXPECT_TRUE(dca::math::util::isLessVector(v1, v2));
  EXPECT_FALSE(dca::math::util::isLessVector(v2, v1));
  EXPECT_TRUE(dca::math::util::isLessVector(v1, v3));
  EXPECT_FALSE(dca::math::util::isLessVector(v1, v1));
}

//
// Has smaller norm
//
TYPED_TEST(VectorOperationsRealTest, HasSmallerNorm) {
  using ScalarType = TypeParam;

  // Squares of the norms differ by less than the default tolerance = 1.e-6.
  const std::vector<ScalarType> v1{0., 1.};
  const std::vector<ScalarType> v2{1., 0.};
  EXPECT_TRUE(dca::math::util::hasSmallerNorm(v1, v2));
  EXPECT_FALSE(dca::math::util::hasSmallerNorm(v1, v1));

  // Squares of the norms differ by more than the default tolerance = 1.e-6.
  const std::vector<ScalarType> v3{0., 2.};
  EXPECT_TRUE(dca::math::util::hasSmallerNorm(v1, v3));
}

//
// Is same vector
//
TYPED_TEST(VectorOperationsRealTest, IsSameVector) {
  using ScalarType = TypeParam;
  const std::vector<ScalarType> v1{0., 1.};
  const std::vector<ScalarType> v2{0., 2.};
  const std::vector<ScalarType> v3{1.e-6, 1.};

  EXPECT_TRUE(dca::math::util::isSameVector(v1, v1));
  EXPECT_FALSE(dca::math::util::isSameVector(v1, v2));
  EXPECT_FALSE(dca::math::util::isSameVector(v2, v1));
  EXPECT_TRUE(dca::math::util::isSameVector(v1, v3));
  EXPECT_TRUE(dca::math::util::isSameVector(v3, v1));
}
