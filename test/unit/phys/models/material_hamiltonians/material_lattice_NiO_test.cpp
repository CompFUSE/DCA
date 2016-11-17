// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the specialization of material_lattice for NiO.

#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

// Helper structs since we can only template tests on types.
namespace dca {
namespace testing {
// dca::testing::

struct NiOSymmetricStruct {
  static constexpr phys::models::material_name_type type = phys::models::NiO_symmetric;
};

struct NiOUnsymmetricStruct {
  static constexpr phys::models::material_name_type type = phys::models::NiO_unsymmetric;
};

}  // testing
}  // dca

template <typename T>
class MaterialLatticeNiOTest : public ::testing::Test {};

typedef ::testing::Types<dca::testing::NiOSymmetricStruct, dca::testing::NiOUnsymmetricStruct> NiOTypes;
TYPED_TEST_CASE(MaterialLatticeNiOTest, NiOTypes);

TYPED_TEST(MaterialLatticeNiOTest, Initialize_H_0) {
  using Lattice =
      dca::phys::models::material_lattice<TypeParam::type, dca::phys::domains::no_symmetry<3>>;

  using BandDmn = dca::func::dmn<8, int>;
  using SpinDmn = dca::func::dmn<2, int>;
  using BandSpinDmn = dca::func::dmn_variadic<dca::func::dmn_0<BandDmn>, dca::func::dmn_0<SpinDmn>>;

  using KDmn = dca::func::dmn<3, std::vector<double>>;
  const double a = Lattice::latticeConstant();
  KDmn::set_elements({{0., 0., 0.}, {0., M_PI / a, 0.}, {M_PI / a, M_PI / a, M_PI / a}});

  dca::func::function<std::complex<double>,
                      dca::func::dmn_variadic<BandSpinDmn, BandSpinDmn, dca::func::dmn_0<KDmn>>>
      H_0;

  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<Lattice>> params;
  params.set_t_ij_file_name(DCA_SOURCE_DIR
                            "/include/dca/phys/models/material_hamiltonians/NiO/t_ij_NiO.txt");

  Lattice::initialize_H_0(params, H_0);

  // All imaginary parts should be smaller than 10^-3.
  for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
    for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1)
      for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
        for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2)
          for (int k = 0; k < KDmn::dmn_size(); ++k)
            EXPECT_LE(std::abs(H_0(b1, s1, b2, s2, k).imag()), 1.e-3);

  // All matrix elements with different spin indices should be zero.
  for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
    for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
      for (int k = 0; k < KDmn::dmn_size(); ++k) {
        EXPECT_DOUBLE_EQ(0., H_0(b1, 0, b2, 1, k).real());
        EXPECT_DOUBLE_EQ(0., H_0(b1, 1, b2, 0, k).real());
      }

  // Check some nonvanishing Hamiltonian matrix elements.
  EXPECT_DOUBLE_EQ(-1.7838558854405486, H_0(0, 0, 0, 0, 0).real());
  EXPECT_DOUBLE_EQ(-4.452512149016615, H_0(5, 1, 5, 1, 2).real());
  EXPECT_DOUBLE_EQ(1.428376402198317, H_0(6, 1, 5, 1, 2).real());
}
