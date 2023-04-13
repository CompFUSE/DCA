// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// This file tests model traits

#include "dca/phys/parameters/model_parameters.hpp"

#include "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"

#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "dca/phys/models/traits.hpp"

namespace dca {
namespace phys {
namespace params {
template <typename Lattice>
class MockParameters
    : public phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> {
public:
  constexpr static int lattice_dimension = Lattice::DIMENSION;
  using Scalar = double;
  using lattice_type = Lattice;

  using CDA = ClusterDomainAliases<lattice_dimension>;
  // DCA cluster domains
  using RClusterDmn = typename CDA::RClusterDmn;
};
}  // namespace params
}  // namespace phys
}  // namespace dca

using namespace dca;

TEST(ModelTraitsTest, HasInitializeNonDensityInteractionMethod) {
  using Lattice = phys::models::HundLattice<phys::domains::D4>;
  dca::phys::params::MockParameters<Lattice> hund_parameters;
  phys::models::HundLattice<phys::domains::D4>::printNonDensityType(hund_parameters);
  bool has_ndim = false;
  if constexpr (phys::models::HasInitializeNonDensityInteractionMethod<decltype(hund_parameters)>::value)
    has_ndim = true;
  EXPECT_TRUE(has_ndim);
  using Lattice2 = phys::models::bilayer_lattice<phys::domains::D4>;
  dca::phys::params::MockParameters<Lattice2> bilayer_parameters;
  has_ndim = false;
  if constexpr (phys::models::HasInitializeNonDensityInteractionMethod<decltype(bilayer_parameters)>::value)
    has_ndim = true;
  EXPECT_FALSE(has_ndim);  
}
