// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests cluster_domain.hpp

#include <vector>
#include <array>

#include "dca/platform/dca_gpu.h"
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/domains/cluster/cluster_domain_symmetry_initializer.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/parameters/analysis_parameters.hpp"
#include "dca/phys/parameters/domains_parameters.hpp"
#include "dca/phys/parameters/model_parameters.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/util/to_string.hpp"
#include "dca/io/json/json_reader.hpp"

class ClusterDomainsTest : public ::testing::Test {
public:
  using Lattice = dca::phys::models::square_lattice<dca::phys::domains::no_symmetry<2>>;
  using Model = dca::phys::models::TightBindingModel<Lattice>;
  using CDA = dca::phys::ClusterDomainAliases<Lattice::DIMENSION>;

protected:
  static void setupTestCase() {}
};


TEST(ClusterDomainsTest, initializeKDmn) {
  dca::io::JSONReader reader;
  using DomainsParameters = dca::phys::params::DomainsParameters;
  using AnalysisParameters = dca::phys::params::AnalysisParameters;
  using Model = ClusterDomainsTest::Model;

  class ClusterDomainsTestParameters : public AnalysisParameters, public DomainsParameters
  {
  public:
    ClusterDomainsTestParameters() : AnalysisParameters(Model::DIMENSION), DomainsParameters(Model::DIMENSION) {}
    void read(dca::io::JSONReader& reader) {
      AnalysisParameters::readWrite(reader);
      DomainsParameters::readWrite(reader);
    }      
  };

  ClusterDomainsTestParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/domains/cluster/input.json");
  pars.read(reader);

  reader.close_file();
  
  using dca::vectorToString;

  using CDA = dca::phys::ClusterDomainAliases<ClusterDomainsTest::Lattice::DIMENSION>;

  std::array<double, 4> basis{{2., 0., 0., 2.}};  // Lattice basis: [1, 0], [0, 1].
  std::vector<std::vector<int>> cluster{{20, 20}, {20, -20}};
  dca::phys::domains::cluster_domain_initializer<CDA::RClusterDmn>::execute(basis.data(), cluster);

  dca::phys::domains::cluster_domain_symmetry_initializer<
      CDA::RClusterDmn, typename Model::lattice_type::DCA_point_group>::execute();

  dca::phys::domains::cluster_domain_initializer<CDA::RSpHostDmn>::execute(Model::get_r_DCA_basis(),
                                                                           pars.get_sp_host());
  dca::phys::domains::cluster_domain_symmetry_initializer<
      CDA::RSpHostDmn, typename Model::lattice_type::DCA_point_group>::execute();
  // side effect of above is KSpHostDmn initialized

  EXPECT_EQ(CDA::KSpHostDmn::dmn_size(), 800);   

  // for (int k_ind = 0; k_ind < CDA::KSpHostDmn::dmn_size(); k_ind++)
  //   std::cout << vectorToString(CDA::KSpHostDmn::get_elements()[k_ind]) << '\n';

  dca::phys::domains::cluster_domain_initializer<CDA::RQHostDmn>::execute(Model::get_r_DCA_basis(),
                                                           pars.get_q_host());
  dca::phys::domains::cluster_domain_symmetry_initializer<
    CDA::RQHostDmn, typename Model::lattice_type::DCA_point_group>::execute();

  EXPECT_EQ(CDA::KQHostDmn::dmn_size(), 18);

  dca::phys::domains::cluster_domain_initializer<CDA::RQFineDmn>::execute(Model::get_r_DCA_basis(),
                                                           pars.get_q_host_fine());
  dca::phys::domains::cluster_domain_symmetry_initializer<
    CDA::RQFineDmn, typename Model::lattice_type::DCA_point_group>::execute();

  EXPECT_EQ(CDA::KQFineDmn::dmn_size(), 100);

}
