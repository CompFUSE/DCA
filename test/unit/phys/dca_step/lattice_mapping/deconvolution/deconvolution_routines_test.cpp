// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests deconvolution_routines.hpp.

#include "dca/phys/dca_step/lattice_mapping/deconvolution/deconvolution_routines.hpp"

#include "gtest/gtest.h"

#include "dca/config/threading.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

using namespace dca;

TEST(DeconvolutionRoutinesTest, ProjectionOperator) {
  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::square_lattice<PointGroup>;
  using Model = phys::models::TightBindingModel<Lattice>;

  using ConcurrencyType = parallel::NoConcurrency;
  using ParametersType =
      phys::params::Parameters<ConcurrencyType, Threading, profiling::NullProfiler, Model,
                               void /*RandomNumberGenerator*/, phys::solver::CT_AUX>;
  using KSourceDmn = func::dmn_0<
      phys::domains::cluster_domain<double, Lattice::DIMENSION, phys::domains::CLUSTER,
                                    phys::domains::MOMENTUM_SPACE, phys::domains::BRILLOUIN_ZONE>>;
  using KTargetDmn = func::dmn_0<
      phys::domains::cluster_domain<double, Lattice::DIMENSION, phys::domains::LATTICE_SP,
                                    phys::domains::MOMENTUM_SPACE, phys::domains::BRILLOUIN_ZONE>>;

  ConcurrencyType concurrency(0, nullptr);

  ParametersType parameters("", concurrency);
  parameters.read_input_and_broadcast<io::JSONReader>(DCA_SOURCE_DIR
                                                      "/test/unit/phys/dca_step/lattice_mapping/"
                                                      "deconvolution/"
                                                      "deconvolution_routines_test_input.json");
  parameters.update_model();
  parameters.update_domains();

  phys::latticemapping::deconvolution_routines<ParametersType, KSourceDmn, KTargetDmn> deconv_routines(
      parameters);

  // Target (lattice) k-domain to target k-domain projection operators
  const auto& projection_op = deconv_routines.get_T();
  const auto& projection_op_symmetrized = deconv_routines.get_T_symmetrized();

  // projection_op.print();
  // projection_op_symmetrized.print();

  // Check the first rows of the projection operator matrices.
  EXPECT_DOUBLE_EQ(0.66963106982612841, projection_op(0, 0));
  EXPECT_DOUBLE_EQ(0.14867881635766231, projection_op(0, 1));
  EXPECT_DOUBLE_EQ(0.14867881635766231, projection_op(0, 2));
  EXPECT_DOUBLE_EQ(0.033011297458547091, projection_op(0, 3));

  EXPECT_DOUBLE_EQ(0.66963106982612841, projection_op_symmetrized(0, 0));
  EXPECT_DOUBLE_EQ(0.14867881635766231, projection_op_symmetrized(0, 1));
  EXPECT_DOUBLE_EQ(0.14867881635766231, projection_op_symmetrized(0, 2));
  EXPECT_DOUBLE_EQ(0.033011297458547091, projection_op_symmetrized(0, 3));

  // Target (lattice) k-domain to source (cluster) k-domain projection operators
  const auto& projection_op_source = deconv_routines.get_T_source();
  const auto& projection_op_source_symmetrized = deconv_routines.get_T_source_symmetrized();

  // projection_op_source.print();
  // projection_op_source_symmetrized.print();

  // Source (cluster) k-domain and target (lattice) k-domain are identical in this example (see
  // input file).
  for (int j = 0; j < projection_op.size().second; j++)
    for (int i = 0; i < projection_op.size().second; i++) {
      EXPECT_DOUBLE_EQ(projection_op(i, j), projection_op_source(i, j));
      EXPECT_DOUBLE_EQ(projection_op_symmetrized(i, j), projection_op_source_symmetrized(i, j));
    }
}
