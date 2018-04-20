// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Jérémie Bouquet (bouquetj@gmail.com).
//
// This class tests the CPU walker used by the ctint cluster solver. The fast updated matrix
// are compared with their direct computation.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "gtest/gtest.h"

#include "walker_wrapper_submatrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

using G0Setup = typename dca::testing::G0Setup<dca::testing::LatticeSquare>;
using namespace dca::phys::solver;
using Walker = testing::phys::solver::ctint::WalkerWrapperSubmatrix<G0Setup::Parameters>;
using Matrix = Walker::Matrix;
using MatrixPair = std::array<Matrix, 2>;

TEST_F(G0Setup, RemoveAndInstertVertex) {
  std::vector<double> rng_values(1000);
  for (double& x : rng_values)
    x = double(std::rand()) / RAND_MAX;
  G0Setup::RngType rng(rng_values);

  ctint::G0Interpolation<dca::linalg::CPU> g0(
      dca::phys::solver::ctint::details::shrinkG0(data->G0_r_t));
  G0Setup::LabelDomain label_dmn;
  ctint::DMatrixBuilder<dca::linalg::CPU> builder(g0, Rdmn::parameter_type::get_subtract_matrix(),
                                                  label_dmn.get_branch_domain_steps(),
                                                  parameters.getAlphas());
  Walker walker(parameters, rng, G0Setup::interaction_vertices, builder);

  // *******************************
  // Test vertex insert/removal ****
  // *******************************
  // Set rng values.
  
  rng.setNewValues(std::vector<double>{0.95, 0.01, 0.2, 0.8, 0.01, 0.98, 0.9, 0.3, 0.32, 0.81, 0.71,
	0.5, 0.1, 0.88, 0.04, 0.1, 0.7, 0.1, 0.9, 0.22, 0.7, 0.76, 0.9, 0.43, 0.4, 0.1, 0.9, 0.2,
	0.4, 0.12, 0.97, 0.3, 0.6, 0.54, 0.87, 0.3, 0.08, 0.12, 0.54, 0.98});
  
  MatrixPair old_M(walker.getM());
  
  walker.doStep(10);
  
  MatrixPair new_M(walker.getM());

  // Compute directly the new M.
  
  walker.setMFromConfig();
  MatrixPair direct_M(walker.getM());
  
  for (int s = 0; s < 2; ++s)
    for (int j = 0; j < new_M[s].nrCols(); ++j)
      for (int i = 0; i < new_M[s].nrRows(); ++i) {
        EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 1e-7);
      }
  
  // *******************************
  // Test removal of all vertices***
  // *******************************

  old_M = walker.getM();

  int size = old_M[0].size().first;
  
  std::vector<double> values;
  for (int i = 0; i < size * 4; ++i) {
    values.push_back(0.99);
  }
  
  rng.setNewValues(values);

  walker.doStep(size - 1);
  new_M = walker.getM();
  
  //walker.setMFromConfig();
  direct_M = walker.getM();

  for (int s = 0; s < 2; ++s)
    EXPECT_NEAR(direct_M[s](0,0), new_M[s](0,0), 1e-7);
}   

