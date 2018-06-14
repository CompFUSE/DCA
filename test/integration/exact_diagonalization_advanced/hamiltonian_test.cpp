// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the construction of the ED Hamiltonian on a two-site square-lattice Hubbard model
// in the one-particle up-down sector.

#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hamiltonian.hpp"

#include <complex>

#include "gtest/gtest.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/fock_space.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/options.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

using namespace dca;

using Matrix = linalg::Matrix<double, linalg::CPU>;

// Returns the Hamiltonian matrix for the one-particle up-down sector.
Matrix referenceHamiltonian(const double t, const double U, const double mu) {
  Matrix H("up-down-Hamiltonian", 4);

  // The basis is {|ud, 0>, |d, u>, |u, d>, |0, ud>}.
  // Note that the factor of 2 in the non-vanishing, off-diagonal matrix elements originates from
  // the periodic boundary conditions. In a two-site system (along x), electrons hopping in the +/-
  // x-direction end up on the same site in this case.
  H(0, 0) = U - mu, H(0, 1) = 2 * t, H(0, 2) = -2 * t, H(0, 3) = 0.;
  H(1, 0) = 2 * t, H(1, 1) = -mu, H(1, 2) = 0., H(1, 3) = 2 * t;
  H(2, 0) = -2 * t, H(2, 1) = 0., H(2, 2) = -mu, H(2, 3) = -2 * t;
  H(3, 0) = 0., H(3, 1) = 2 * t, H(3, 2) = -2 * t, H(3, 3) = U - mu;

  return H;
}

TEST(HamiltonianTest, ConstructHamiltonian) {
  using Lattice = phys::models::square_lattice<phys::domains::D4>;
  using Model = phys::models::TightBindingModel<Lattice>;
  using Parameters = phys::params::Parameters<parallel::NoConcurrency, void, void, Model, void,
                                              phys::solver::CT_AUX>;  // CT_AUX is a placeholder.
  using EdOptions = phys::solver::ed::Options<Parameters>;

  using OrbitalSpinDmn = func::dmn_variadic<func::dmn_0<phys::domains::electron_band_domain>,
                                            func::dmn_0<phys::domains::electron_spin_domain>>;
  using RDmn = func::dmn_0<phys::domains::cluster_domain<
      double, 2, phys::domains::CLUSTER, phys::domains::REAL_SPACE, phys::domains::BRILLOUIN_ZONE>>;

  parallel::NoConcurrency concurrency(0, nullptr);

  Parameters params("", concurrency);
  params.read_input_and_broadcast<io::JSONReader>(
      DCA_SOURCE_DIR
      "/test/integration/exact_diagonalization_advanced/hamiltonian_test_input.json");
  params.update_domains();

  const double t = 1.;  // Hopping parameter.
  const double U = 5.;  // On-site interaction.

  // Non-interacting Hamiltonian
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, RDmn>> H_0(
      "H_0(r)");
  H_0(0, 0, 1) = H_0(1, 1, 1) = -2.;  // Index 1 corresponds to the lattice vector that connects
                                      // nearest neighors.

  // Interaction Hamiltonian
  func::function<double, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, RDmn>> H_int;
  H_int(0, 1, 0) = H_int(1, 0, 0) = U;  // Index 0 corresponds to the origin.

  phys::solver::ed::Fock_space<Parameters, EdOptions> fock_obj(
      true, true);  // true, true: construct with occupation and magnetization symmetry.
  fock_obj.initialize_rep();

  phys::solver::ed::Hamiltonian<Parameters, EdOptions> hamiltonian_obj(params);
  hamiltonian_obj.initialize(H_0, H_int);
  hamiltonian_obj.construct_Hamiltonians(true);  // true: include the interaction part.

  Matrix ham_matrix_ref = referenceHamiltonian(t, U, params.get_chemical_potential());

  const int test_sector = 4;  // Index of the one-particle up-down sector in the Fock space.
  const auto& ham_matrix = hamiltonian_obj.get_Hamiltonian(test_sector);
  ham_matrix.print();

  for (int j = 0; j < ham_matrix_ref.nrCols(); ++j) {
    for (int i = 0; i < ham_matrix_ref.nrRows(); ++i) {
      EXPECT_DOUBLE_EQ(ham_matrix_ref(i, j), ham_matrix(i, j).real());
      EXPECT_DOUBLE_EQ(0., ham_matrix(i, j).imag());
    }
  }
}
