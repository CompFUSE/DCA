// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the insertion of spin flip and pair hopping terms

#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hamiltonian.hpp"

#include "gtest/gtest.h"
#include <string>

#include "dca/function/function_utils.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/options.hpp"
#include "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp"
#include "test/integration/exact_diagonalization_advanced/ed_cluster_solver_test_helper.hpp"

using Parameters = dca::testing::Parameters;
using Options = dca::phys::solver::ed::Options<Parameters>;
using Matrix = dca::linalg::Matrix<double, dca::linalg::CPU>;

Matrix referenceHamiltonian(const double U, const double mu, const double jh);

TEST(HamiltonianTest, BuildHamiltonian) {
  namespace test = dca::testing;
  using Model = dca::phys::models::HundLattice<dca::phys::domains::D4>;
  using ModelWrp = dca::phys::models::TightBindingModel<Model>;
  using Parameters = dca::phys::params::Parameters<dca::parallel::NoConcurrency, void, void,
                                                   ModelWrp, void, dca::phys::solver::CT_INT>;

  dca::parallel::NoConcurrency concurrency(0, nullptr);
  Parameters pars("", concurrency);
  const std::string input_dir = DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/exact_diagonalization_advanced/";
  pars.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "nd_test_input.json");
  pars.update_domains();

  using dca::func::function;
  using dca::func::dmn_variadic;
  using Nu = test::OrbitalSpinDmn;
  using test::RDmn;
  using test::KDmn;
  function<double, dmn_variadic<Nu, Nu, RDmn>> H_int("H_int");
  function<std::complex<double>, dmn_variadic<Nu, Nu, KDmn>> H0_k("H0_k");
  function<std::complex<double>, dmn_variadic<Nu, Nu, RDmn>> H0_r("H0_r");
  function<double, dmn_variadic<Nu, Nu, Nu, Nu, RDmn>> H_nd("H_nd");

  // Initialize Hamiltonians
  Model::initialize_H_interaction(H_int, pars);
  Model::initializeNonDensityInteraction(H_nd, pars);
  Model::initialize_H_0(pars, H0_k);
  dca::math::transform::FunctionTransform<KDmn, RDmn>::execute(H0_k, H0_r);

  dca::phys::solver::ed::Fock_space<Parameters, Options> fock_obj(true, true);
  //  Don't apply symmetry for test: fock_obj.apply_translation_symmetry();
  fock_obj.initialize_rep();
  using Hamiltonian = dca::phys::solver::ed::Hamiltonian<Parameters, Options>;
  Hamiltonian hamiltonian_obj(pars);

  hamiltonian_obj.initialize(dca::func::utils::real(H0_r), H_int, H_nd);
  hamiltonian_obj.construct_Hamiltonians(true);

  // Build reference Hamiltonian
  const int test_sector = 4;  // Index of up-down sector of the Fock space.
  Matrix ham_check = referenceHamiltonian(pars.get_U(), pars.get_chemical_potential(), pars.get_Jh());

  std::cout << "Solver matrix:\n";
  const auto& solver_matrix = hamiltonian_obj.getHamiltonian(test_sector);
  solver_matrix.print();

  for (int j = 0; j < ham_check.nrCols(); j++)
    for (int i = 0; i < ham_check.nrRows(); i++)
      EXPECT_NEAR(ham_check(i, j), solver_matrix(i, j).real(), 1e-12);
}

// Build the expected matrix for the up-down sector.
Matrix referenceHamiltonian(const double U, const double mu, const double jh) {
  Matrix H(4, 4); // up-down Hamiltonian.
  const int np = 2, nb = 2;  // number of particles and bands.
  const double shift = -np * mu - U / 2. * np + nb * U / 4.;
  // The basis is: {|ud,0>, |d,u>, |u,d>, |0,ud>}
  // Note: with p.b.c the electrons can hop in two directions on the same site, hence the 2*t.
  H(0, 0) = U + shift, H(0, 1) = 0., H(0, 2) = 0., H(0, 3) = jh;
  H(1, 0) = 0, H(1, 1) = shift, H(1, 2) = jh, H(1, 3) = 0.;
  H(2, 0) = 0, H(2, 1) = jh, H(2, 2) = shift, H(2, 3) = 0.;
  H(3, 0) = jh, H(3, 1) = 0., H(3, 2) = 0., H(3, 3) = U + shift;

  return H;
}
