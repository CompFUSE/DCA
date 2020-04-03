// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class computes all single-particle functions.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_SP_GREENS_FUNCTION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_SP_GREENS_FUNCTION_HPP

#include <complex>
#include <iostream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/fermionic_overlap_matrices.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/fock_space.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/greens_functions/c_operator.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/greens_functions/sp_greens_function_data.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hamiltonian.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hilbert_spaces/hilbert_space.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameters_type, typename ed_options>
class SpGreensFunction {
public:
  using BDmn = typename ed_options::b_dmn;
  using SDmn = typename ed_options::s_dmn;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  using Profiler = typename ed_options::profiler_t;
  using Concurrency = typename ed_options::concurrency_type;

  using Scalar = typename ed_options::scalar_type;
  using Complex = typename ed_options::complex_type;

  using Vector = typename ed_options::vector_type;
  using Matrix = typename ed_options::matrix_type;
  using IntMatrix = typename ed_options::int_matrix_type;

  using NuDmn = typename ed_options::nu_dmn;
  using Bsr = typename ed_options::b_s_r;

  using BsDmn = typename ed_options::bs_dmn_type;
  using BsrDmn = typename ed_options::bsr_dmn_type;

  using NuNuRDmn = typename ed_options::nu_nu_r_dmn_type;

  using FermionicHamiltonian = Hamiltonian<parameters_type, ed_options>;
  using FermionicOverlap = fermionic_overlap_matrices<parameters_type, ed_options>;

  using FermionicFockSpace = Fock_space<parameters_type, ed_options>;
  using HilbertSpace = Hilbert_space<parameters_type, ed_options>;

  using FermionicFockDmn = func::dmn_0<FermionicFockSpace>;

  using SpGreensFunctionData = sp_Greens_function_data<ed_options>;

  using TDmn = func::dmn_0<domains::time_domain>;
  using WDmn = func::dmn_0<domains::frequency_domain>;
  using WRealDmn = func::dmn_0<domains::frequency_domain_real_axis>;
  using WVertexDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;

public:
  SpGreensFunction(parameters_type& parameters_ref, FermionicHamiltonian& Hamiltonian_ref,
                   FermionicOverlap& overlap_ref);

  template <typename MOMS_w_imag_type, typename MOMS_WRealDmn_type>
  void compute_all_sp_functions(MOMS_w_imag_type& MOMS_imag, MOMS_WRealDmn_type& MOMS_real,
                                bool interacting);

  template <typename MOMS_w_imag_type, typename MOMS_WRealDmn_type>
  void compute_all_sp_functions_slow(MOMS_w_imag_type& MOMS_imag, MOMS_WRealDmn_type& MOMS_real,
                                     bool interacting);

  template <typename w_dmn>
  void compute_S_k_w(
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, w_dmn>>& G_k_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, w_dmn>>& G0_k_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, w_dmn>>& S_k_w);

private:
  /*!
   *  new functions ...
   */
  void compute_real_space_Greens_functions(
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w_im,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>& G_r_w_re,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, RClusterDmn, WVertexDmn,
                                                 WVertexDmn>>& G_nu_nu_r_r_w_w,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, KClusterDmn, WVertexDmn,
                                                 WVertexDmn>>& G_nu_nu_k_k_w_w);

  void renormalize_real_space_Greens_functions(
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w_im,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>& G_r_w_re,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, RClusterDmn, WVertexDmn,
                                                 WVertexDmn>>& G_nu_nu_r_r_w_w,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, KClusterDmn, WVertexDmn,
                                                 WVertexDmn>>& G_nu_nu_k_k_w_w);

  void compute_Greens_functions_ac_slow(std::vector<SpGreensFunctionData>& data_vec);  // ,

  void compute_Greens_functions_ca_slow(std::vector<SpGreensFunctionData>& data_vec);  // ,

  void compute_sp_Greens_function(int nu_i_nu_j_delta_r, Scalar E_0, Scalar E_1, Complex factor,
                                  SpGreensFunctionData& data);

  /*!
   *   Here we try to compute sp-Greens-function as defined by its definition. It might be slower,
   *   but it is a good test to see that we do the right thing for the tp-Greens function
   */
  int has_nonzero_overlap(int HS_i, int HS_j, bool is_creation, int bsr_ind);

  void get_nonzero_overlap(int HS_i, int HS_j, bool is_creation, int bsr_ind, Matrix& matrix,
                           Matrix& tmp);

  void compute_sp_permutations(int bsr_0, int bsr_1, std::vector<std::vector<c_operator>>& sp_perms);

  Complex phi_e_ikr(int r0, int r1, int k0, int k1);

  Complex phi_e_iwt(Scalar E1, Scalar E2, int i1, int i2);

  void compute_sp_Greens_function_slow(int nu_i_nu_j_delta_r, Scalar E_i, Scalar E_j,
                                       Complex factor, SpGreensFunctionData& data);

  void compute_Greens_functions_slow(std::vector<SpGreensFunctionData>& data_vec);

  /*!
   *  old functions ...
   */
  void compute_Greens_functions_main(
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>&
          G_r_WRealDmn,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>>& G_k_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WRealDmn>>&
          G_k_WRealDmn,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>>& G_k_t);

  void compute_Greens_functions_st(
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>&
          G_r_WRealDmn,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t);

  void compute_Greens_functions_ac_translation(
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>&
          G_r_WRealDmn,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t);

  void compute_Greens_functions_ca_translation(
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>&
          G_r_WRealDmn,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t);

  void compute_sp_Greens_function(
      int nu_i, int nu_j, int delta_r, Scalar E_0, Scalar E_1, Complex factor,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w,
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>&
          G_r_WRealDmn,
      func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t);

private:
  parameters_type& parameters;
  Concurrency& concurrency;

  double CUT_OFF;

  FermionicHamiltonian& hamiltonian;
  FermionicOverlap& overlap;

  func::function<Vector, FermionicFockDmn>& eigen_energies;
  func::function<Matrix, FermionicFockDmn>& eigen_states;

  func::function<int, func::dmn_variadic<FermionicFockDmn, FermionicFockDmn, Bsr>>& creation_set_all;
  func::function<int, func::dmn_variadic<FermionicFockDmn, FermionicFockDmn, Bsr>>& annihilation_set_all;

  func::function<int, func::dmn_variadic<RClusterDmn, RClusterDmn>> rj_minus_ri;

  func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, RClusterDmn, WVertexDmn, WVertexDmn>>
      G0_nonlocal_nu_nu_r_r_w_w;
  func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, KClusterDmn, WVertexDmn, WVertexDmn>>
      G0_nonlocal_nu_nu_k_k_w_w;

  func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, RClusterDmn, WVertexDmn, WVertexDmn>>
      G_nonlocal_nu_nu_r_r_w_w;
  func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, KClusterDmn, WVertexDmn, WVertexDmn>>
      G_nonlocal_nu_nu_k_k_w_w;
};

template <typename parameters_type, typename ed_options>
SpGreensFunction<parameters_type, ed_options>::SpGreensFunction(parameters_type& parameters_ref,
                                                                FermionicHamiltonian& Hamiltonian_ref,
                                                                FermionicOverlap& overlap_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CUT_OFF(parameters.get_eigenvalue_cut_off()),

      hamiltonian(Hamiltonian_ref),
      overlap(overlap_ref),

      //       ac_Hilbert_space_indices(0),
      //       ca_Hilbert_space_indices(0),

      eigen_energies(hamiltonian.get_eigen_energies()),
      eigen_states(hamiltonian.get_eigen_states()),

      creation_set_all(overlap.get_creation_set_all()),
      annihilation_set_all(overlap.get_annihilation_set_all()),

      rj_minus_ri("rj_minus_ri"),

      G0_nonlocal_nu_nu_r_r_w_w("G0_nonlocal_nu_nu_r_r_w_w"),
      G0_nonlocal_nu_nu_k_k_w_w("G0_nonlocal_nu_nu_k_k_w_w"),

      G_nonlocal_nu_nu_r_r_w_w("G_nonlocal_nu_nu_r_r_w_w"),
      G_nonlocal_nu_nu_k_k_w_w("G_nonlocal_nu_nu_k_k_w_w") {
  for (int ri = 0; ri < RClusterDmn::dmn_size(); ri++)
    for (int rj = 0; rj < RClusterDmn::dmn_size(); rj++)
      rj_minus_ri(ri, rj) = RClusterDmn::parameter_type::subtract(ri, rj);
}

template <typename parameters_type, typename ed_options>
template <typename w_dmn>
void SpGreensFunction<parameters_type, ed_options>::compute_S_k_w(
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, w_dmn>>& G_k_w,
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, w_dmn>>& G0_k_w,
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, w_dmn>>& S_k_w) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t" << __FUNCTION__ << std::endl;

  int matrix_dim = BDmn::dmn_size() * SDmn::dmn_size();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inverted_matrix(matrix_dim);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_cluster_excluded_inverted_matrix(
      matrix_dim);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> sigma_matrix(matrix_dim);

  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w_dmn::dmn_size(); w_ind++) {
      dca::linalg::matrixop::copyArrayToMatrix(
          matrix_dim, matrix_dim, &G_k_w(0, 0, 0, 0, k_ind, w_ind), matrix_dim, G_inverted_matrix);
      dca::linalg::matrixop::inverse(G_inverted_matrix, ipiv, work);

      dca::linalg::matrixop::copyArrayToMatrix(matrix_dim, matrix_dim,
                                               &G0_k_w(0, 0, 0, 0, k_ind, w_ind), matrix_dim,
                                               G0_cluster_excluded_inverted_matrix);
      dca::linalg::matrixop::inverse(G0_cluster_excluded_inverted_matrix, ipiv, work);

      for (int j = 0; j < sigma_matrix.nrCols(); ++j)
        for (int i = 0; i < sigma_matrix.nrRows(); ++i)
          sigma_matrix(i, j) = G0_cluster_excluded_inverted_matrix(i, j) - G_inverted_matrix(i, j);

      dca::linalg::matrixop::copyMatrixToArray(sigma_matrix, &S_k_w(0, 0, 0, 0, k_ind, w_ind),
                                               matrix_dim);
    }
  }

  // Print out the low Matsubara frequency values of the self-energy.
  if (concurrency.id() == concurrency.first()) {
    int freqs = std::min(8, w_dmn::dmn_size() / 2);

    std::cout << "\n\t Sigma\n\n";
    for (int w_i = w_dmn::dmn_size() / 2 - freqs; w_i < w_dmn::dmn_size() / 2 + freqs; w_i++) {
      std::cout << w_dmn::get_elements()[w_i] << "\t";
      for (int k_i = 0; k_i < KClusterDmn::dmn_size(); k_i++)
        std::cout << real(S_k_w(0, 0, 0, 0, k_i, w_i)) << "\t" << imag(S_k_w(0, 0, 0, 0, k_i, w_i))
                  << "\t";
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

template <typename parameters_type, typename ed_options>
template <typename MOMS_w_imag_type, typename MOMS_WRealDmn_type>
void SpGreensFunction<parameters_type, ed_options>::compute_all_sp_functions_slow(
    MOMS_w_imag_type& MOMS_imag, MOMS_WRealDmn_type& MOMS_real, bool interacting) {
  if (interacting) {
    compute_real_space_Greens_functions(MOMS_imag.G_r_w, MOMS_real.G_r_w, MOMS_imag.G_r_t,
                                        G_nonlocal_nu_nu_r_r_w_w, G_nonlocal_nu_nu_k_k_w_w);
    renormalize_real_space_Greens_functions(MOMS_imag.G_r_w, MOMS_real.G_r_w, MOMS_imag.G_r_t,
                                            G_nonlocal_nu_nu_r_r_w_w, G_nonlocal_nu_nu_k_k_w_w);

    math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(MOMS_imag.G_r_w,
                                                                          MOMS_imag.G_k_w);
    math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(MOMS_real.G_r_w,
                                                                          MOMS_real.G_k_w);
    math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(MOMS_imag.G_r_t,
                                                                          MOMS_imag.G_k_t);
  }
  else {
    compute_real_space_Greens_functions(MOMS_imag.G0_r_w, MOMS_real.G0_r_w, MOMS_imag.G0_r_t,
                                        G0_nonlocal_nu_nu_r_r_w_w, G0_nonlocal_nu_nu_k_k_w_w);
    renormalize_real_space_Greens_functions(MOMS_imag.G0_r_w, MOMS_real.G0_r_w, MOMS_imag.G0_r_t,
                                            G0_nonlocal_nu_nu_r_r_w_w, G0_nonlocal_nu_nu_k_k_w_w);

    math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(MOMS_imag.G0_r_w,
                                                                          MOMS_imag.G0_k_w);
    math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(MOMS_real.G0_r_w,
                                                                          MOMS_real.G0_k_w);
    math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(MOMS_imag.G0_r_t,
                                                                          MOMS_imag.G0_k_t);
  }
}

template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::compute_real_space_Greens_functions(
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w_im,
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>& G_r_w_re,
    func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t,
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, RClusterDmn, WVertexDmn,
                                               WVertexDmn>>& /*G_nu_nu_r_r_w_w*/,
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, KClusterDmn, WVertexDmn,
                                               WVertexDmn>>& /*G_nu_nu_k_k_w_w*/) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t" << __FUNCTION__ << std::endl;

  G_r_w_im = 0;
  G_r_w_re = 0;
  G_r_t = 0;

  int start = clock();

  {
    int n_threads = 1;

    std::vector<SpGreensFunctionData> data_vec(n_threads);

    for (int l = 0; l < n_threads; l++)
      data_vec[l].initialize(parameters);

    compute_Greens_functions_slow(data_vec);

    for (int l = 0; l < n_threads; l++) {
      data_vec[l].sum_to(G_r_w_im, G_r_w_re, G_r_t);
    }
  }

  {
    for (int t_i = TDmn::dmn_size() / 2; t_i < TDmn::dmn_size(); t_i++)
      for (int r_i = 0; r_i < RClusterDmn::dmn_size(); r_i++)
        for (int nu_j = 0; nu_j < 2 * BDmn::dmn_size(); nu_j++)
          for (int nu_i = 0; nu_i < 2 * BDmn::dmn_size(); nu_i++)
            G_r_t(nu_i, nu_j, r_i, t_i - TDmn::dmn_size() / 2) = -G_r_t(nu_i, nu_j, r_i, t_i);
  }

  int end = clock();

  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__
              << " total time : " << double(end - start) / double(CLOCKS_PER_SEC) << std::endl;
}

template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::renormalize_real_space_Greens_functions(
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>& G_r_w_im,
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn>>& G_r_w_re,
    func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>>& G_r_t,
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, RClusterDmn, WVertexDmn,
                                               WVertexDmn>>& G_nu_nu_r_r_w_w,
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, KClusterDmn, WVertexDmn,
                                               WVertexDmn>>& G_nu_nu_k_k_w_w) {
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\t" << __FUNCTION__ << std::endl;
  }

  const std::vector<HilbertSpace>& Hilbert_spaces = FermionicFockDmn::get_elements();

  double beta = parameters.get_beta();

  double Z = 0;
  for (int HS_i = 0; HS_i < Hilbert_spaces.size(); ++HS_i)
    for (int n = 0; n < Hilbert_spaces[HS_i].size(); ++n)
      Z += std::exp(-beta * eigen_energies(HS_i)[n]);

  double factor = 1. / Z;

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\tfactor = " << factor << std::endl;
  }

  G_r_w_im *= factor;
  G_r_w_re *= factor;
  G_r_t *= -factor;

  G_nu_nu_r_r_w_w *= factor;
  G_nu_nu_k_k_w_w *= factor;
}

template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::compute_Greens_functions_ac_slow(
    std::vector<SpGreensFunctionData>& data_vec) {
  std::cout << "\n\n\t" << __FUNCTION__ << "\n\n";

  std::vector<HilbertSpace>& Hilbert_spaces = FermionicFockDmn::get_elements();

  Scalar beta = parameters.get_beta();

  double zero = 0, nonzero = 0;

  for (int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0) {
    for (int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1) {
      for (int ind = 0; ind < NuNuRDmn::dmn_size(); ind++) {
        SpGreensFunctionData& data = data_vec[0];

        data.set_indices(ind);

        bool done = false;

        if (annihilation_set_all(HS_0, HS_1, data.bsr_i) != -1 &&
            creation_set_all(HS_1, HS_0, data.bsr_j) != -1) {
          for (int n0 = 0; n0 < Hilbert_spaces[HS_0].size(); ++n0) {
            Scalar E_n0 = eigen_energies(HS_0)[n0];
            Scalar w_e = std::exp(-beta * E_n0);

            if (w_e > CUT_OFF) {
              for (int n1 = 0; n1 < Hilbert_spaces[HS_1].size(); ++n1) {
                if (not done) {
                  overlap.compute_annihilation_matrix_fast(HS_0, HS_1, data.bsr_i,
                                                           data.annihilation_bsr_i, data.tmp);
                  overlap.compute_creation_matrix_fast(HS_1, HS_0, data.bsr_j, data.creation_bsr_j,
                                                       data.tmp);

                  done = true;
                }

                Scalar E_n1 = eigen_energies(HS_1)[n1];

                Complex factor = w_e * data.annihilation_bsr_i(n0, n1) * data.creation_bsr_j(n1, n0);

                if (abs(factor) > CUT_OFF) {
                  nonzero += 1;

                  //                                       compute_sp_Greens_function(data.nu_i,
                  //                                       data.nu_j, data.delta_r,
                  //                                                                  E_n0, E_n1,
                  //                                                                  factor,
                  //                                                                  G_r_w,
                  //                                                                  G_r_WRealDmn,
                  //                                                                  G_r_t);

                  compute_sp_Greens_function(ind, E_n0, E_n1, factor, data);
                }
                else {
                  zero += 1;
                }
              }
            }
          }
        }
      }
    }
  }

  std::cout << "\n\t nonzero/total : " << nonzero / (zero + nonzero) << "\n";
}

template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::compute_Greens_functions_ca_slow(
    std::vector<SpGreensFunctionData>& data_vec)  // ,
//                                                                                                     func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn     > >& G_r_w,
//                                                                                                     func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WRealDmn> >& G_r_WRealDmn,
//                                                                                                     func::function<             double , func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, t     > >& G_r_t)
{
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t" << __FUNCTION__ << "\n\n";

  std::vector<HilbertSpace>& Hilbert_spaces = FermionicFockDmn::get_elements();

  Scalar beta = parameters.get_beta();

  double zero = 0, nonzero = 0;

  for (int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0) {
    for (int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1) {
      for (int ind = 0; ind < NuNuRDmn::dmn_size(); ind++) {
        SpGreensFunctionData& data = data_vec[0];

        data.set_indices(ind);

        bool done = false;

        if (annihilation_set_all(HS_1, HS_0, data.bsr_i) != -1 &&
            creation_set_all(HS_0, HS_1, data.bsr_j) != -1) {
          for (int n0 = 0; n0 < Hilbert_spaces[HS_0].size(); ++n0) {
            Scalar E_n0 = eigen_energies(HS_0)[n0];
            Scalar w_e = std::exp(-beta * E_n0);

            if (w_e > CUT_OFF) {
              for (int n1 = 0; n1 < Hilbert_spaces[HS_1].size(); ++n1) {
                if (not done) {
                  overlap.compute_annihilation_matrix_fast(HS_1, HS_0, data.bsr_i,
                                                           data.annihilation_bsr_i, data.tmp);
                  overlap.compute_creation_matrix_fast(HS_0, HS_1, data.bsr_j, data.creation_bsr_j,
                                                       data.tmp);

                  done = true;
                }

                Scalar E_n1 = eigen_energies(HS_1)[n1];

                Complex factor = w_e * data.creation_bsr_j(n0, n1) * data.annihilation_bsr_i(n1, n0);

                if (abs(factor) > CUT_OFF) {
                  nonzero += 1;

                  //                                       compute_sp_Greens_function(data.nu_i,
                  //                                       data.nu_j, data.delta_r,
                  //                                                                  E_n1, E_n0,
                  //                                                                  factor,
                  //                                                                  G_r_w,
                  //                                                                  G_r_WRealDmn,
                  //                                                                  G_r_t);

                  compute_sp_Greens_function(ind, E_n1, E_n0, factor, data);
                }
                else {
                  zero += 1;
                }
              }
            }
          }
        }
      }
    }
  }

  std::cout << "\n\t nonzero/total : " << nonzero / (zero + nonzero) << "\n";
}

/*!
 *   G(\tau) = e^(\tau*(E_0-E_1))/(e^(\beta*(E_0-E_1))+1)
 *   G(\omega) = 1/(i*\omega + (E_0-E_1))
 *   G(\tau) = 1/\beta \sum_{m=-\infty}^{\infty} G(\omega_m) e^{i*\omega_m*\tau}
 */
template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::compute_sp_Greens_function(
    int nu_i_nu_j_delta_r, Scalar E_0, Scalar E_1, Complex factor, SpGreensFunctionData& data) {
  Scalar beta = parameters.get_beta();

  Scalar ONE(1);

  {
    for (int t_i = TDmn::dmn_size() / 2; t_i < TDmn::dmn_size(); ++t_i)
      data.G_tau(t_i) = 1. / (std::exp((beta - data.tau(t_i)) * (E_0 - E_1)) +
                              std::exp(-data.tau(t_i) * (E_0 - E_1)));

    for (int w_i = 0; w_i < WDmn::dmn_size(); ++w_i)
      data.G_w_im(w_i) = ONE / (data.w_im(w_i) + E_0 - E_1);

    for (int w_i = 0; w_i < WRealDmn::dmn_size(); ++w_i)
      data.G_w_re(w_i) = ONE / (data.w_re(w_i) + E_0 - E_1);
  }

  {
    Scalar factor_re = real(factor);

    for (int t_i = TDmn::dmn_size() / 2; t_i < TDmn::dmn_size(); ++t_i)
      data.G_tau__nu_nu_r(t_i, nu_i_nu_j_delta_r) += factor_re * data.G_tau(t_i);

    for (int w_i = 0; w_i < WDmn::dmn_size(); ++w_i)
      data.G_w_im__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_im(w_i);

    for (int w_i = 0; w_i < WRealDmn::dmn_size(); ++w_i)
      data.G_w_re__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_re(w_i);
  }
}

template <typename parameters_type, typename ed_options>
int SpGreensFunction<parameters_type, ed_options>::has_nonzero_overlap(int HS_i, int HS_j,
                                                                       bool is_creation, int bsr_ind) {
  if (is_creation)
    return creation_set_all(HS_i, HS_j, bsr_ind);
  else
    return annihilation_set_all(HS_i, HS_j, bsr_ind);
}

template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::get_nonzero_overlap(int HS_i, int HS_j,
                                                                        bool is_creation,
                                                                        int bsr_ind, Matrix& matrix,
                                                                        Matrix& tmp) {
  if (is_creation)
    overlap.compute_creation_matrix_fast(HS_i, HS_j, bsr_ind, matrix, tmp);
  else
    overlap.compute_annihilation_matrix_fast(HS_i, HS_j, bsr_ind, matrix, tmp);
}

template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::compute_sp_permutations(
    int bsr_0, int bsr_1, std::vector<std::vector<c_operator>>& sp_perms) {
  sp_perms.resize(0);

  std::vector<c_operator> c_operators(2);

  {
    c_operators[0].index = 0;
    c_operators[1].index = 1;

    c_operators[0].bsr_ind = bsr_0;
    c_operators[1].bsr_ind = bsr_1;

    c_operators[0].creation = false;
    c_operators[1].creation = true;
  }

  sp_perms.push_back(c_operators);
}

/*!
 *   G(\tau)   = e^{(E_i-E_j)\tau}
 *   G(\omega) =
 *   G(\tau)   =
 */
template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::compute_sp_Greens_function_slow(
    int nu_i_nu_j_delta_r, Scalar E_i, Scalar E_j, Complex factor, SpGreensFunctionData& data) {
  Scalar beta = parameters.get_beta();

  Complex e_min_beta_E_i = std::exp(-beta * E_i);
  Complex e_min_beta_E_j = std::exp(-beta * E_j);

  if (abs(e_min_beta_E_i) > CUT_OFF or abs(e_min_beta_E_j) > CUT_OFF) {
    Scalar ONE(1);

    {
      data.G_tau = 0;
      for (int t_i = TDmn::dmn_size() / 2; t_i < TDmn::dmn_size(); ++t_i) {
        data.G_tau(t_i) =
            std::exp(-E_i * (beta - data.tau(t_i))) * std::exp(-E_j * (data.tau(t_i) - 0.));
      }

      data.G_w_im = 0;
      for (int w_i = 0; w_i < WDmn::dmn_size(); ++w_i) {
        Complex tmp = ONE / (data.w_im(w_i) + E_i - E_j);
        data.G_w_im(w_i) += std::exp(-beta * E_i) * tmp;
        data.G_w_im(w_i) += std::exp(-beta * E_j) * tmp;
      }

      data.G_w_re = 0;
      for (int w_i = 0; w_i < WRealDmn::dmn_size(); ++w_i) {
        Complex tmp = ONE / (data.w_re(w_i) + E_i - E_j);
        data.G_w_re(w_i) += std::exp(-beta * E_i) * tmp;
        data.G_w_re(w_i) += std::exp(-beta * E_j) * tmp;
      }
    }

    {
      Scalar factor_re = real(factor);

      for (int t_i = TDmn::dmn_size() / 2; t_i < TDmn::dmn_size(); ++t_i)
        data.G_tau__nu_nu_r(t_i, nu_i_nu_j_delta_r) += factor_re * data.G_tau(t_i);

      for (int w_i = 0; w_i < WDmn::dmn_size(); ++w_i)
        data.G_w_im__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_im(w_i);

      for (int w_i = 0; w_i < WRealDmn::dmn_size(); ++w_i)
        data.G_w_re__nu_nu_r(w_i, nu_i_nu_j_delta_r) += factor * data.G_w_re(w_i);
    }
  }
}

/*!
 *   G_{\nu, \mu}(R, \omega) = \frac{1}{Z} \sum_{i,j} e^{-\beta E_i} \phi(E_i, E_j, \omega) *
 * \langle i | c_{\nu, \vec{R}} | j \rangle * \langle j | c^{\dagger}_{\mu, \vec{0}} | j \rangle
 *   \phi(E_i, E_j, \omega)  = \int_0^{\beta} d\tau e^{((E_i-E_j)+i \omega)\tau}
 *                           = \frac{e^{((E_i-E_j)+i \omega)\beta}-1}{((E_i-E_j)+i \omega)}
 */
template <typename parameters_type, typename ed_options>
void SpGreensFunction<parameters_type, ed_options>::compute_Greens_functions_slow(
    std::vector<SpGreensFunctionData>& data_vec) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << std::endl;
  ;

  int origin = RClusterDmn::parameter_type::origin_index();

  const std::vector<HilbertSpace>& Hilbert_spaces = FermionicFockDmn::get_elements();

  for (int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0) {
    for (int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1) {
      for (int nu_0 = 0; nu_0 < NuDmn::dmn_size(); nu_0++) {
        for (int nu_1 = 0; nu_1 < NuDmn::dmn_size(); nu_1++) {
          for (int r_0 = 0; r_0 < RClusterDmn::dmn_size(); r_0++) {
            int r_1 = origin;

            SpGreensFunctionData& data = data_vec[0];

            int nu_0_nu_1_r_0 = data.nu_nu_r_dmn(nu_0, nu_1, r_0);

            int bsr_0 = data.nu_r_dmn(nu_0, r_0);
            int bsr_1 = data.nu_r_dmn(nu_1, r_1);

            std::vector<std::vector<c_operator>> sp_perms;

            compute_sp_permutations(bsr_0, bsr_1, sp_perms);

            for (int prm_ind = 0; prm_ind < sp_perms.size(); prm_ind++) {
              Scalar sign = 1;  //-((prm_ind%2)-0.5)*2.0;

              std::vector<c_operator>& operators = sp_perms[prm_ind];

              if (has_nonzero_overlap(HS_0, HS_1, operators[0].creation, operators[0].bsr_ind) != -1 &&
                  has_nonzero_overlap(HS_1, HS_0, operators[1].creation, operators[1].bsr_ind) != -1) {
                bool done = false;

                for (int l_0 = 0; l_0 < Hilbert_spaces[HS_0].size(); ++l_0) {
                  Scalar E_0 = eigen_energies(HS_0)[l_0];

                  for (int l_1 = 0; l_1 < Hilbert_spaces[HS_1].size(); ++l_1) {
                    Scalar E_1 = eigen_energies(HS_1)[l_1];

                    if (not done) {
                      get_nonzero_overlap(HS_0, HS_1, operators[0].creation, operators[0].bsr_ind,
                                          data.overlap_0, data.tmp);
                      get_nonzero_overlap(HS_1, HS_0, operators[1].creation, operators[1].bsr_ind,
                                          data.overlap_1, data.tmp);

                      done = true;
                    }

                    Complex factor = sign;  //*std::exp(-beta*E_0);

                    factor *= data.overlap_0(l_0, l_1);
                    factor *= data.overlap_1(l_1, l_0);

                    if (abs(factor) > CUT_OFF) {
                      compute_sp_Greens_function_slow(nu_0_nu_1_r_0, E_0, E_1, factor, data);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace ed
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_SP_GREENS_FUNCTION_HPP
