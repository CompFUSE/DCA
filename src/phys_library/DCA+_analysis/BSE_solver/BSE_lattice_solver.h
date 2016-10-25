// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Mi Jiang (jiangmi@itp.phys.ethz.ch)
//
// Description

/*
  Modifications compared with BSE_lattice_solver_original.h
  ==========================================================
  BSE_lattice_solver_Peter1:
  1. N_LAMBDAS = 4096; 2. only symmetrize 5 leading eigenvalues
  2a. print out Gamma_lattice

  ==========================================================
  BSE_lattice_solver_Peter2:
  3. comment out characterize_leading_eigenvectors()
  4. symmetrize Gamma_lattice (peoject away odd-frequency part)

  ==========================================================
  BSE_lattice_solver_Peter3:
  5. use sqrt(chi0)*gamma*sqrt(chi0) instead of gamma*chi0

  ==========================================================
  BSE_lattice_solver_Peter4:
  6. use real sqrt(chi0)*gamma*sqrt(chi0) symmetric matrix for particle-particle SC channel
  7. delete some useless routines
  8. always diagonalize full sqrt(chi0)*gamma*sqrt(chi0)
  9. For real symmetric matrix, do not symmetrize_leading_eigenvectors()
  10. print out Gamma_sym
  11. print out sigma_lattice(k,w)

  ==========================================================
  BSE_lattice_solver_Peter5:
  12. only store N_LAMBDA leading eigenvalues and eigenvectors to decrease size of BSE.hdf5.
  Maybe need all eigenpairs for susceptibility
  this way do not need change N_LAMBDA for different cluster sizes for regular DCA, keep using
  N_LAMBDAS = 10

  ==========================================================
  BSE_lattice_solver_Peter6:
  13. only do lattice mapping for Gamma, chi0 and diagonalization done via python code
  by comment out // BSE_lattice_solver_obj.diagonolize_Gamma_chi_0(Gamma_lattice, chi_0) in
  BSE_solver.h
  14. move symmetrize Gamma_lattice for even-frequency part for Gamma_sym into
  compute_Gamma_lattice_3,
  instead of diagonolize_full_Gamma_chi_0
  15. only write Gamma and sigma_lattice data into hdf5
*/

#ifndef PHYS_LIBRARY_DCA_ANALYSIS_BSE_SOLVER_BSE_LATTICE_SOLVER_H
#define PHYS_LIBRARY_DCA_ANALYSIS_BSE_SOLVER_BSE_LATTICE_SOLVER_H

#include <algorithm>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/util/comparison_methods.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/util/print_time.hpp"
#include "comp_library/linalg/linalg.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/high_temperature_series_expansion_solver.h"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_tp.h"
#include "phys_library/DCA+_step/lattice_mapping/lattice_mapping_tp.h"
#include "phys_library/DCA+_step/lattice_mapping/lattice_mapping_sp.h"
#include "phys_library/DCA+_step/symmetrization/diagrammatic_symmetries.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize.h"
#include "phys_library/domains/cluster/centered_cluster_domain.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"

using namespace dca;

namespace DCA {

template <class parameters_type, class MOMS_type>
class BSE_lattice_solver {
public:
  using scalartype = double;

  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  const static int N_LAMBDAS = 10;
  using lambda_dmn_type = func::dmn_0<func::dmn<N_LAMBDAS, int>>;

  const static int N_CUBIC = 3;
  using cubic_harmonics_dmn_type = func::dmn_0<func::dmn<N_CUBIC, int>>;

  using w_VERTEX = func::dmn_0<DCA::vertex_frequency_domain<DCA::COMPACT>>;
  using b = func::dmn_0<electron_band_domain>;
  using b_b = func::dmn_variadic<b, b>;

  using k_DCA = func::dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                           CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST = func::dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                            LATTICE_SP, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST_VERTEX = func::dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                                   LATTICE_TP, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  using host_vertex_r_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                                    LATTICE_TP, REAL_SPACE, BRILLOUIN_ZONE>;
  using crystal_harmonics_expansion = centered_cluster_domain<host_vertex_r_cluster_type>;
  using crystal_harmonics_expansion_dmn_t = func::dmn_0<crystal_harmonics_expansion>;

  using chi_vector_dmn_t = func::dmn_variadic<b, b, crystal_harmonics_expansion_dmn_t>;

  using cluster_eigenvector_dmn_t = func::dmn_variadic<b, b, k_DCA, w_VERTEX>;
  using lattice_eigenvector_dmn_t = func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>;
  using crystal_eigenvector_dmn_t =
      func::dmn_variadic<b, b, crystal_harmonics_expansion_dmn_t, w_VERTEX>;
  using cubic_eigenvector_dmn_t = func::dmn_variadic<b, b, cubic_harmonics_dmn_type, w_VERTEX>;

  using DCA_matrix_dmn_t = func::dmn_variadic<cluster_eigenvector_dmn_t, cluster_eigenvector_dmn_t>;
  using HOST_matrix_dmn_t = func::dmn_variadic<lattice_eigenvector_dmn_t, lattice_eigenvector_dmn_t>;
  using crystal_matrix_dmn_t =
      func::dmn_variadic<crystal_eigenvector_dmn_t, crystal_eigenvector_dmn_t>;

public:
  BSE_lattice_solver(parameters_type& parameters, MOMS_type& MOMS);

  func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& get_Gamma_lattice() {
    return Gamma_lattice;
  }

  template <typename Writer>
  void write(Writer& writer);

  void compute_chi_0_lattice(func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0);

  void compute_Gamma_lattice_1(func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster);
  void compute_Gamma_lattice_3(func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster);

  void diagonolize_Gamma_chi_0(
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

  func::function<std::complex<scalartype>, lambda_dmn_type>& get_leading_eigenvalues() {
    return leading_eigenvalues;
  };

private:
  void initialize();

  void set_chi_0_matrix(func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0);

  void diagonolize_full_Gamma_chi_0_real(
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

  void diagonolize_full_Gamma_chi_0(
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

  void record_eigenvalues_and_eigenvectors(dca::linalg::Vector<scalartype, dca::linalg::CPU>& L,
                                           dca::linalg::Matrix<scalartype, dca::linalg::CPU>& VR);

  void record_eigenvalues_and_eigenvectors(
      dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU>& L,
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VL,
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VR);

  void diagonolize_folded_Gamma_chi_0(
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice);

  void record_eigenvalues_and_folded_eigenvectors(
      dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU>& L,
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VL,
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VR);

  void compute_folded_susceptibility(
      func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice,
      dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU>& L,
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VL,
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VR);

  void symmetrize_leading_eigenvectors();

  void characterize_leading_eigenvectors();

  void print_on_shell();
  void print_on_shell_ppSC();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  MOMS_type& MOMS;

  func::function<std::complex<scalartype>, HOST_matrix_dmn_t> Gamma_lattice;
  func::function<std::complex<scalartype>, HOST_matrix_dmn_t> Gamma_sym;
  func::function<std::complex<scalartype>, func::dmn_variadic<b_b, b_b, k_HOST_VERTEX, w_VERTEX>>
      chi_0_function;

  func::function<std::complex<scalartype>, lambda_dmn_type> leading_eigenvalues;
  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, cubic_eigenvector_dmn_t>>
      leading_symmetry_decomposition;
  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, lattice_eigenvector_dmn_t>>
      leading_eigenvectors;

  func::function<scalartype, lambda_dmn_type> leading_eigenvalues_real;
  func::function<scalartype, func::dmn_variadic<lambda_dmn_type, lattice_eigenvector_dmn_t>>
      leading_eigenvectors_real;

  func::function<std::complex<scalartype>, func::dmn_variadic<chi_vector_dmn_t, chi_vector_dmn_t>> chi_q;

  func::function<std::complex<scalartype>,
                 func::dmn_variadic<k_HOST_VERTEX, crystal_harmonics_expansion_dmn_t>>
      psi_k;
  func::function<std::complex<scalartype>,
                 func::dmn_variadic<lattice_eigenvector_dmn_t, crystal_eigenvector_dmn_t>>
      crystal_harmonics;

  func::function<std::complex<scalartype>, func::dmn_variadic<k_HOST_VERTEX, cubic_harmonics_dmn_type>>
      leading_symmetry_functions;
};

template <class parameters_type, class MOMS_type>
BSE_lattice_solver<parameters_type, MOMS_type>::BSE_lattice_solver(parameters_type& parameters_ref,
                                                                   MOMS_type& MOMS_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      MOMS(MOMS_ref),

      Gamma_lattice("Gamma_lattice"),
      Gamma_sym("Gamma_sym"),
      chi_0_function("chi-0-function"),

      leading_eigenvalues("leading-eigenvalues"),
      leading_symmetry_decomposition("leading-symmetry-decomposition"),
      leading_eigenvectors("leading-eigenvectors"),

      leading_eigenvalues_real("leading-eigenvalues-real"),
      leading_eigenvectors_real("leading-eigenvectors-real"),

      chi_q("chi(q)"),

      psi_k("psi_k"),
      crystal_harmonics("crystal_harmonics"),

      leading_symmetry_functions("leading-symmetry-functions")
//     leading_U_K("leading_U_K"),
//     leading_Vt_K("leading_Vt_K"),

//     leading_U_k("leading_U_k_interpolated"),
//     leading_Vt_k("leading_Vt_k_interpolated")
{
  initialize();
}

template <class parameters_type, class MOMS_type>
template <typename Writer>
void BSE_lattice_solver<parameters_type, MOMS_type>::write(Writer& writer) {
  if (true) {
    writer.execute(leading_eigenvalues);
    writer.execute(leading_eigenvectors);

    writer.execute(leading_symmetry_decomposition);
    writer.execute(leading_symmetry_functions);

    writer.execute(Gamma_lattice);
    writer.execute(chi_0_function);
  }

  else {
    writer.execute(leading_eigenvalues_real);
    writer.execute(leading_eigenvectors_real);

    writer.execute(Gamma_sym);
    writer.execute(chi_0_function);
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::initialize() {
  {
    double r_cut_off = parameters.get_BSE_cut_off_radius();

    crystal_harmonics_expansion::initialize();

    std::vector<std::vector<double>> r_vecs(0);

    for (int l = 0; l < crystal_harmonics_expansion::get_size(); l++)
      if (math::util::l2Norm(crystal_harmonics_expansion::get_elements()[l]) < r_cut_off)
        r_vecs.push_back(crystal_harmonics_expansion::get_elements()[l]);

    sort(r_vecs.begin(), r_vecs.end(), math::util::hasSmallerNorm<double>);

    crystal_harmonics_expansion::get_size() = r_vecs.size();
    crystal_harmonics_expansion::get_elements() = r_vecs;

    if (concurrency.id() == concurrency.last()) {
      std::cout << "\n\n\t crystal-vectors : \n";
      for (int l = 0; l < crystal_harmonics_expansion::get_size(); l++) {
        std::cout << "\t" << l << "\t";
        math::util::print(crystal_harmonics_expansion::get_elements()[l]);
        std::cout << std::endl;
      }
    }
  }

  {
    psi_k.reset();
    crystal_harmonics.reset();

    const std::complex<double> I(0, 1);

    for (int l = 0; l < crystal_harmonics_expansion::get_size(); l++) {
      std::vector<double> r_vec = crystal_harmonics_expansion::get_elements()[l];

      for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
        psi_k(k_ind, l) =
            std::exp(I * math::util::innerProduct(r_vec, k_HOST_VERTEX::get_elements()[k_ind])) /
            std::sqrt(double(k_HOST_VERTEX::dmn_size()));

      for (int w_ind = 0; w_ind < w_VERTEX::dmn_size(); w_ind++)
        for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
          for (int m_ind = 0; m_ind < b::dmn_size(); m_ind++)
            for (int n_ind = 0; n_ind < b::dmn_size(); n_ind++)
              crystal_harmonics(n_ind, m_ind, k_ind, w_ind, n_ind, m_ind, l, w_ind) = psi_k(k_ind, l);
    }
  }

  {
    leading_symmetry_functions.reset();
    leading_symmetry_decomposition.reset();

    for (int l = 0; l < cubic_harmonics_dmn_type::dmn_size(); l++) {
      std::complex<double> norm = 0;

      for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++) {
        std::complex<double> value = 0;

        double kx = k_HOST_VERTEX::get_elements()[k_ind][0];
        double ky = k_HOST_VERTEX::get_elements()[k_ind][1];

        switch (l) {
          case 0:  // s-wave
            value = 1;
            break;

          case 1:  // p-wave
            value = cos(kx) + cos(ky);
            break;

          case 2:  // d-wave
            value = cos(kx) - cos(ky);
            break;

          default:
            throw std::logic_error(__FUNCTION__);
        }

        norm += conj(value) * value;

        leading_symmetry_functions(k_ind, l) = value;
      }

      for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
        leading_symmetry_functions(k_ind, l) /= std::sqrt(1.e-16 + real(norm));
    }
  }
}

// Peter's modify May 4
template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::compute_chi_0_lattice(
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0) {
  profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << std::endl << std::endl;

  using HTS_solver_type =
      DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, dca::linalg::CPU, parameters_type, MOMS_type>;

  using lattice_map_sp_type = DCA::lattice_mapping_sp<parameters_type, k_DCA, k_HOST>;

  using coarsegraining_sp_type = DCA::coarsegraining_sp<parameters_type, k_DCA>;
  using coarsegraining_tp_type = DCA::coarsegraining_tp<parameters_type, k_HOST_VERTEX>;

  lattice_map_sp_type lattice_mapping_obj(parameters);

  MOMS.Sigma_lattice_interpolated = 0.;
  MOMS.Sigma_lattice_coarsegrained = 0.;

  MOMS.Sigma_lattice = 0.;

  if (parameters.do_DCA_plus()) {
    // in case we do the analysis with the DCA+

    if (parameters.use_HTS_approximation()) {
      coarsegraining_sp_type coarsegraining_sp_obj(parameters);

      MOMS_type MOMS_HTS(parameters);

      MOMS_HTS.H_HOST = MOMS.H_HOST;
      MOMS_HTS.H_interactions = MOMS.H_interactions;
      HTS_solver_type HTS_solver(parameters, MOMS_HTS);

      lattice_mapping_obj.execute_with_HTS_approximation(
          MOMS_HTS, HTS_solver, coarsegraining_sp_obj, MOMS.Sigma, MOMS.Sigma_lattice_interpolated,
          MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
    }
    else {
      lattice_mapping_obj.execute(MOMS.Sigma, MOMS.Sigma_lattice_interpolated,
                                  MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
    }

    {
      coarsegraining_tp_type coarsegraining_tp_obj(parameters);

      coarsegraining_tp_obj.execute(MOMS.H_HOST, MOMS.Sigma_lattice, chi_0_function);

      set_chi_0_matrix(chi_0);
    }
  }
  else {
    // in case we do the analysis with the DCA
    coarsegraining_tp_type coarsegraining_tp_obj(parameters);
    coarsegraining_tp_obj.execute(MOMS.H_HOST, MOMS.Sigma, chi_0_function);
    set_chi_0_matrix(chi_0);
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::set_chi_0_matrix(
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0) {
  scalartype renorm = 1. / (parameters.get_beta() * k_HOST_VERTEX::dmn_size());

  for (int w_ind = 0; w_ind < w_VERTEX::dmn_size(); w_ind++)
    for (int K_ind = 0; K_ind < k_HOST_VERTEX::dmn_size(); K_ind++)

      for (int m2 = 0; m2 < b::dmn_size(); m2++)
        for (int n2 = 0; n2 < b::dmn_size(); n2++)

          for (int m1 = 0; m1 < b::dmn_size(); m1++)
            for (int n1 = 0; n1 < b::dmn_size(); n1++)
              chi_0(n1, m1, K_ind, w_ind, n2, m2, K_ind, w_ind) =
                  renorm * chi_0_function(n1, m1, n2, m2, K_ind, w_ind);
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::compute_Gamma_lattice_3(
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& Gamma_cluster) {
  profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

  // Need following condition lines for regular DCA; otherwise execute(Gamma_cluster, Gamma_lattice)
  // still do SVD etc.
  if (!parameters.do_DCA_plus()) {
    int N = lattice_eigenvector_dmn_t::dmn_size();
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        Gamma_lattice(i, j) = Gamma_cluster(i, j);
  }
  else {
    lattice_mapping_tp<parameters_type, k_DCA, k_HOST_VERTEX> lattice_map_tp_obj(parameters);

    lattice_map_tp_obj.execute(Gamma_cluster, Gamma_lattice);
  }

  if (parameters.do_symmetrization_of_Gamma()) {
    if (true) {
      if (concurrency.id() == concurrency.first())
        std::cout << "symmetrize Gamma_lattice according to the symmetry-group \n" << std::endl;

      symmetrize::execute(Gamma_lattice, parameters.get_q_channel_vec());
    }

    if (true) {
      if (concurrency.id() == concurrency.first())
        std::cout << "symmetrize Gamma_lattice according to diagrammatic symmetries \n"
                  << std::endl;

      diagrammatic_symmetries<parameters_type> diagrammatic_symmetries_obj(parameters);

      diagrammatic_symmetries_obj.execute(Gamma_lattice);
    }
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_Gamma_chi_0(
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice) {
  if (parameters.do_diagonalization_on_folded_Gamma_chi_0()) {
    diagonolize_folded_Gamma_chi_0(Gamma_lattice, chi_0_lattice);
  }
  else {
    diagonolize_full_Gamma_chi_0(Gamma_lattice, chi_0_lattice);
  }

  characterize_leading_eigenvectors();
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_full_Gamma_chi_0_real(
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& /*Gamma_lattice*/,
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice) {
  profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

  int N = lattice_eigenvector_dmn_t::dmn_size();

  dca::linalg::Matrix<scalartype, dca::linalg::CPU> chi_0_Gamma_chi_0("chi_0_Gamma_chi_0", N);
  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> chi_0_Gamma_chi_0_temp(
      "chi_0_Gamma_chi_0_temp", N);
  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> chi_0_Gamma("chi_0_Gamma", N);
  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> Gamma("Gamma", N);
  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> chi_0("chi_0", N);

  {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t compute Gamma_chi_0_lattice " << dca::util::print_time() << " ...";

    for (int j = 0; j < N; j++)
      for (int i = 0; i < N; i++)
        Gamma(i, j) = Gamma_sym(i, j);

    for (int j = 0; j < N; j++)
      for (int i = 0; i < N; i++)
        chi_0(i, j) = sqrt(chi_0_lattice(i, j));

    dca::linalg::matrixop::gemm(chi_0, Gamma, chi_0_Gamma);
    dca::linalg::matrixop::gemm(chi_0_Gamma, chi_0, chi_0_Gamma_chi_0_temp);

    for (int j = 0; j < N; j++)
      for (int i = 0; i < N; i++)
        chi_0_Gamma_chi_0(i, j) = real(chi_0_Gamma_chi_0_temp(i, j));

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";
  }

  {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t diagonalize Gamma_chi_0 in pp SC channel " << dca::util::print_time()
                << " ...";

    dca::linalg::Vector<scalartype, dca::linalg::CPU> L("L (BSE_lattice_solver)", N);

    dca::linalg::Matrix<scalartype, dca::linalg::CPU> VR("VR (BSE_lattice_solver)", N);

    dca::linalg::matrixop::eigensolverHermitian('V', 'U', chi_0_Gamma_chi_0, L, VR);

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";

    record_eigenvalues_and_eigenvectors(L, VR);

    print_on_shell_ppSC();

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_full_Gamma_chi_0(
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice) {
  profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

  int N = lattice_eigenvector_dmn_t::dmn_size();

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> Gamma_chi_0("Gamma_chi_0", N);

  {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t compute Gamma_chi_0_lattice " << dca::util::print_time() << " ...";

    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> Gamma("Gamma", N);
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> chi_0("chi_0", N);

    for (int j = 0; j < N; j++)
      for (int i = 0; i < N; i++)
        Gamma(i, j) = Gamma_lattice(i, j);

    for (int j = 0; j < N; j++)
      for (int i = 0; i < N; i++)
        chi_0(i, j) = chi_0_lattice(i, j);

    dca::linalg::matrixop::gemm(Gamma, chi_0, Gamma_chi_0);

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";
  }

  {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t diagonalize Gamma_chi_0 " << dca::util::print_time() << " ...";

    dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU> L("L (BSE_lattice_solver)", N);

    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> VR("VR (BSE_lattice_solver)", N);
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> VL("VL (BSE_lattice_solver)", N);

    dca::linalg::matrixop::eigensolver('N', 'V', Gamma_chi_0, L, VL, VR);

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";

    record_eigenvalues_and_eigenvectors(L, VL, VR);

    print_on_shell();

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::record_eigenvalues_and_eigenvectors(
    dca::linalg::Vector<scalartype, dca::linalg::CPU>& L,
    dca::linalg::Matrix<scalartype, dca::linalg::CPU>& VR) {
  int N = lattice_eigenvector_dmn_t::dmn_size();
  std::vector<std::pair<scalartype, int>> eigenvals_mod(N);

  for (int i = 0; i < N; i++) {
    eigenvals_mod[i].first = std::abs(L[i] - 1.);
    eigenvals_mod[i].second = i;
  }

  std::stable_sort(eigenvals_mod.begin(), eigenvals_mod.end(), math::util::pairLess<scalartype, int>);

  for (int i = 0; i < N_LAMBDAS; i++) {
    int index = eigenvals_mod[i].second;

    leading_eigenvalues_real(i) = L[index];

    for (int j = 0; j < N; j++)
      leading_eigenvectors_real(i, j) = VR(j, index);
  }

  if (concurrency.id() == concurrency.last())
    std::cout << "\n\n\t recording eigenvalues and eigenvectors finished! "
              << dca::util::print_time() << "\n";
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::record_eigenvalues_and_eigenvectors(
    dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU>& L,
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& /*VL*/,
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VR) {
  int N = lattice_eigenvector_dmn_t::dmn_size();
  int M = crystal_eigenvector_dmn_t::dmn_size();

  std::vector<std::pair<std::complex<scalartype>, int>> eigenvals(M);

  for (int i = 0; i < M; i++) {
    eigenvals[i].first = L[i];
    eigenvals[i].second = i;
  }

  stable_sort(eigenvals.begin(), eigenvals.end(),
              math::util::susceptibilityPairGreater<scalartype, int>);

  for (int i = 0; i < N_LAMBDAS; i++) {
    int index = eigenvals[eigenvals.size() - 1 - i].second;

    leading_eigenvalues(i) = L[index];

    for (int j = 0; j < N; j++)
      leading_eigenvectors(i, j) = VR(j, index);
  }

  symmetrize_leading_eigenvectors();
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::diagonolize_folded_Gamma_chi_0(
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& Gamma_lattice,
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice) {
  if (concurrency.id() == concurrency.last())
    std::cout << __FUNCTION__ << std::endl;

  profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);
  if (false) {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t symmetrize Gamma_lattice for even-frequency part (only need to "
                   "symmetrize w2 argument in Gamma(k1,w1,k2,w2)), then compute "
                   "Gamma_chi_0_lattice "
                << dca::util::print_time() << " ...";

    for (int w2 = 0; w2 < w_VERTEX::dmn_size(); w2++)
      for (int K2 = 0; K2 < k_HOST_VERTEX::dmn_size(); K2++)
        for (int m2 = 0; m2 < b::dmn_size(); m2++)
          for (int n2 = 0; n2 < b::dmn_size(); n2++)

            for (int w1 = 0; w1 < w_VERTEX::dmn_size(); w1++)
              for (int K1 = 0; K1 < k_HOST_VERTEX::dmn_size(); K1++)
                for (int m1 = 0; m1 < b::dmn_size(); m1++)
                  for (int n1 = 0; n1 < b::dmn_size(); n1++) {
                    Gamma_sym(n1, m1, K1, w1, n2, m2, K2, w2) =
                        0.5 *
                        (Gamma_lattice(n1, m1, K1, w1, n2, m2, K2, w2) +
                         Gamma_lattice(n1, m1, K1, w1, n2, m2, K2, w_VERTEX::dmn_size() - 1 - w2));
                  }
  }

  int N = lattice_eigenvector_dmn_t::dmn_size();
  int M = crystal_eigenvector_dmn_t::dmn_size();

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> Gamma_chi_0_crystal("Gamma_chi_0",
                                                                                      M);

  {
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> Gamma_chi_0_lattice(
        "Gamma_chi_0", N);

    {
      if (concurrency.id() == concurrency.last())
        std::cout << "\n\n\t compute Gamma_chi_0_lattice " << dca::util::print_time() << " ...";

      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> Gamma("Gamma", N);
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> chi_0("chi_0", N);

      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          Gamma(i, j) = Gamma_lattice(i, j);

      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          chi_0(i, j) = chi_0_lattice(i, j);

      dca::linalg::matrixop::gemm(Gamma, chi_0, Gamma_chi_0_lattice);

      if (concurrency.id() == concurrency.last())
        std::cout << " finished " << dca::util::print_time() << "\n";
    }

    {
      if (concurrency.id() == concurrency.last())
        std::cout << "\n\n\t compute P_Gamma_chi_0_lattice_P " << dca::util::print_time() << " ...";

      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> P("P",
                                                                        std::pair<int, int>(N, M));
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> tmp(
          "tmp", std::pair<int, int>(N, M));

      for (int j = 0; j < M; j++)
        for (int i = 0; i < N; i++)
          P(i, j) = crystal_harmonics(i, j);

      dca::linalg::matrixop::gemm('N', 'N', Gamma_chi_0_lattice, P, tmp);
      dca::linalg::matrixop::gemm('C', 'N', P, tmp, Gamma_chi_0_crystal);

      if (concurrency.id() == concurrency.last())
        std::cout << " finished " << dca::util::print_time() << "\n";
    }
  }

  {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t diagonalize P_Gamma_chi_0_lattice_P " << dca::util::print_time()
                << " ...";

    dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU> L("L (BSE_lattice_solver)", M);

    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> VL("VL (BSE_lattice_solver)", M);
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> VR("VR (BSE_lattice_solver)", M);

    dca::linalg::matrixop::eigensolver('N', 'V', Gamma_chi_0_crystal, L, VL, VR);

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";

    record_eigenvalues_and_folded_eigenvectors(L, VL, VR);

    print_on_shell();

    if (parameters.compute_P_q_lattice()) {
    }
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::record_eigenvalues_and_folded_eigenvectors(
    dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU>& L,
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& /*VL*/,
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VR) {
  int N = lattice_eigenvector_dmn_t::dmn_size();
  int M = crystal_eigenvector_dmn_t::dmn_size();

  std::vector<std::pair<std::complex<scalartype>, int>> eigenvals(M);

  for (int i = 0; i < M; i++) {
    eigenvals[i].first = L[i];
    eigenvals[i].second = i;
  }

  stable_sort(eigenvals.begin(), eigenvals.end(),
              math::util::susceptibilityPairGreater<scalartype, int>);

  for (int i = 0; i < N_LAMBDAS; i++) {
    int index = eigenvals[eigenvals.size() - 1 - i].second;

    leading_eigenvalues(i) = L[index];

    for (int j = 0; j < N; j++)
      for (int l = 0; l < M; l++)
        leading_eigenvectors(i, j) += crystal_harmonics(j, l) * VR(l, index);
  }

  symmetrize_leading_eigenvectors();
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::compute_folded_susceptibility(
    func::function<std::complex<scalartype>, HOST_matrix_dmn_t>& chi_0_lattice,
    dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU>& L,
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VL,
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& VR) {
  int N = lattice_eigenvector_dmn_t::dmn_size();
  int M = crystal_eigenvector_dmn_t::dmn_size();

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> P_1_min_Gamma_chi_0_P(
      "P_1_min_Gamma_chi_0_P", M);
  {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t invert VR  " << dca::util::print_time() << " ...";

    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> D_inv("D_inv", M);
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> VR_D_inv("VR_D_inv", M);

    for (int l = 0; l < M; l++)
      D_inv(l, l) = 1. / (1. - L[l]);

    VL = VR;
    dca::linalg::matrixop::inverse(VL);

    dca::linalg::matrixop::gemm('N', 'N', VR, D_inv, VR_D_inv);
    dca::linalg::matrixop::gemm('N', 'N', VR_D_inv, VL, P_1_min_Gamma_chi_0_P);

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";
  }

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> P_chi_0_P("P_chi_0_P", M);
  {
    if (concurrency.id() == concurrency.last())
      std::cout << "\n\n\t compute P_chi_0_P " << dca::util::print_time() << " ...";

    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> chi_0("chi_0", N);

    for (int j = 0; j < N; j++)
      for (int i = 0; i < N; i++)
        chi_0(i, j) = chi_0_lattice(i, j);

    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> P("P", std::pair<int, int>(N, M));
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> tmp("tmp",
                                                                        std::pair<int, int>(N, M));

    for (int j = 0; j < M; j++)
      for (int i = 0; i < N; i++)
        P(i, j) = crystal_harmonics(i, j);

    dca::linalg::matrixop::gemm('N', 'N', chi_0, P, tmp);
    dca::linalg::matrixop::gemm('C', 'N', P, tmp, P_chi_0_P);

    if (concurrency.id() == concurrency.last())
      std::cout << " finished " << dca::util::print_time() << "\n";
  }

  {
    func::function<std::complex<scalartype>,
                   func::dmn_variadic<crystal_eigenvector_dmn_t, crystal_eigenvector_dmn_t>>
        chi_q_tmp("chi_q_tmp");

    {
      dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> chi_matrix("chi", M);

      dca::linalg::matrixop::gemm('N', 'N', P_chi_0_P, P_1_min_Gamma_chi_0_P, chi_matrix);

      for (int j = 0; j < M; j++)
        for (int i = 0; i < M; i++)
          chi_q_tmp(i, j) = chi_matrix(i, j);
    }

    chi_q = 0.;

    for (int w2 = 0; w2 < w_VERTEX::dmn_size(); w2++)
      for (int K2 = 0; K2 < M; K2++)
        for (int m2 = 0; m2 < b::dmn_size(); m2++)
          for (int n2 = 0; n2 < b::dmn_size(); n2++)

            for (int w1 = 0; w1 < w_VERTEX::dmn_size(); w1++)
              for (int K1 = 0; K1 < M; K1++)
                for (int m1 = 0; m1 < b::dmn_size(); m1++)
                  for (int n1 = 0; n1 < b::dmn_size(); n1++)
                    chi_q(n1, m1, K1, n2, m2, K2) += chi_q_tmp(n1, m1, K1, w1, n2, m2, K2, w2);
  }

  if (concurrency.id() == concurrency.last()) {
    std::cout << "\n\n\t real(chi_q) \n\n";
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < M; j++)
        std::cout << real(chi_q(i, j)) << "\t";
      std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "\n\n\t imag(chi_q) \n\n";
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < M; j++)
        std::cout << imag(chi_q(i, j)) << "\t";
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::print_on_shell() {
  if (concurrency.id() == concurrency.last())
    std::cout << __FUNCTION__ << std::endl;

  int N = k_HOST_VERTEX::dmn_size();
  int M = crystal_harmonics_expansion_dmn_t::dmn_size();

  if (concurrency.id() == concurrency.last()) {
    std::cout.precision(6);
    std::cout << std::scientific;

    {
      std::cout << "\n\n\t\t leading eigenvalues : ( T=" << 1. / parameters.get_beta() << " )\n\n";
      for (int i = 0; i < N_LAMBDAS; i++)
        std::cout << "\t" << i << "\t[" << real(leading_eigenvalues(i)) << ", "
                  << imag(leading_eigenvalues(i)) << "]\n";
    }

    {
      std::cout << "\n\n\t\t leading eigenvectors : \n\n";

      {
        std::cout << "\t" << -1 << "\t[ " << 0. << ", " << 0. << "]";
        for (int i = 0; i < N_LAMBDAS; i++)
          std::cout << "\t[ " << real(leading_eigenvalues(i)) << ", "
                    << imag(leading_eigenvalues(i)) << "]";
        std::cout << "\n\t========================================\n";
      }

      for (int l = 0; l < M; l++) {
        std::cout << "\t" << l << "\t[ " << crystal_harmonics_expansion::get_elements()[l][0]
                  << ", " << crystal_harmonics_expansion::get_elements()[l][1] << "]";

        for (int i = 0; i < N_LAMBDAS; i++) {
          std::complex<scalartype> result = 0;
          std::complex<scalartype> norm = 0;

          for (int j = 0; j < N; j++) {
            result += conj(psi_k(j, l)) * leading_eigenvectors(i, 0, 0, j, w_VERTEX::dmn_size() / 2);
            norm += conj(leading_eigenvectors(i, 0, 0, j, w_VERTEX::dmn_size() / 2)) *
                    leading_eigenvectors(i, 0, 0, j, w_VERTEX::dmn_size() / 2);
          }

          std::cout << "\t[ " << real(result / std::sqrt(real(norm) + 1.e-16)) << ", "
                    << imag(result / std::sqrt(real(norm) + 1.e-16)) << "]";
        }

        std::cout << "\n";
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::print_on_shell_ppSC() {
  if (concurrency.id() == concurrency.last())
    std::cout << __FUNCTION__ << std::endl;

  int N = k_HOST_VERTEX::dmn_size();
  //  std::cout << "\t N=" << k_HOST_VERTEX::dmn_size() << "\n";

  if (concurrency.id() == concurrency.last()) {
    std::cout.precision(6);
    std::cout << std::scientific;
    std::cout << "\n\n\t\t 10 leading eigenvalues : ( T=" << 1. / parameters.get_beta() << " )\n\n";
    for (int i = 0; i < N_LAMBDAS; i++)
      std::cout << "\t" << i << "\t[" << leading_eigenvalues_real(i) << "]\n";

    {
      int ind0pi = 0;
      int indpi0 = 0;
      // record the index for k=[0,pi] and [pi,0]
      for (int k_ind = 0; k_ind < N; k_ind++) {
        double kx = k_HOST_VERTEX::get_elements()[k_ind][0];
        double ky = k_HOST_VERTEX::get_elements()[k_ind][1];
        if (std::abs(kx) < 1.e-3 && std::abs(ky - 3.141592) < 1.e-3) {
          ind0pi = k_ind;
        }
        if (std::abs(ky) < 1.e-3 && std::abs(kx - 3.141592) < 1.e-3) {
          indpi0 = k_ind;
        }
      }

      std::cout << "\n\n\t\t Phi_(:,k) for 10 leading eigenvectors: \n\n";
      std::cout << "  w         k=[" << k_HOST_VERTEX::get_elements()[ind0pi][0] << ", "
                << k_HOST_VERTEX::get_elements()[ind0pi][1] << "]         k=["
                << k_HOST_VERTEX::get_elements()[indpi0][0] << ", "
                << k_HOST_VERTEX::get_elements()[indpi0][1] << "] \n\n";
      for (int i = 0; i < N_LAMBDAS; i++) {
        for (int w = 0; w < w_VERTEX::dmn_size(); w++) {
          std::cout << i << "   " << w_VERTEX::get_elements()[w] << "   "
                    << leading_eigenvectors_real(i, 0, 0, ind0pi, w) << "   "
                    << leading_eigenvectors_real(i, 0, 0, indpi0, w) << "\n";
        }
        std::cout << "----------------------------------------------------------------------"
                  << "\n";
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::symmetrize_leading_eigenvectors() {
  profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t" << __FUNCTION__ << " " << dca::util::print_time() << std::endl;

  int N = lattice_eigenvector_dmn_t::dmn_size();
  ;

  for (int i = 0; i < N_LAMBDAS; i++) {
    scalartype N_phases = 1.e4;

    std::complex<scalartype> alpha_min = 0;
    scalartype norm = 1.e6;

    for (int l = 0; l < N_phases; l++) {
      std::complex<scalartype> alpha =
          std::complex<scalartype>(cos(2. * M_PI * l / N_phases), sin(2. * M_PI * l / N_phases));

      scalartype result = 0;

      for (int w1 = 0; w1 < w_VERTEX::dmn_size() / 2; w1++)
        for (int K1 = 0; K1 < k_HOST_VERTEX::dmn_size(); K1++)
          for (int m1 = 0; m1 < b::dmn_size(); m1++)
            for (int n1 = 0; n1 < b::dmn_size(); n1++)
              result += std::abs(
                  alpha * leading_eigenvectors(i, n1, m1, K1, w1) -
                  conj(alpha * leading_eigenvectors(i, n1, m1, K1, w_VERTEX::dmn_size() - 1 - w1)));

      if (result < norm) {
        norm = result;
        alpha_min = alpha;
      }
    }

    for (int l = 0; l < N; l++)
      leading_eigenvectors(i, l) *= alpha_min;
  }

  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << " finished " << dca::util::print_time() << std::endl;
}

template <class parameters_type, class MOMS_type>
void BSE_lattice_solver<parameters_type, MOMS_type>::characterize_leading_eigenvectors() {
  profiler_type prof(__FUNCTION__, "BSE_lattice_solver", __LINE__);

  for (int n_lam = 0; n_lam < N_LAMBDAS; n_lam++) {
    for (int w_ind = 0; w_ind < w_VERTEX::dmn_size(); w_ind++) {
      for (int m_ind = 0; m_ind < b::dmn_size(); m_ind++) {
        for (int n_ind = 0; n_ind < b::dmn_size(); n_ind++) {
          for (int n_cub = 0; n_cub < N_CUBIC; n_cub++) {
            std::complex<scalartype> scal_prod = 0;

            std::complex<scalartype> norm_phi = 0;
            std::complex<scalartype> norm_psi = 0;

            for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++) {
              std::complex<scalartype> phi_k_val =
                  leading_eigenvectors(n_lam, m_ind, n_ind, k_ind, w_ind);
              std::complex<scalartype> psi_k_val = leading_symmetry_functions(k_ind, n_cub);

              scal_prod += conj(psi_k_val) * phi_k_val;

              norm_phi += conj(phi_k_val) * phi_k_val;
              norm_psi += conj(psi_k_val) * psi_k_val;
            }

            leading_symmetry_decomposition(n_lam, m_ind, n_ind, n_cub, w_ind) =
                scal_prod / std::sqrt(1.e-16 + norm_phi * norm_psi);
          }
        }
      }
    }
  }
}

}  // DCA

#endif  // PHYS_LIBRARY_DCA_ANALYSIS_BSE_SOLVER_BSE_LATTICE_SOLVER_H
