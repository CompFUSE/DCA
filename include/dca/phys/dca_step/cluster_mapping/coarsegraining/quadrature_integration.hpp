// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// Quadrature integration class.
//
// Template parameters:
// IntegrationDmn: Domain/space of the integration.
// OtherDmn: Other dimensions of the considered functions.
// Threading: Threading library for parallelization over the integration domain.
//
// TODO: Rename class and methods as there is no integration done here, but Green's functions
//       computed.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_QUADRATURE_INTEGRATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_QUADRATURE_INTEGRATION_HPP

#include <cmath>  // std::exp, std::isnan
#include <complex>
#include <iostream>
#include <stdexcept>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/parallel/util/get_bounds.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename ScalarType, typename IntegrationDmn, typename OtherDmn, typename Threading = void>
class quadrature_integration {
private:
  using Function =
      func::function<std::complex<ScalarType>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>;

public:
  static void quadrature_integration_G_q_w_st(const Function& I_q, const Function& H_q,
                                              const Function& S_q, Function& G_q);

  static void quadrature_integration_G_q_w_mt(int nr_threads, const Function& I_q,
                                              const Function& H_q, const Function& S_q,
                                              Function& G_q);

  static void quadrature_integration_G_q_t_st(ScalarType beta, ScalarType f_val, ScalarType t_val,
                                              const Function& I_q, const Function& H_q,
                                              Function& G_q);

  static void quadrature_integration_G_q_t_mt(int nr_threads, ScalarType beta, ScalarType f_val,
                                              ScalarType t_val, const Function& I_q,
                                              const Function& H_q, Function& G_q);

private:
  static void quadrature_integration_G_q_w_mt_impl(int id, int nr_threads, const Function& I_q,
                                                   const Function& H_q, const Function& S_q,
                                                   Function& G_q);

  static void quadrature_integration_G_q_t_mt_impl(int id, int nr_threads, ScalarType beta,
                                                   ScalarType f_val, ScalarType t_val,
                                                   const Function& I_q, const Function& H_q,
                                                   Function& G_q);
};

template <typename ScalarType, typename IntegrationDmn, typename OtherDmn, typename Threading>
void quadrature_integration<ScalarType, IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_w_st(
    const Function& I_q, const Function& H_q, const Function& S_q, Function& G_q) {
  linalg::Matrix<std::complex<ScalarType>, linalg::CPU> G_inv("G_inv", OtherDmn::dmn_size());

  // Allocate the work space for inverse only once.
  linalg::Vector<int, linalg::CPU> ipiv;
  linalg::Vector<std::complex<ScalarType>, linalg::CPU> work;

  for (int q_ind = 0; q_ind < IntegrationDmn::dmn_size(); ++q_ind) {
    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        G_inv(i, j) = I_q(i, j, q_ind) - H_q(i, j, q_ind) - S_q(i, j, q_ind);

    linalg::matrixop::inverse(G_inv, ipiv, work);

    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        G_q(i, j, q_ind) = G_inv(i, j);
  }
}

template <typename ScalarType, typename IntegrationDmn, typename OtherDmn, typename Threading>
void quadrature_integration<ScalarType, IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_w_mt(
    const int nr_threads, const Function& I_q, const Function& H_q, const Function& S_q,
    Function& G_q) {
  G_q = 0.;

  Threading parallelization_obj;
  parallelization_obj.execute(nr_threads, quadrature_integration_G_q_w_mt_impl, std::ref(I_q),
                              std::ref(H_q), std::ref(S_q), std::ref(G_q));
}

template <typename ScalarType, typename IntegrationDmn, typename OtherDmn, typename Threading>
void quadrature_integration<ScalarType, IntegrationDmn, OtherDmn,
                            Threading>::quadrature_integration_G_q_w_mt_impl(int id, int nr_threads,
                                                                             const Function& I_q,
                                                                             const Function& H_q,
                                                                             const Function& S_q,
                                                                             Function& G_q) {
  const IntegrationDmn q_dmn;
  const std::pair<int, int> q_bounds = parallel::util::getBounds(id, nr_threads, q_dmn);

  linalg::Matrix<std::complex<ScalarType>, linalg::CPU> G_inv("G_inv", OtherDmn::dmn_size());

  // Allocate the work space for inverse only once.
  linalg::Vector<int, linalg::CPU> ipiv;
  linalg::Vector<std::complex<ScalarType>, linalg::CPU> work;

  for (int q_ind = q_bounds.first; q_ind < q_bounds.second; ++q_ind) {
    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        G_inv(i, j) = I_q(i, j, q_ind) - H_q(i, j, q_ind) - S_q(i, j, q_ind);

    linalg::matrixop::inverse(G_inv, ipiv, work);

    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        G_q(i, j, q_ind) = G_inv(i, j);
  }
}

template <typename ScalarType, typename IntegrationDmn, typename OtherDmn, typename Threading>
void quadrature_integration<ScalarType, IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_t_st(
    const ScalarType beta, const ScalarType f_val, const ScalarType t_val, const Function& I_q,
    const Function& H_q, Function& G_q) {
  G_q = 0.;

  linalg::Matrix<std::complex<ScalarType>, linalg::CPU> H_m("H_m", OtherDmn::dmn_size());

  linalg::Vector<ScalarType, linalg::CPU> L("e_l", OtherDmn::dmn_size());
  linalg::Matrix<std::complex<ScalarType>, linalg::CPU> V("V_l", OtherDmn::dmn_size());

  linalg::Vector<ScalarType, linalg::CPU> G_t("e_l", OtherDmn::dmn_size());

  for (int q_ind = 0; q_ind < IntegrationDmn::dmn_size(); ++q_ind) {
    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        H_m(i, j) = H_q(i, j, q_ind) - I_q(i, j, q_ind);

    linalg::matrixop::eigensolverGreensFunctionMatrix('V', 'U', H_m, L, V);

    for (int i = 0; i < OtherDmn::dmn_size(); ++i) {
      if (L[i] < 0)
        G_t[i] = f_val * std::exp(L[i] * (beta - t_val)) / (std::exp(L[i] * beta) + 1.);
      else
        G_t[i] = f_val * std::exp(-L[i] * t_val) / (std::exp(-L[i] * beta) + 1.);

      if (std::isnan(G_t[i])) {
        std::cout << "\n"
                  << "\tWarning in " << __FUNCTION__ << ": G_t[i] = " << G_t[i] << "\n"
                  << "\n"
                  << "\tL[i] = " << L[i] << "\n"
                  << "\tbeta = " << beta << "\n"
                  << "\ttau = " << t_val << "\n"
                  << "\tstd::exp(L[i]*beta) = " << std::exp(L[i] * beta) << "\n"
                  << "\tstd::exp(L[i]*(beta-t_val)) = " << std::exp(L[i] * (beta - t_val)) << "\n"
                  << std::endl;

        throw std::logic_error(__FUNCTION__);
      }
    }
    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        for (int l = 0; l < OtherDmn::dmn_size(); ++l)
          G_q(i, j, q_ind) += V(i, l) * G_t[l] * conj(V(j, l));  // G_t[l]*real(conj(V(l,i))*V(l,j));
  }
}

template <typename ScalarType, typename IntegrationDmn, typename OtherDmn, typename Threading>
void quadrature_integration<ScalarType, IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_t_mt(
    const int nr_threads, const ScalarType beta, const ScalarType f_val, const ScalarType t_val,
    const Function& I_q, const Function& H_q, Function& G_q) {
  G_q = 0.;

  Threading parallelization_obj;
  parallelization_obj.execute(nr_threads, quadrature_integration_G_q_t_mt_impl, beta, f_val, t_val,
                              std::ref(I_q), std::ref(H_q), std::ref(G_q));
}

template <typename ScalarType, typename IntegrationDmn, typename OtherDmn, typename Threading>
void quadrature_integration<ScalarType, IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_t_mt_impl(
    const int id, const int nr_threads, const ScalarType beta, const ScalarType f_val,
    const ScalarType t_val, const Function& I_q, const Function& H_q, Function& G_q) {
  const IntegrationDmn q_dmn;
  const std::pair<int, int> q_bounds = parallel::util::getBounds(id, nr_threads, q_dmn);

  linalg::Matrix<std::complex<ScalarType>, linalg::CPU> H_m("H_m", OtherDmn::dmn_size());

  linalg::Vector<ScalarType, linalg::CPU> L("e_l", OtherDmn::dmn_size());
  linalg::Matrix<std::complex<ScalarType>, linalg::CPU> V("V_l", OtherDmn::dmn_size());

  linalg::Vector<ScalarType, linalg::CPU> G_t("e_l", OtherDmn::dmn_size());

  for (int q_ind = q_bounds.first; q_ind < q_bounds.second; ++q_ind) {
    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        H_m(i, j) = H_q(i, j, q_ind) - I_q(i, j, q_ind);

    linalg::matrixop::eigensolverGreensFunctionMatrix('V', 'U', H_m, L, V);

    for (int i = 0; i < OtherDmn::dmn_size(); ++i) {
      if (L[i] < 0)
        G_t[i] = f_val * std::exp(L[i] * (beta - t_val)) / (std::exp(L[i] * beta) + 1.);
      else
        G_t[i] = f_val * std::exp(-L[i] * t_val) / (std::exp(-L[i] * beta) + 1.);

      if (std::isnan(G_t[i])) {
        std::cout << "\n"
                  << "\tWarning in " << __FUNCTION__ << ": G_t[i] = " << G_t[i] << "\n"
                  << "\n"
                  << "\tL[i] = " << L[i] << "\n"
                  << "\tbeta = " << beta << "\n"
                  << "\ttau = " << t_val << "\n"
                  << "\tstd::exp(L[i]*beta) = " << std::exp(L[i] * beta) << "\n"
                  << "\tstd::exp(L[i]*(beta-t_val)) = " << std::exp(L[i] * (beta - t_val)) << "\n"
                  << std::endl;

        throw std::logic_error(__FUNCTION__);
      }
    }

    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < OtherDmn::dmn_size(); ++i)
        for (int l = 0; l < OtherDmn::dmn_size(); ++l)
          G_q(i, j, q_ind) += V(i, l) * G_t[l] * conj(V(j, l));
  }
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_QUADRATURE_INTEGRATION_HPP
