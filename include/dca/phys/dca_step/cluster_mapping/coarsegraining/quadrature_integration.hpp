// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
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

#include <complex>
#include <iostream>
#include <stdexcept>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/parallel/util/get_bounds.hpp"
#include "dca/parallel/util/threading_data.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename IntegrationDmn, typename OtherDmn, typename Threading = void>
class quadrature_integration {
public:
  template <typename scalar_type>
  static void quadrature_integration_G_q_w_st(
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& S_q,
      func::function<std::complex<scalar_type>,
                     func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q);

  template <typename scalar_type>
  static void quadrature_integration_G_q_w_mt(
      int nr_threads, func::function<std::complex<scalar_type>,
                                     func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& S_q,
      func::function<std::complex<scalar_type>,
                     func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q);

  template <typename scalar_type>
  static void quadrature_integration_G_q_t_st(
      scalar_type beta, scalar_type f_val, scalar_type t_val,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
      func::function<std::complex<scalar_type>,
                     func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q);

  template <typename scalar_type>
  static void quadrature_integration_G_q_t_mt(
      int nr_threads, scalar_type beta, scalar_type f_val, scalar_type t_val,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
      func::function<std::complex<scalar_type>,
                     func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q);

private:
  template <typename scalar_type>
  static void* quadrature_integration_G_q_w_mt(void* data);

  template <typename scalar_type>
  static void* quadrature_integration_G_q_t_mt(void* data);

  template <typename scalar_type>
  struct quadrature_integration_functions {
    quadrature_integration_functions()
        : I_q_ptr(nullptr), H_q_ptr(nullptr), S_q_ptr(nullptr), G_q_ptr(nullptr) {}

    scalar_type beta;

    scalar_type f_val;
    scalar_type t_val;

    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>* I_q_ptr;
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>* H_q_ptr;
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>* S_q_ptr;
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>* G_q_ptr;
  };
};

template <typename IntegrationDmn, typename OtherDmn, typename Threading>
template <typename scalar_type>
void quadrature_integration<IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_w_st(
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& S_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q) {
  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> G_inv("G_inv",
                                                                         OtherDmn::dmn_size());

  // Allocate the work space for inverse only once.
  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU> work;

  for (int q_ind = 0; q_ind < IntegrationDmn::dmn_size(); q_ind++) {
    for (int j = 0; j < OtherDmn::dmn_size(); j++)
      for (int i = 0; i < OtherDmn::dmn_size(); i++)
        G_inv(i, j) = I_q(i, j, q_ind) - H_q(i, j, q_ind) - S_q(i, j, q_ind);

    dca::linalg::matrixop::inverse(G_inv, ipiv, work);

    for (int j = 0; j < OtherDmn::dmn_size(); j++)
      for (int i = 0; i < OtherDmn::dmn_size(); i++)
        G_q(i, j, q_ind) = G_inv(i, j);
  }
}

template <typename IntegrationDmn, typename OtherDmn, typename Threading>
template <typename scalar_type>
void quadrature_integration<IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_w_mt(
    int nr_threads,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& S_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q) {
  G_q = 0.;

  quadrature_integration_functions<scalar_type> quadrature_integration_functions_obj;

  quadrature_integration_functions_obj.I_q_ptr = &I_q;
  quadrature_integration_functions_obj.H_q_ptr = &H_q;
  quadrature_integration_functions_obj.S_q_ptr = &S_q;
  quadrature_integration_functions_obj.G_q_ptr = &G_q;

  Threading parallelization_obj;

  parallelization_obj.execute(nr_threads, quadrature_integration_G_q_w_mt<scalar_type>,
                              (void*)&quadrature_integration_functions_obj);
}

template <typename IntegrationDmn, typename OtherDmn, typename Threading>
template <typename scalar_type>
void* quadrature_integration<IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_w_mt(
    void* void_ptr) {
  typedef quadrature_integration_functions<scalar_type> quadrature_functions_type;

  dca::parallel::ThreadingData* data_ptr = static_cast<dca::parallel::ThreadingData*>(void_ptr);
  quadrature_functions_type* functions_ptr = static_cast<quadrature_functions_type*>(data_ptr->arg);

  int id = data_ptr->id;
  int nr_threads = data_ptr->num_threads;

  func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q =
      *(functions_ptr->I_q_ptr);
  func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q =
      *(functions_ptr->H_q_ptr);
  func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& S_q =
      *(functions_ptr->S_q_ptr);
  func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q =
      *(functions_ptr->G_q_ptr);

  IntegrationDmn q_dmn;
  std::pair<int, int> q_bounds = dca::parallel::util::getBounds(id, nr_threads, q_dmn);

  {
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> G_inv("G_inv",
                                                                           OtherDmn::dmn_size());

    // Allocate the work space for inverse only once.
    dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
    dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU> work;

    for (int q_ind = q_bounds.first; q_ind < q_bounds.second; q_ind += 1) {
      for (int j = 0; j < OtherDmn::dmn_size(); j++)
        for (int i = 0; i < OtherDmn::dmn_size(); i++)
          G_inv(i, j) = I_q(i, j, q_ind) - H_q(i, j, q_ind) - S_q(i, j, q_ind);

      dca::linalg::matrixop::inverse(G_inv, ipiv, work);

      for (int j = 0; j < OtherDmn::dmn_size(); j++)
        for (int i = 0; i < OtherDmn::dmn_size(); i++)
          G_q(i, j, q_ind) = G_inv(i, j);
    }
  }

  return 0;
}

template <typename IntegrationDmn, typename OtherDmn, typename Threading>
template <typename scalar_type>
void quadrature_integration<IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_t_st(
    scalar_type beta, scalar_type f_val, scalar_type t_val,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q) {
  G_q = 0.;

  {
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> H_m("H_m", OtherDmn::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> L("e_l", OtherDmn::dmn_size());
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> V("V_l", OtherDmn::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> G_t("e_l", OtherDmn::dmn_size());

    for (int q_ind = 0; q_ind < IntegrationDmn::dmn_size(); q_ind++) {
      for (int j = 0; j < OtherDmn::dmn_size(); j++)
        for (int i = 0; i < OtherDmn::dmn_size(); i++)
          H_m(i, j) = H_q(i, j, q_ind) - I_q(i, j, q_ind);

      if (false)
        dca::linalg::matrixop::eigensolverHermitian('V', 'U', H_m, L, V);
      else
        dca::linalg::matrixop::eigensolverGreensFunctionMatrix('V', 'U', H_m, L, V);

      for (int i = 0; i < OtherDmn::dmn_size(); i++) {
        if (L[i] < 0)
          G_t[i] = f_val * std::exp(L[i] * (beta - t_val)) / (std::exp(L[i] * beta) + 1.);
        else
          G_t[i] = f_val * std::exp(-L[i] * t_val) / (std::exp(-L[i] * beta) + 1.);

        if (G_t[i] != G_t[i]) {
          std::cout << "\n\t warning in compute_G_q_t --> G_t[i] : " << G_t[i] << "\n";
          std::cout << "\n\tL[i] : " << L[i] << "\n";
          std::cout << "\n\tbeta : " << beta << "\n";
          std::cout << "\n\ttau  : " << t_val << "\n";
          std::cout << "\n\tstd::exp(L[i]*beta)  : " << std::exp(L[i] * beta) << "\n";
          std::cout << "\n\tstd::exp(L[i]*(beta-t_val)) : " << std::exp(L[i] * (beta - t_val))
                    << "\n";

          throw std::logic_error(__FUNCTION__);
        }
      }
      for (int j = 0; j < OtherDmn::dmn_size(); j++)
        for (int i = 0; i < OtherDmn::dmn_size(); i++)
          for (int l = 0; l < OtherDmn::dmn_size(); l++)
            G_q(i, j, q_ind) += V(i, l) * G_t[l] * conj(V(j, l));  // G_t[l]*real(conj(V(l,i))*V(l,j));
    }
  }
}

template <typename IntegrationDmn, typename OtherDmn, typename Threading>
template <typename scalar_type>
void quadrature_integration<IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_t_mt(
    int nr_threads, scalar_type beta, scalar_type f_val, scalar_type t_val,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q) {
  G_q = 0.;

  quadrature_integration_functions<scalar_type> quadrature_integration_functions_obj;

  quadrature_integration_functions_obj.beta = beta;

  quadrature_integration_functions_obj.t_val = t_val;
  quadrature_integration_functions_obj.f_val = f_val;

  quadrature_integration_functions_obj.I_q_ptr = &I_q;
  quadrature_integration_functions_obj.H_q_ptr = &H_q;
  quadrature_integration_functions_obj.G_q_ptr = &G_q;

  Threading parallelization_obj;

  parallelization_obj.execute(nr_threads, quadrature_integration_G_q_t_mt<scalar_type>,
                              (void*)&quadrature_integration_functions_obj);
}

template <typename IntegrationDmn, typename OtherDmn, typename Threading>
template <typename scalar_type>
void* quadrature_integration<IntegrationDmn, OtherDmn, Threading>::quadrature_integration_G_q_t_mt(
    void* void_ptr) {
  typedef quadrature_integration_functions<scalar_type> quadrature_functions_type;

  dca::parallel::ThreadingData* data_ptr = static_cast<dca::parallel::ThreadingData*>(void_ptr);
  quadrature_functions_type* functions_ptr = static_cast<quadrature_functions_type*>(data_ptr->arg);

  int id = data_ptr->id;
  int nr_threads = data_ptr->num_threads;

  double beta = functions_ptr->beta;

  double t_val = functions_ptr->t_val;
  double f_val = functions_ptr->f_val;

  func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& I_q =
      *(functions_ptr->I_q_ptr);
  func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& H_q =
      *(functions_ptr->H_q_ptr);
  func::function<std::complex<scalar_type>, func::dmn_variadic<OtherDmn, OtherDmn, IntegrationDmn>>& G_q =
      *(functions_ptr->G_q_ptr);

  IntegrationDmn q_dmn;
  std::pair<int, int> q_bounds = dca::parallel::util::getBounds(id, nr_threads, q_dmn);

  {
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> H_m("H_m", OtherDmn::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> L("e_l", OtherDmn::dmn_size());
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> V("V_l", OtherDmn::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> G_t("e_l", OtherDmn::dmn_size());

    for (int q_ind = q_bounds.first; q_ind < q_bounds.second; q_ind += 1) {
      for (int j = 0; j < OtherDmn::dmn_size(); j++)
        for (int i = 0; i < OtherDmn::dmn_size(); i++)
          H_m(i, j) = H_q(i, j, q_ind) - I_q(i, j, q_ind);

      if (false)
        dca::linalg::matrixop::eigensolverHermitian('V', 'U', H_m, L, V);
      else
        dca::linalg::matrixop::eigensolverGreensFunctionMatrix('V', 'U', H_m, L, V);

      for (int i = 0; i < OtherDmn::dmn_size(); i++) {
        if (L[i] < 0)
          G_t[i] = f_val * std::exp(L[i] * (beta - t_val)) / (std::exp(L[i] * beta) + 1.);
        else
          G_t[i] = f_val * std::exp(-L[i] * t_val) / (std::exp(-L[i] * beta) + 1.);

        if (G_t[i] != G_t[i]) {
          std::cout << "\n\t warning in compute_G_q_t --> G_t[i] : " << G_t[i] << "\n";
          std::cout << "\n\tL[i] : " << L[i] << "\n";
          std::cout << "\n\tbeta : " << beta << "\n";
          std::cout << "\n\ttau  : " << t_val << "\n";
          std::cout << "\n\tstd::exp(L[i]*beta)  : " << std::exp(L[i] * beta) << "\n";
          std::cout << "\n\tstd::exp(L[i]*(beta-t_val)) : " << std::exp(L[i] * (beta - t_val))
                    << "\n";

          throw std::logic_error(__FUNCTION__);
        }
      }

      for (int j = 0; j < OtherDmn::dmn_size(); j++)
        for (int i = 0; i < OtherDmn::dmn_size(); i++)
          for (int l = 0; l < OtherDmn::dmn_size(); l++)
            G_q(i, j, q_ind) += V(i, l) * G_t[l] * conj(V(j, l));
    }
  }

  return 0;
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_QUADRATURE_INTEGRATION_HPP
