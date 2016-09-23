// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_QUADRATURE_INTEGRATION_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_QUADRATURE_INTEGRATION_H

#include <complex>
#include <iostream>
#include <stdexcept>

#include "dca/concurrency/util/get_bounds.hpp"
#include "dca/concurrency/parallelization_pthreads.h"
#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg.hpp"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace DCA {

template <typename parameters_type, typename q_dmn_t>
class quadrature_integration {
public:
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;

public:
  template <typename scalar_type>
  static void quadrature_integration_G_q_w_st(
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

  template <typename scalar_type>
  static void quadrature_integration_G_q_w_mt(
      int nr_threads, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

  template <typename scalar_type>
  static void quadrature_integration_G_q_t_st(
      scalar_type beta, scalar_type f_val, scalar_type t_val,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

  template <typename scalar_type>
  static void quadrature_integration_G_q_t_mt(
      int nr_threads, scalar_type beta, scalar_type f_val, scalar_type t_val,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

private:
  template <typename scalar_type>
  static void* quadrature_integration_G_q_w_mt(void* data);

  template <typename scalar_type>
  static void* quadrature_integration_G_q_t_mt(void* data);

  template <typename scalar_type>
  struct quadrature_integration_functions {
    quadrature_integration_functions()
        : I_q_ptr(NULL), H_q_ptr(NULL), S_q_ptr(NULL), G_q_ptr(NULL) {}

    scalar_type beta;

    scalar_type f_val;
    scalar_type t_val;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>* I_q_ptr;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>* H_q_ptr;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>* S_q_ptr;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>* G_q_ptr;
  };
};

template <typename parameters_type, typename q_dmn_t>
template <typename scalar_type>
void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_st(
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> G_inv("G_inv", nu::dmn_size());
  LIN_ALG::GEINV<dca::linalg::CPU>::plan<std::complex<scalar_type>> geinv_obj(G_inv);

  for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++) {
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        G_inv(i, j) = I_q(i, j, q_ind) - H_q(i, j, q_ind) - S_q(i, j, q_ind);

    geinv_obj.execute(G_inv);

    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        G_q(i, j, q_ind) = G_inv(i, j);
  }
}

template <typename parameters_type, typename q_dmn_t>
template <typename scalar_type>
void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_mt(
    int nr_threads, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  G_q = 0.;

  quadrature_integration_functions<scalar_type> quadrature_integration_functions_obj;

  quadrature_integration_functions_obj.I_q_ptr = &I_q;
  quadrature_integration_functions_obj.H_q_ptr = &H_q;
  quadrature_integration_functions_obj.S_q_ptr = &S_q;
  quadrature_integration_functions_obj.G_q_ptr = &G_q;

  dca::concurrency::parallelization<dca::concurrency::POSIX_LIBRARY> parallelization_obj;

  parallelization_obj.execute(nr_threads, quadrature_integration_G_q_w_mt<scalar_type>,
                              (void*)&quadrature_integration_functions_obj);
}

template <typename parameters_type, typename q_dmn_t>
template <typename scalar_type>
void* quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_mt(void* void_ptr) {
  typedef quadrature_integration_functions<scalar_type> quadrature_functions_type;

  dca::concurrency::posix_data* data_ptr = static_cast<dca::concurrency::posix_data*>(void_ptr);
  quadrature_functions_type* functions_ptr = static_cast<quadrature_functions_type*>(data_ptr->args);

  int id = data_ptr->id;
  int nr_threads = data_ptr->nr_threads;

  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q =
      *(functions_ptr->I_q_ptr);
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q =
      *(functions_ptr->H_q_ptr);
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q =
      *(functions_ptr->S_q_ptr);
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q =
      *(functions_ptr->G_q_ptr);

  q_dmn_t q_dmn;
  std::pair<int, int> q_bounds = dca::concurrency::util::getBounds(id, nr_threads, q_dmn);

  {
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> G_inv("G_inv", nu::dmn_size());
    LIN_ALG::GEINV<dca::linalg::CPU>::plan<std::complex<scalar_type>> geinv_obj(G_inv);

    for (int q_ind = q_bounds.first; q_ind < q_bounds.second; q_ind += 1) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_inv(i, j) = I_q(i, j, q_ind) - H_q(i, j, q_ind) - S_q(i, j, q_ind);

      geinv_obj.execute(G_inv);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_q(i, j, q_ind) = G_inv(i, j);
    }
  }

  return 0;
}

template <typename parameters_type, typename q_dmn_t>
template <typename scalar_type>
void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_st(
    scalar_type beta, scalar_type f_val, scalar_type t_val,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  G_q = 0.;

  {
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> H_m("H_m", nu::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> L("e_l", nu::dmn_size());
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> V("V_l", nu::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> G_t("e_l", nu::dmn_size());

    for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          H_m(i, j) = H_q(i, j, q_ind) - I_q(i, j, q_ind);

      if (false)
        LIN_ALG::GEEV<dca::linalg::CPU>::execute('V', 'U', H_m, L, V);
      else
        LIN_ALG::GEEV<dca::linalg::CPU>::execute_on_Greens_function_matrix('V', 'U', H_m, L, V);

      for (int i = 0; i < nu::dmn_size(); i++) {
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
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          for (int l = 0; l < nu::dmn_size(); l++)
            G_q(i, j, q_ind) += V(i, l) * G_t[l] * conj(V(j, l));  // G_t[l]*real(conj(V(l,i))*V(l,j));
    }
  }
}

template <typename parameters_type, typename q_dmn_t>
template <typename scalar_type>
void quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_mt(
    int nr_threads, scalar_type beta, scalar_type f_val, scalar_type t_val,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  G_q = 0.;

  quadrature_integration_functions<scalar_type> quadrature_integration_functions_obj;

  quadrature_integration_functions_obj.beta = beta;

  quadrature_integration_functions_obj.t_val = t_val;
  quadrature_integration_functions_obj.f_val = f_val;

  quadrature_integration_functions_obj.I_q_ptr = &I_q;
  quadrature_integration_functions_obj.H_q_ptr = &H_q;
  quadrature_integration_functions_obj.G_q_ptr = &G_q;

  dca::concurrency::parallelization<dca::concurrency::POSIX_LIBRARY> parallelization_obj;

  parallelization_obj.execute(nr_threads, quadrature_integration_G_q_t_mt<scalar_type>,
                              (void*)&quadrature_integration_functions_obj);
}

template <typename parameters_type, typename q_dmn_t>
template <typename scalar_type>
void* quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_mt(void* void_ptr) {
  typedef quadrature_integration_functions<scalar_type> quadrature_functions_type;

  dca::concurrency::posix_data* data_ptr = static_cast<dca::concurrency::posix_data*>(void_ptr);
  quadrature_functions_type* functions_ptr = static_cast<quadrature_functions_type*>(data_ptr->args);

  int id = data_ptr->id;
  int nr_threads = data_ptr->nr_threads;

  double beta = functions_ptr->beta;

  double t_val = functions_ptr->t_val;
  double f_val = functions_ptr->f_val;

  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q =
      *(functions_ptr->I_q_ptr);
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q =
      *(functions_ptr->H_q_ptr);
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q =
      *(functions_ptr->G_q_ptr);

  q_dmn_t q_dmn;
  std::pair<int, int> q_bounds = dca::concurrency::util::getBounds(id, nr_threads, q_dmn);

  {
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> H_m("H_m", nu::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> L("e_l", nu::dmn_size());
    dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> V("V_l", nu::dmn_size());

    dca::linalg::Vector<scalar_type, dca::linalg::CPU> G_t("e_l", nu::dmn_size());

    for (int q_ind = q_bounds.first; q_ind < q_bounds.second; q_ind += 1) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          H_m(i, j) = H_q(i, j, q_ind) - I_q(i, j, q_ind);

      if (false)
        LIN_ALG::GEEV<dca::linalg::CPU>::execute('V', 'U', H_m, L, V);
      else
        LIN_ALG::GEEV<dca::linalg::CPU>::execute_on_Greens_function_matrix('V', 'U', H_m, L, V);

      for (int i = 0; i < nu::dmn_size(); i++) {
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

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          for (int l = 0; l < nu::dmn_size(); l++)
            G_q(i, j, q_ind) += V(i, l) * G_t[l] * conj(V(j, l));
    }
  }

  return 0;
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_QUADRATURE_INTEGRATION_H
