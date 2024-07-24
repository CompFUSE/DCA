// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements ct_aux_walker_tools.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/ct_aux_walker_tools.hpp"

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <complex>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

//
// Definitions of template specialization for CPU
//
template <typename Scalar>
CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::CT_AUX_WALKER_TOOLS(int k_ph)
    : r(k_ph), c(k_ph), d(k_ph) {}

template <typename Scalar>
Scalar CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::consistentScalarDiv(Scalar x, Scalar y) {
  if constexpr (dca::util::IsComplex_t<Scalar>::value) {
    Scalar quot;
    using LocalReal = typename dca::util::RealAlias<Scalar>;
    LocalReal s = (std::fabs(y.real()) + (std::fabs(y.imag())));
    LocalReal oos = 1.0 / s;
    LocalReal ars = x.real() * oos;
    LocalReal ais = x.imag() * oos;
    LocalReal brs = y.real() * oos;
    LocalReal bis = y.imag() * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0 / s;
    quot = {((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos};
    return quot;
  }
  else
    return x / y;
}

template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::compute_Gamma(
    dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma,
    dca::linalg::Matrix<Scalar, dca::linalg::CPU>& N,
    dca::linalg::Matrix<Scalar, dca::linalg::CPU>& G_precomputed,
    dca::linalg::Vector<int, dca::linalg::CPU>& random_vertex_vector,
    dca::linalg::Vector<Scalar, dca::linalg::CPU>& exp_V,
    dca::linalg::Vector<Scalar, dca::linalg::CPU>& exp_delta_V, int /*thread_id*/, int /*stream_id*/) {
  Gamma.resize(random_vertex_vector.size());

  assert(Gamma.nrRows() == Gamma.nrCols());

  int vertex_index = N.nrCols() - G_precomputed.nrCols();

  for (int i = 0; i < Gamma.nrRows(); i++) {
    for (int j = 0; j < Gamma.nrCols(); j++) {
      int configuration_e_spin_index_i = random_vertex_vector[i];
      int configuration_e_spin_index_j = random_vertex_vector[j];

      if (configuration_e_spin_index_j < vertex_index) {
        Real delta = (configuration_e_spin_index_i == configuration_e_spin_index_j) ? 1. : 0.;

        Scalar N_ij = N(configuration_e_spin_index_i, configuration_e_spin_index_j);

        Gamma(i, j) = consistentScalarDiv((N_ij * exp_V[j] - delta), (exp_V[j] - Real(1.)));
      }
      else {
        Gamma(i, j) =
            G_precomputed(configuration_e_spin_index_i, configuration_e_spin_index_j - vertex_index);
      }

      if (i == j) {
        // This significantly improves agreement of complex gamma elements and phase between
        // cpu and gpu.  Both cuda and magma use a similar implementation for complex division.
        // if constexpr (dca::util::IsComplex_t<Scalar>::value) {
        //   auto consistent_complex_div = [](auto gamma_k) {
        //     Scalar quot;
        //     auto y = gamma_k - dca::util::RealAlias<Scalar>(1.0);
        //     using LocalReal = typename dca::util::RealAlias<Scalar>;
        //     LocalReal s = (std::fabs(y.real()) + (std::fabs(y.imag())));
        //     LocalReal oos = 1.0 / s;
        //     LocalReal ars = gamma_k.real() * oos;
        //     LocalReal ais = gamma_k.imag() * oos;
        //     LocalReal brs = y.real() * oos;
        //     LocalReal bis = y.imag() * oos;
        //     s = (brs * brs) + (bis * bis);

        //     oos = 1.0 / s;
        //     quot = {((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos};
        //     return quot;
        //   };

        Scalar gamma_k = exp_delta_V[j];
        auto y = gamma_k - dca::util::RealAlias<Scalar>(1.0);
        Scalar inter_gamma = consistentScalarDiv(gamma_k, y);  //(gamma_k) / (gamma_k - Real(1.));
        Gamma(i, j) -= inter_gamma;                            //(gamma_k) / (gamma_k - Real(1.));
        // }
        // else {
        //   Scalar gamma_k = exp_delta_V[j];
        //   Scalar inter_gamma = (gamma_k) / (gamma_k - Real(1.));
        //   Gamma(i, j) -= inter_gamma;  //(gamma_k) / (gamma_k - Real(1.));
        // }
      }
    }
  }
}

template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(
    dca::linalg::Matrix<Scalar, dca::linalg::CPU>& M, int index) {
  int LD_i = M.leadingDimension();
  int LD_j = M.capacity().second;

  {
    Scalar* M_ptr = &M(0, index);
    for (int i = 0; i < LD_i; ++i)
      M_ptr[i] = 0.;
  }

  {
    Scalar* M_ptr = &M(index, 0);
    for (int j = 0; j < LD_j; ++j)
      M_ptr[j * LD_i] = 0.;
  }

  M(index, index) = 1.;
}

template <typename Scalar>
bool CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::test_max_min(
    int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU, Real max_ref, Real min_ref) {
  Real Gamma_val = std::abs(Gamma_LU(0, 0));
  Real max = Gamma_val;
  Real min = Gamma_val;

  for (int i = 1; i < n + 1; i++) {
    Gamma_val = std::abs(Gamma_LU(i, i));
    max = Gamma_val > max ? Gamma_val : max;
    min = Gamma_val < min ? Gamma_val : min;
  }

  // It's hard to get the min test right and it's not actually a failure state
  // and std::fabs(min_ref - min) < 1.e-12)
  if (std::abs(max_ref - max) < 1.e-12)
    return true;
  else {
    std::cout << __FUNCTION__ << " for Gamma_LU has failed.!\n";
    std::cout.precision(16);
    std::cout << "\n\t n : " << n << "\n";
    for (int i = 1; i < n + 1; i++) {
      Gamma_val = std::abs(Gamma_LU(i, i));
      std::cout << Gamma_val << ", ";
    }
    std::cout << '\n';
    std::cout << std::scientific;
    std::cout << "max"
              << "\t\t"
              << "max_ref"
              << "\t\t"
              << "std::fabs(max_ref - max)" << '\n';
    std::cout << max_ref << "\t" << max << "\t" << std::fabs(max_ref - max) << '\n';
    std::cout << "min"
              << "\t\t"
              << "min_ref"
              << "\t\t"
              << "std::fabs(min_ref - min)" << '\n';
    std::cout << min_ref << "\t" << min << "\t" << std::fabs(min_ref - min) << '\n';
    std::cout << std::endl;

    // Gamma_LU.print();

    return false;
  }
}

template <typename Scalar>
auto CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma(
    int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU, Scalar exp_delta_V, Real& max,
    Real& min) -> Real {
  // solve_Gamma_slow(n, Gamma_LU);
  // solve_Gamma_fast(n, Gamma_LU);
  solve_Gamma_BLAS(n, Gamma_LU);

  Scalar Gamma_LU_n_n = Gamma_LU(n, n);
  Real Gamma_val = std::abs(Gamma_LU_n_n);

  if (n > 0) {
    Real new_max = (Gamma_val > max) ? Gamma_val : max;
    Real new_min = (Gamma_val < min) ? Gamma_val : min;

    if ((new_max / new_min) > 1.e6)
      return 1.e-16;
    else {
      max = new_max;
      min = new_min;
    }
  }
  else {
    max = Gamma_val;
    min = Gamma_val;
  }

  // Here this is fine since we don't reset Gamma_LU to identity as in the blocked solve
#ifndef NDEBUG
  if (!test_max_min(n, Gamma_LU, max, min)) {
    std::cerr << "solve_Gamma test_max_min on Gamma_LU failed!\n";
    throw std::runtime_error("solve_Gamma test_max_min on Gamma_LU failed!");
  }
#endif

  Scalar phani_gamma = exp_delta_V - Real(1.);
  Scalar determinant_ratio = -phani_gamma * Gamma_LU_n_n;

  if (std::imag(determinant_ratio) > std::numeric_limits<Real>::epsilon() * 100) {
    throw(std::logic_error("The determinant is complex."));
  }
  return std::real(determinant_ratio);
}

/*!                 /            |   \          /            |      \ /            |            \
 *                  |            |   |          |            |      | |            |      |
 *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
 *                  |            |   |  ---/    |            |      | |            |      |
 *                  |------------|---|          |------------|------| |------------|------|
 *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
 *
 *   \Gamma_n = L_n*U_n
 *          s = L_n*y
 *          w = x*U_n
 *          d = -x*y+\beta
 */
template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_slow(
    int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU) {
  int LD = Gamma_LU.leadingDimension();

  {
    Scalar* y = Gamma_LU.ptr(0, n);
    Scalar* x = Gamma_LU.ptr(n, 0);

    {
      if (false) {  // serial
        for (int i = 0; i < n; i++)
          for (int j = 0; j < i; j++)
            y[i] -= Gamma_LU(i, j) * y[j];
      }
      else {  // parallell
        for (int j = 0; j < n; j++)
          for (int i = j + 1; i < n; i++)
            y[i] -= Gamma_LU(i, j) * y[j];
      }
    }

    {
      if (true) {  // serial
        for (int j = 0; j < n; j++) {
          for (int i = 0; i < j; i++)
            x[j * LD] -= x[i * LD] * Gamma_LU(i, j);
          x[j * LD] /= Gamma_LU(j, j);
        }
      }
      else {  // parallell
        for (int i = 0; i < n; i++) {
          x[i * LD] /= Gamma_LU(i, i);
          for (int j = i + 1; j < n; j++)
            x[j * LD] -= x[i * LD] * Gamma_LU(i, j);
        }
      }
    }

    for (int i = 0; i < n; i++)
      Gamma_LU(n, n) -= x[i * LD] * y[i];
  }
}

template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_slow(int n, Scalar* Gamma_LU, int LD) {
  {
    Scalar* y = &Gamma_LU[0 + n * LD];  // Gamma_LU.ptr(0, n);
    Scalar* x = &Gamma_LU[n];           // Gamma_LU.ptr(n, 0);

    {
      if (false) {  // serial
        for (int i = 0; i < n; i++)
          for (int j = 0; j < i; j++)
            y[i] -= Gamma_LU[i + j * LD] * y[j];
      }
      else {  // parallell
        for (int j = 0; j < n; j++)
          for (int i = j + 1; i < n; i++)
            y[i] -= Gamma_LU[i + j * LD] * y[j];
      }
    }

    {
      if (true) {  // serial
        for (int j = 0; j < n; j++) {
          for (int i = 0; i < j; i++)
            x[j * LD] -= x[i * LD] * Gamma_LU[i + j * LD];
          x[j * LD] /= Gamma_LU[j + j * LD];
        }
      }
      else {  // parallell
        for (int i = 0; i < n; i++) {
          x[i * LD] /= Gamma_LU[i + i * LD];
          for (int j = i + 1; j < n; j++)
            x[j * LD] -= x[i * LD] * Gamma_LU[i + j * LD];
        }
      }
    }

    for (int i = 0; i < n; i++)
      Gamma_LU[n + n * LD] -= x[i * LD] * y[i];
  }
}

/*!                 /            |   \          /            |      \ /            |            \
 *                  |            |   |          |            |      | |            |      |
 *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
 *                  |            |   |  ---/    |            |      | |            |      |
 *                  |------------|---|          |------------|------| |------------|------|
 *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
 *
 *   \Gamma_n = L_n*U_n
 *          s = L_n*y
 *          w = x*U_n
 *          d = -x*y+\beta
 */
template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_fast(
    int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU) {
  solve_Gamma_fast(n, &Gamma_LU(0, 0), Gamma_LU.leadingDimension());
}

template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_fast(int n, Scalar* A, int LD) {
  {
    Scalar* y = &A[0 + n * LD];

    {  // parallell
      Scalar y_val = 0;

      Scalar* y_ptr = NULL;
      Scalar* G_ptr = NULL;

      for (int j = 0; j < n; j++) {
        y_val = y[j];

        y_ptr = &A[0 + n * LD];
        G_ptr = &A[0 + j * LD];

        for (int i = j + 1; i < n; i++)
          y_ptr[i] -= G_ptr[i] * y_val;
      }
    }

    Scalar* x = &r[0];

    assert(r.size() >= n);
    {
      for (int j = 0; j < n; j++)
        x[j] = A[n + j * LD];
    }

    {  // serial
      Scalar x_val = 0;

      Scalar* x_ptr = NULL;
      Scalar* G_ptr = NULL;

      for (int j = 0; j < n; j++) {
        x_val = x[j];

        x_ptr = x;
        G_ptr = &A[0 + j * LD];

        for (int i = 0; i < j; i++)
          x_val -= x_ptr[i] * G_ptr[i];

        x[j] = consistentScalarDiv(x_val, G_ptr[j]);
      }
    }

    {
      Scalar xy = 0;

      for (int i = 0; i < n; i++)
        xy += x[i] * y[i];

      A[n + n * LD] -= xy;
    }

    {
      for (int j = 0; j < n; j++)
        A[n + j * LD] = x[j];
    }
  }
}

/*!                 /            |   \          /            |      \ /            |            \
 *                  |            |   |          |            |      | |            |      |
 *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
 *                  |            |   |  ---/    |            |      | |            |      |
 *                  |------------|---|          |------------|------| |------------|------|
 *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
 *
 *   \Gamma_n = L_n*U_n
 *          s = L_n*y
 *          w = x*U_n
 *          d = -x*y+\beta
 */
template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_BLAS(
    int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU /*, Scalar exp_delta_V*/) {
  int lda = Gamma_LU.leadingDimension();

  {
    dca::linalg::blas::trsv("L", "N", "U", n, Gamma_LU.ptr(0, 0), lda, Gamma_LU.ptr(0, n), 1);

    dca::linalg::blas::trsv("U", "T", "N", n, Gamma_LU.ptr(0, 0), lda, Gamma_LU.ptr(n, 0), lda);

    {
      Scalar xy = 0;
      for (int i = 0; i < n; i++)
        xy += Gamma_LU(n, i) * Gamma_LU(i, n);

      Gamma_LU(n, n) -= xy;
    }
  }
}

template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_BLAS(
    int n, Scalar* Gamma_LU /*, Scalar exp_delta_V*/, int lda) {
  {
    dca::linalg::blas::trsv("L", "N", "U", n, Gamma_LU, lda, Gamma_LU, 1);
    dca::linalg::blas::trsv("U", "T", "N", n, Gamma_LU, lda, Gamma_LU, lda);
  }
}

template <typename Scalar>
auto CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_blocked(
    int const n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU, Scalar exp_delta_V,
    Real& max, Real& min) -> Scalar {
  // std::cout << "\t(" << min << ", " << max << " ) ";

  solve_Gamma_blocked(n, Gamma_LU);

  Scalar Gamma_LU_n_n = Gamma_LU(n, n);

  Real Gamma_val = std::abs(Gamma_LU_n_n);

  // std::cout << " --> " << Gamma_val << " --> (";

  if (n > 0) {
    Real new_max = (Gamma_val > max) ? Gamma_val : max;
    Real new_min = (Gamma_val < min) ? Gamma_val : min;

    // The Gamma matrix is too ill-conditioned, don't accept this move.
    if ((new_max / new_min) > 1.e6) {
      // For interactions between the same spin type we might do another update with the same Gamma
      // matrix.
      // Since the current diagonal element should not be considered for max/min, we need to already
      // update the Gamma matrix (which will set it to 1).
      // CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_LU, n);
      return Scalar(1.e-16);
    }
    else {
      max = new_max;
      min = new_min;

#ifndef NDEBUG
      // // This has to be here since it will fail almost always when we set Gamma_LU to identity.
      // if (!test_max_min(n, Gamma_LU, max, min)) {
      // 	std::cerr << "solve_Gamma_blocked test_max_min on Gamma_LU failed!\n";
      //   //throw std::runtime_error("solve_Gamma_blocked test_max_min on Gamma_LU failed!");
      // }
#endif
    }
  }
  else {
    max = Gamma_val;
    min = Gamma_val;
  }

  auto phani_gamma = exp_delta_V - Scalar(1.);
  auto determinant_ratio = -phani_gamma * Gamma_LU_n_n;

  return determinant_ratio;
}

/*!                 /            |   \          /            |      \ /            |            \
 *                  |            |   |          |            |      | |            |      |
 *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
 *                  |            |   |  ---/    |            |      | |            |      |
 *                  |------------|---|          |------------|------| |------------|------|
 *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
 *
 *   \Gamma_n = L_n*U_n
 *          s = L_n*y
 *          w = x*U_n
 *          d = -x*y+\beta
 *
 *
 *
 *   /            |      \ /     \       /     \
 *   |            |      | |     |       |     |
 *   |    L_00    |   0  | | y_0 |  ---  | s_0 |
 *   |            |      | |     |  ---  |     |
 *   |------------|------| |-----|       |-----|
 *   \    L_10    | L_11 / \ y_1 /       \ s_1 /
 *
 *   L_00 y_0            = s_0  -->  y_0 = L_00^{-1} * s_0             // L_00^{-1}*s term is
 * precomputed in dtrsm
 *   L_10 y_0 + L_11 y_1 = s_1  -->  y_1 = L_11^{-1} *(s_1 - L_10 y_0)
 *
 *
 *                         /            |      \
 *                         |            |      |
 *                         |    U_00    | U_01 |  ---
 *                         |            |      |  ---
 *                         |------------|------|
 *  [   x_0    |  x_1  ] * \     0      | U_11 /      [   w_0    | w_1 ]
 *
 *   x_0 U_00            = w_0  --> x_0 =  w_0             * U_00^{-1}  //  w_0 * U_00^{-1} term is
 * precomputed in dtrsm
 *   x_0 U_01 + x_1 U_11 = w_1  --> x_1 = (w_1 - x_0 U_01) * U_11^{-1}
 *
 *
 */
template <typename Scalar>
void CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::solve_Gamma_blocked(
    int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU) {
  assert(n > -1 and n < Gamma_LU.size().first);

  int Nk = BLOCK_SIZE;

  int N = Gamma_LU.size().first;
  int LD = Gamma_LU.leadingDimension();

  Scalar* A = &Gamma_LU(0, 0);

  {
    int l = n % Nk;
    int bl = (n - l) / Nk;

    int Ic = bl * Nk;

    assert(n == l + bl * Nk);

    Scalar* A_00 = &A[(bl + 0) * Nk + (bl + 0) * Nk * LD];

    // update diagonal block
    if (Ic > 0 and l > 0) {
      dca::linalg::blas::gemv("N", l, Ic, -1., &A[Ic + 0 * LD], LD, &A[0 + n * LD], 1, 1.,
                              &A[Ic + n * LD], 1);
      dca::linalg::blas::gemv("T", Ic, l, -1., &A[0 + Ic * LD], LD, &A[n + 0 * LD], LD, 1.,
                              &A[n + Ic * LD], LD);
    }

    solve_Gamma_fast(l, A_00, LD);
    // solve_Gamma_slow(l, A_00, LD);

    {
      Scalar xy = 0;
      for (int i = 0; i < Ic; i++)
        xy += A[n + i * LD] * A[i + n * LD];

      A[n + n * LD] -= xy;
    }

    // update non-diagonal block
    if (l > 0 and ((l + 1) % Nk) == 0 and N - (bl + 1) * Nk > 0) {
      {
        Scalar* A_10 = &A[(bl + 1) * Nk + (bl + 0) * Nk * LD];

        for (int l = 0; l < bl; ++l) {
          Scalar* L_il = &A[(bl + 1) * Nk + (l + 0) * Nk * LD];
          Scalar* U_li = &A[(l + 0) * Nk + (bl + 0) * Nk * LD];

          dca::linalg::blas::gemm("N", "N", N - (bl + 1) * Nk, Nk, Nk, -1., L_il, LD, U_li, LD, 1.,
                                  A_10, LD);
        }

        dca::linalg::blas::trsm("R", "U", "N", "N", N - (bl + 1) * Nk, Nk, 1., A_00, LD, A_10, LD);
      }

      {
        Scalar* A_01 = &A[(bl + 0) * Nk + (bl + 1) * Nk * LD];

        for (int l = 0; l < bl; ++l) {
          Scalar* L_il = &A[(bl + 0) * Nk + (l + 0) * Nk * LD];
          Scalar* U_li = &A[(l + 0) * Nk + (bl + 1) * Nk * LD];

          dca::linalg::blas::gemm("N", "N", Nk, N - (bl + 1) * Nk, Nk, -1., L_il, LD, U_li, LD, 1.,
                                  A_01, LD);
        }

        dca::linalg::blas::trsm("L", "L", "N", "U", Nk, N - (bl + 1) * Nk, 1., A_00, LD, A_01, LD);
      }
    }
  }
}

template <typename Scalar>
auto CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::apply_bennett_on_Gamma(
    int k, int n, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& Gamma_LU, Scalar exp_delta_V,
    Real& max, Real& min) -> Real {
  int ld = Gamma_LU.leadingDimension();

  Scalar* r_ptr = &r[0];
  Scalar* c_ptr = &c[0];
  Scalar* d_ptr = &d[0];

  {  // store previous diagonal
    for (int i = 0; i < n; ++i)
      d_ptr[i] = Gamma_LU(i, i);
  }

  {  // remove the column
    for (int i = 0; i < n; ++i)
      r_ptr[i] = 0.;

    for (int i = 0; i < n; ++i)
      c_ptr[i] = 0.;

    r_ptr[k] = -1.;
    c_ptr[k] = -1.;

    for (int i = 0; i < n; ++i) {
      for (int j = 0; i > j and j < k; ++j)
        c_ptr[i] += Gamma_LU(i, j) * Gamma_LU(j, k);

      if (i <= k)
        c_ptr[i] += Gamma_LU(i, k);

      if (i > k)
        c_ptr[i] += Gamma_LU(i, k) * Gamma_LU(k, k);
    }

    dca::linalg::lapack::standardBennet(n, ld, &Gamma_LU(0, 0), c_ptr, r_ptr);
  }

  {  // remove the row
    for (int i = 0; i < n; ++i)
      r_ptr[i] = 0.;

    for (int i = 0; i < n; ++i)
      c_ptr[i] = 0.;

    r_ptr[k] = -1.;
    c_ptr[k] = -1.;

    for (int i = 0; i < n; ++i) {
      for (int j = 0; k > j and j < i; ++j)
        r_ptr[i] += Gamma_LU(k, j) * Gamma_LU(j, i);

      if (k <= i)
        r_ptr[i] += Gamma_LU(k, i);

      if (k > i)
        r_ptr[i] += Gamma_LU(k, i) * Gamma_LU(i, i);
    }

    dca::linalg::lapack::standardBennet(n, ld, &Gamma_LU(0, 0), c_ptr, r_ptr);
  }

  Scalar ratio = 1.;
  for (int i = 0; i < n; ++i)
    ratio *= consistentScalarDiv(Gamma_LU(i, i), d_ptr[i]);

  {
    Real Gamma_val = std::abs(Gamma_LU(0, 0));

    max = Gamma_val;
    min = Gamma_val;

    for (int i = 1; i < n; i++) {
      Gamma_val = std::abs(Gamma_LU(i, i));

      max = Gamma_val > max ? Gamma_val : max;
      min = Gamma_val < min ? Gamma_val : min;
    }

    if ((max / min) > 1.e6)
      return 1.e-16;
  }

  const Scalar phani_gamma = exp_delta_V - Real(1.);
  const Scalar det_ratio = consistentScalarDiv(-ratio, phani_gamma);

  if (std::abs(std::imag(det_ratio)) > std::numeric_limits<Real>::epsilon() * 1000) {
    throw(std::logic_error("The determinant is complex."));
  }

  return std::real(det_ratio);
}

// Template instantiation.
template class CT_AUX_WALKER_TOOLS<dca::linalg::CPU, double>;
template class CT_AUX_WALKER_TOOLS<dca::linalg::CPU, std::complex<double>>;
template class CT_AUX_WALKER_TOOLS<dca::linalg::CPU, float>;
template class CT_AUX_WALKER_TOOLS<dca::linalg::CPU, std::complex<float>>;

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca
