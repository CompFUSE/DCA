// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a class template (i.a. templated on the dimension) for cached NFT.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_CACHED_NFT_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_CACHED_NFT_HPP

#include <cassert>
#include <complex>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/util/ignore.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace detail {
// dca::phys::solver::ctaux::detail::

struct triple {
public:
  static bool less(const triple& l, const triple& r) {
    if (l.first == r.first)
      if (l.second == r.second)
        return (l.third < r.third);
      else
        return (l.second < r.second);
    else
      return (l.first < r.first);
  }

public:
  int first;
  int second;
  int third;
};

}  // detail
// dca::phys::solver::ctaux::

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
class cached_nft {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using r_cluster_type = typename r_dmn_t::parameter_type;

  typedef func::dmn_variadic<b, r_dmn_t> b_r_dmn_t;
  typedef func::dmn_variadic<b, b, r_dmn_t, r_dmn_t, w_vertex_dmn_t, w_vertex_dmn_t>
      b_b_r_dmn_r_dmn_w_vertex_dmn_t_w_vertex_dmn_t;
  typedef func::dmn_variadic<b, b, r_dmn_t, r_dmn_t, w_vertex_pos_dmn_t, w_vertex_dmn_t>
      b_b_r_dmn_r_dmn_w_vertex_pos_dmn_t_w_vertex_dmn_t;

public:
  cached_nft(std::vector<double>& matsubara_frequencies_ref);

  template <typename configuration_t, typename matrix_t>
  double execute(configuration_t& configuration, matrix_t& M,
                 func::function<std::complex<scalar_type>,
                                b_b_r_dmn_r_dmn_w_vertex_dmn_t_w_vertex_dmn_t>& M_r_r_w_w);

  template <typename configuration_t, typename matrix_t>
  double execute(configuration_t& configuration, matrix_t& M,
                 func::function<std::complex<scalar_type>,
                                b_b_r_dmn_r_dmn_w_vertex_pos_dmn_t_w_vertex_dmn_t>& M_r_r_w_w);

private:
  void initialize() {}

  template <typename configuration_t>
  void resize_all(configuration_t& configuration);

  template <typename configuration_t>
  void sort_configuration(configuration_t& configuration);

  template <typename configuration_t>
  void compute_T(configuration_t& configuration);

  template <typename configuration_t, typename matrix_t>
  void compute_M_matrix(configuration_t& configuration, matrix_t& M, int b_i, int r_i, int b_j,
                        int r_j);

  void compute_T_matrices(int b_i, int r_i, int b_j, int r_j);

  double execute_FT();
  double execute_trimmed_FT();

private:
  std::vector<scalar_type> W;

  std::vector<detail::triple> p;

  func::function<int, b_r_dmn_t> start_index;
  func::function<int, b_r_dmn_t> end_index;

  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> T_l_times_M_ij_times_T_r;
  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> M_ij;
  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> T;
  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> T_l;
  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> T_r;
  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> T_l_times_M_ij;
};

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t, w_vertex_pos_dmn_t>::cached_nft(
    std::vector<double>& matsubara_frequencies_ref)
    : W(matsubara_frequencies_ref.size(), 0.),

      p(0),

      start_index("start_index"),
      end_index("end_index") {
  for (size_t l = 0; l < matsubara_frequencies_ref.size(); l++)
    W[l] = matsubara_frequencies_ref[l];

  initialize();
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
template <typename configuration_t, typename matrix_t>
double cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t, w_vertex_pos_dmn_t>::execute(
    configuration_t& configuration, matrix_t& M,
    func::function<std::complex<scalar_type>, b_b_r_dmn_r_dmn_w_vertex_dmn_t_w_vertex_dmn_t>& M_r_r_w_w) {
  double FLOPS = 0;

  int r0_index = r_cluster_type::origin_index();
  int N_r = r_cluster_type::get_size();
  int N_b = b::dmn_size();
  int N_w = w_vertex_dmn_t::dmn_size();

  resize_all(configuration);

  compute_T(configuration);

  for (int b_i = 0; b_i < N_b; b_i++) {
    for (int r_i = 0; r_i < N_r; r_i++) {
      int n_I = end_index(b_i, r_i) - start_index(b_i, r_i);

      for (int b_j = 0; b_j < N_b; b_j++) {
        for (int r_j = 0; r_j < N_r; r_j++) {
          int min_r_j = r_cluster_type::subtract(r_j, r0_index);

          int n_J = end_index(b_j, r_j) - start_index(b_j, r_j);

          if (n_I > 0 && n_J > 0) {
            compute_M_matrix(configuration, M, b_i, r_i, b_j, r_j);

            compute_T_matrices(b_i, r_i, b_j, r_j);

            FLOPS += execute_FT();

            for (int w2 = 0; w2 < N_w; w2++)
              for (int w1 = 0; w1 < N_w; w1++)
                M_r_r_w_w(b_i, b_j, r_i, min_r_j, w1, w2) = T_l_times_M_ij_times_T_r(w1, w2);
          }
          else {
            for (int w2 = 0; w2 < N_w; w2++)
              for (int w1 = 0; w1 < N_w; w1++)
                M_r_r_w_w(b_i, b_j, r_i, min_r_j, w1, w2) = 0;
          }
        }
      }
    }
  }

  return (FLOPS * (1.e-9));
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
template <typename configuration_t, typename matrix_t>
double cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t, w_vertex_pos_dmn_t>::execute(
    configuration_t& configuration, matrix_t& M,
    func::function<std::complex<scalar_type>, b_b_r_dmn_r_dmn_w_vertex_pos_dmn_t_w_vertex_dmn_t>&
        M_r_r_w_w) {
  double FLOPS = 0.;

  int r0_index = r_cluster_type::origin_index();
  int N_r = r_cluster_type::get_size();
  int N_b = b::dmn_size();
  int N_w = w_vertex_dmn_t::dmn_size();
  int N_w_pos = w_vertex_pos_dmn_t::dmn_size();

  resize_all(configuration);

  compute_T(configuration);

  for (int b_i = 0; b_i < N_b; b_i++) {
    for (int r_i = 0; r_i < N_r; r_i++) {
      int n_I = end_index(b_i, r_i) - start_index(b_i, r_i);

      for (int b_j = 0; b_j < N_b; b_j++) {
        for (int r_j = 0; r_j < N_r; r_j++) {
          int min_r_j = r_cluster_type::subtract(r_j, r0_index);

          int n_J = end_index(b_j, r_j) - start_index(b_j, r_j);

          if (n_I > 0 && n_J > 0) {
            compute_M_matrix(configuration, M, b_i, r_i, b_j, r_j);

            compute_T_matrices(b_i, r_i, b_j, r_j);

            FLOPS += execute_trimmed_FT();

            for (int w2 = 0; w2 < N_w; w2++)
              for (int w1 = 0; w1 < N_w_pos; w1++)
                M_r_r_w_w(b_i, b_j, r_i, min_r_j, w1, w2) = T_l_times_M_ij_times_T_r(w1, w2);
          }
          else {
            for (int w2 = 0; w2 < N_w; w2++)
              for (int w1 = 0; w1 < N_w_pos; w1++)
                M_r_r_w_w(b_i, b_j, r_i, min_r_j, w1, w2) = 0;
          }
        }
      }
    }
  }

  return (FLOPS * (1.e-9));
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
template <typename configuration_t>
void cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t, w_vertex_pos_dmn_t>::resize_all(
    configuration_t& configuration) {
  p.resize(configuration.size());

  start_index = 0;
  end_index = configuration.size();

  sort_configuration(configuration);
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
template <typename configuration_t>
void cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t,
                w_vertex_pos_dmn_t>::sort_configuration(configuration_t& configuration) {
  int N_b = b::dmn_size();
  int N_r = r_dmn_t::dmn_size();  // DCA_cluster_type::get_cluster_size();

  int N_v = configuration.size();
  for (int l = 0; l < N_v; l++) {
    p[l].first = configuration[l].get_band();
    p[l].second = configuration[l].get_r_site();
    p[l].third = l;
  }

  sort(p.begin(), p.end(), detail::triple::less);

  for (int b = 0; b < N_b; b++) {
    for (int r = 0; r < N_r; r++) {
      detail::triple t;
      t.first = b;
      t.second = r;

      t.third = 0;

      start_index(b, r) = lower_bound(p.begin(), p.end(), t, detail::triple::less) - p.begin();

      t.third = N_v;

      end_index(b, r) = upper_bound(p.begin(), p.end(), t, detail::triple::less) - p.begin();
    }
  }
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
template <typename configuration_t>
void cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t, w_vertex_pos_dmn_t>::compute_T(
    configuration_t& configuration) {
  int N_v = configuration.size();
  int N_w = W.size();

  T.resizeNoCopy(std::pair<int, int>(N_w, N_v));

  scalar_type x;

  for (int j = 0; j < N_v; j++) {
    for (int i = 0; i < N_w; i++) {
      x = (configuration[j].get_tau() * W[i]);

      T(i, j).real(cos(x));
      T(i, j).imag(sin(x));
    }
  }
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
template <typename configuration_t, typename matrix_t>
void cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t, w_vertex_pos_dmn_t>::compute_M_matrix(
    configuration_t& configuration, matrix_t& M, int b_i, int r_i, int b_j, int r_j) {
  // In release mode 'configuration' is an unused parameter.
  dca::util::ignoreUnused(configuration);

  M_ij.resizeNoCopy(std::pair<int, int>(end_index(b_i, r_i) - start_index(b_i, r_i),
                                        end_index(b_j, r_j) - start_index(b_j, r_j)));

  for (int l_i = start_index(b_i, r_i); l_i < end_index(b_i, r_i); l_i++) {
    assert(p[l_i].first == b_i);
    assert(configuration[p[l_i].third].get_band() == b_i);

    assert(p[l_i].second == r_i);
    assert(configuration[p[l_i].third].get_r_site() == r_i);

    int I = l_i - start_index(b_i, r_i);

    for (int l_j = start_index(b_j, r_j); l_j < end_index(b_j, r_j); l_j++) {
      assert(p[l_j].first == b_j);
      assert(configuration[p[l_j].third].get_band() == b_j);

      assert(p[l_j].second == r_j);
      assert(configuration[p[l_j].third].get_r_site() == r_j);

      int J = l_j - start_index(b_j, r_j);

      M_ij(I, J) = M(p[l_i].third, p[l_j].third);
    }
  }
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
void cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t,
                w_vertex_pos_dmn_t>::compute_T_matrices(int b_i, int r_i, int b_j, int r_j) {
  // T_l matrix
  T_l.resizeNoCopy(
      std::pair<int, int>(w_vertex_dmn_t::dmn_size(), end_index(b_i, r_i) - start_index(b_i, r_i)));

  for (int l_i = start_index(b_i, r_i); l_i < end_index(b_i, r_i); l_i++) {
    int I = l_i - start_index(b_i, r_i);
    memcpy(&T_l(0, I), &T(0, (p[l_i].third)),
           sizeof(std::complex<scalar_type>) * w_vertex_dmn_t::dmn_size());
  }

  // T_r matrix
  T_r.resizeNoCopy(
      std::pair<int, int>(w_vertex_dmn_t::dmn_size(), end_index(b_j, r_j) - start_index(b_j, r_j)));

  for (int l_j = start_index(b_j, r_j); l_j < end_index(b_j, r_j); l_j++) {
    int J = l_j - start_index(b_j, r_j);
    memcpy(&T_r(0, J), &T(0, (p[l_j].third)),
           sizeof(std::complex<scalar_type>) * w_vertex_dmn_t::dmn_size());
  }
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
double cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t, w_vertex_pos_dmn_t>::execute_FT() {
  double FLOPS = 0.;

  assert(T_l.size().first == w_vertex_dmn_t::dmn_size());
  assert(T_l.size().second == M_ij.size().first);

  assert(T_r.size().first == w_vertex_dmn_t::dmn_size());
  assert(T_r.size().second == M_ij.size().second);

  T_l_times_M_ij.resizeNoCopy(std::pair<int, int>(w_vertex_dmn_t::dmn_size(), M_ij.size().second));
  T_l_times_M_ij_times_T_r.resizeNoCopy(
      std::pair<int, int>(w_vertex_dmn_t::dmn_size(), w_vertex_dmn_t::dmn_size()));

  {
    dca::linalg::matrixop::gemm(T_l, M_ij, T_l_times_M_ij);

    {
      int M = T_l.size().first;    // w_VERTEX::dmn_size();//N_w;
      int K = T_l.size().second;   // n_I;
      int N = M_ij.size().second;  // n_J;

      FLOPS += 4. * (M) * (K) * (N);
    }
  }

  {
    dca::linalg::matrixop::gemm('N', 'C', T_l_times_M_ij, T_r, T_l_times_M_ij_times_T_r);

    {
      int M = T_l_times_M_ij.size().first;             // N_w;
      int K = T_l_times_M_ij.size().second;            // n_J;
      int N = T_l_times_M_ij_times_T_r.size().second;  // N_w;

      FLOPS += 4. * (M) * (K) * (N);
    }
  }

  return FLOPS;
}

template <int dimension, class scalar_type, class r_dmn_t, class w_vertex_dmn_t, class w_vertex_pos_dmn_t>
double cached_nft<dimension, scalar_type, r_dmn_t, w_vertex_dmn_t,
                  w_vertex_pos_dmn_t>::execute_trimmed_FT() {
  const static std::complex<scalar_type> ONE(1., 0.);
  const static std::complex<scalar_type> ZERO(0., 0.);

  double FLOPS = 0.;

  assert(w_vertex_pos_dmn_t::dmn_size() == w_vertex_dmn_t::dmn_size() / 2);

  assert(T_l.size().first == w_vertex_dmn_t::dmn_size());
  assert(T_l.size().second == M_ij.size().first);

  assert(T_r.size().first == w_vertex_dmn_t::dmn_size());
  assert(T_r.size().second == M_ij.size().second);

  T_l_times_M_ij.resizeNoCopy(std::pair<int, int>(w_vertex_pos_dmn_t::dmn_size(), M_ij.size().second));
  T_l_times_M_ij_times_T_r.resizeNoCopy(
      std::pair<int, int>(w_vertex_pos_dmn_t::dmn_size(), w_vertex_dmn_t::dmn_size()));

  {
    std::complex<scalar_type>* A = &T_l(w_vertex_pos_dmn_t::dmn_size(), 0);
    std::complex<scalar_type>* B = &M_ij(0, 0);
    std::complex<scalar_type>* C = &T_l_times_M_ij(0, 0);

    int M = T_l.size().first / 2;
    int K = T_l.size().second;
    int N = M_ij.size().second;

    int LDA = T_l.leadingDimension();             // N_w;
    int LDB = M_ij.leadingDimension();            // MAX;
    int LDC = T_l_times_M_ij.leadingDimension();  // N_w;

    dca::linalg::blas::gemm("N", "N", M, N, K, ONE, A, LDA, B, LDB, ZERO, C, LDC);

    FLOPS += 4 * (M) * (K) * (N);
  }

  {
    dca::linalg::matrixop::gemm('N', 'C', T_l_times_M_ij, T_r, T_l_times_M_ij_times_T_r);

    {
      int M = T_l_times_M_ij.size().first;             // N_w;
      int K = T_l_times_M_ij.size().second;            // n_J;
      int N = T_l_times_M_ij_times_T_r.size().second;  // N_w;

      FLOPS += 4. * (M) * (K) * (N);
    }
  }

  return FLOPS;
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_CACHED_NFT_HPP
