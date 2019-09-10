// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Staar (taa@zurich.ibm.com)
//
// This class symmetrizes single-particle Greens functions according to cluster symmetries,
// matsubara frequencies and band-index symmetries.

/*
 *  \section tau imaginary-time domain
 *
 *   \f{eqnarray*}{
 *     G(\tau) &=& -G(\tau+\beta)
 *   \f}
 *
 *  \section omega matsubara-frequency domain
 *
 *   \f{eqnarray*}{
 *     G(\varpi) &=& \overline{G(-\varpi)}
 *   \f}
 *
 *  \section r_and_k cluster domain
 *
 *   For each symmetry operation \f$\mathcal{S}\f$ of the cluster-domain, we have
 *
 *   \f{eqnarray*}{
 *     G(\vec{r}) &=& G(\mathcal{S}(\vec{r})) \\
 *     G(\vec{k}) &=& G(\mathcal{S}(\vec{k})) \\
 *   \f}
 */

#ifndef DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_SINGLE_PARTICLE_FUNCTION_HPP
#define DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_SINGLE_PARTICLE_FUNCTION_HPP

#include <cmath>
#include <complex>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <class T>
std::vector<T> operator-(const std::vector<T>& v) {
  auto v_cpy = v;
  for (auto& x : v_cpy)
    x = -x;
  return v_cpy;
}

class symmetrize_single_particle_function {
public:
  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;
  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using w_VERTEX_EXTENDED = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;
  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

protected:
  template <typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
  static void execute(
      func::function<scalartype, func::dmn_variadic<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1>>& f,
      func::function<int, func::dmn_variadic<nu_dmn_t, nu_dmn_t>>& H_symmetry, bool do_diff = false);

  template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
            domains::CLUSTER_SHAPE S>
  static void execute(
      func::function<scalartype,
                     func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S>>>& f,
      bool /*do_diff*/ = false) {
    return executeCluster(f);
  }

  template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
            domains::CLUSTER_SHAPE S>
  static void execute(
      func::function<scalartype,
                     func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>& f,
      bool /*do_diff*/ = false) {
    return executeCluster(f);
  }

  template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(func::function<scalartype, func::dmn_variadic<b, b, f_dmn_0, f_dmn_1>>& f,
                      bool do_diff = false);

  template <typename scalartype, typename f_dmn_0>
  static void execute(func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0>>& f,
                      bool do_diff = false);

  template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0, f_dmn_1>>& f,
                      bool do_diff = false);

  inline static const auto& getInterbandSym() {
    static std::vector<std::vector<int>> symmetries(16);
    static std::once_flag flag;
    std::call_once(flag, [&] {
      symmetries[0] = std::vector<int>{0};
      symmetries[1] = std::vector<int>{1, 3, -4, -12};
      symmetries[2] = std::vector<int>{2, -8};
      symmetries[3] = symmetries[1];

      symmetries[4] = -symmetries[1];
      symmetries[5] = std::vector<int>{};
      symmetries[6] = std::vector<int>{6, 14, -11, -9};
      symmetries[7] = std::vector<int>{};

      symmetries[8] = -symmetries[2];
      symmetries[9] = -symmetries[6];
      symmetries[10] = std::vector<int>{};
      symmetries[11] = -symmetries[6];

      symmetries[12] = -symmetries[1];
      symmetries[13] = std::vector<int>{};
      symmetries[14] = symmetries[6];
      symmetries[15] = std::vector<int>{};
    });

    return symmetries;
  }

public:
  static bool differenceDetected() {
    return difference_detected_;
  }

private:
  template <typename scalartype>
  static void difference(scalartype val, std::string function_name, std::string dmn_name);

  template <typename scalartype>
  static void difference(scalartype val0, scalartype val1, std::string function_name,
                         std::string dmn_name);

  template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void symmetrize_over_electron_spin(
      func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0, f_dmn_1>>& f, bool do_diff);

  template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
            domains::CLUSTER_SHAPE S>
  static void executeCluster(
      func::function<scalartype,
                     func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S>>>& f,
      bool do_diff = false);

  template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
            domains::CLUSTER_SHAPE S>
  static void executeCluster(
      func::function<scalartype,
                     func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>& f,
      bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, t>& f, bool do_diff = false);

  template <typename scalartype, class ClusterDmn>
  static void executeTimeOrFreq(func::function<scalartype, func::dmn_variadic<b, b, ClusterDmn, t>>& f,
                                bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, w>& f, bool do_diff = false);

  template <typename scalartype, typename ClusterDmn>
  static void executeTimeOrFreq(func::function<scalartype, func::dmn_variadic<b, b, ClusterDmn, w>>& f,
                                bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, w_REAL>& f, bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, func::dmn_variadic<b, b, w_REAL>>& f,
                      bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, w_VERTEX>& f, bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, w_VERTEX_EXTENDED>& f, bool do_diff = false);

  template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
            domains::CLUSTER_SHAPE S>
  static void executeCluster(
      func::function<
          scalartype,
          func::dmn_variadic<
              b, b, func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S>>>>& f,
      bool do_diff = false);

  template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
            domains::CLUSTER_SHAPE S>
  static void executeCluster(
      func::function<
          scalartype,
          func::dmn_variadic<
              b, b, func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>>& f,
      bool do_diff = false);

  template <typename ClusterDmn>
  static int oppositeSite(int idx);

  static bool difference_detected_;
};

bool symmetrize_single_particle_function::difference_detected_ = false;

template <typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1>>& f,
    func::function<int, func::dmn_variadic<nu_dmn_t, nu_dmn_t>>& /*H_symmetry*/, bool do_diff) {
  execute(f, do_diff);
}

template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<b, b, f_dmn_0, f_dmn_1>>& f, bool do_diff) {
  {
    func::function<scalartype, f_dmn_0> f0(f.get_name());

    for (int nu_0 = 0; nu_0 < b::dmn_size(); ++nu_0) {
      for (int nu_1 = 0; nu_1 < b::dmn_size(); ++nu_1) {
        for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1) {
          for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
            f0(ind_0) = f(nu_0, nu_1, ind_0, ind_1);

          symmetrize_single_particle_function::execute(f0, do_diff);

          for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
            f(nu_0, nu_1, ind_0, ind_1) = f0(ind_0);
        }
      }
    }
  }

  {
    func::function<scalartype, f_dmn_1> f1(f.get_name());

    for (int nu_0 = 0; nu_0 < b::dmn_size(); ++nu_0) {
      for (int nu_1 = 0; nu_1 < b::dmn_size(); ++nu_1) {
        for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0) {
          for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1)
            f1(ind_1) = f(nu_0, nu_1, ind_0, ind_1);

          symmetrize_single_particle_function::execute(f1, do_diff);

          for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1)
            f(nu_0, nu_1, ind_0, ind_1) = f1(ind_1);
        }
      }
    }
  }
}

template <typename scalartype, typename f_dmn_0>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0>>& f, bool do_diff) {
  func::function<scalartype, func::dmn_variadic<b, b, f_dmn_0>> f0(f.get_name());

  for (int spin_ind = 0; spin_ind < s::dmn_size(); ++spin_ind) {
    for (int b_0 = 0; b_0 < b::dmn_size(); ++b_0)
      for (int b_1 = 0; b_1 < b::dmn_size(); ++b_1)
        for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
          f0(b_0, b_1, ind_0) = f(b_0, spin_ind, b_1, spin_ind, ind_0);

    symmetrize_single_particle_function::execute(f0, do_diff);

    for (int b_0 = 0; b_0 < b::dmn_size(); ++b_0)
      for (int b_1 = 0; b_1 < b::dmn_size(); ++b_1)
        for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
          f(b_0, spin_ind, b_1, spin_ind, ind_0) = f0(b_0, b_1, ind_0);
  }
}

template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0, f_dmn_1>>& f, bool do_diff) {
  symmetrize_over_electron_spin(f, do_diff);

  // Symmetrize over real space or momentum.
  func::function<scalartype, func::dmn_variadic<b, b, f_dmn_0>> f0(f.get_name());

  for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1) {
    for (int spin_ind = 0; spin_ind < s::dmn_size(); ++spin_ind) {
      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_0 = 0; b_0 < b::dmn_size(); ++b_0)
          for (int b_1 = 0; b_1 < b::dmn_size(); ++b_1)

            f0(b_0, b_1, ind_0) = f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1);

      symmetrize_single_particle_function::executeCluster(f0, do_diff);

      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_0 = 0; b_0 < b::dmn_size(); ++b_0)
          for (int b_1 = 0; b_1 < b::dmn_size(); ++b_1)
            f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1) = f0(b_0, b_1, ind_0);
    }
  }

  // Symmetrize over time or frequency.
  func::function<scalartype, func::dmn_variadic<b, b, f_dmn_0, f_dmn_1>> f1(f.get_name());

  for (int spin_ind = 0; spin_ind < s::dmn_size(); ++spin_ind) {
    for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1)
      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_1 = 0; b_1 < b::dmn_size(); ++b_1)
          for (int b_0 = 0; b_0 < b::dmn_size(); ++b_0)
            f1(b_0, b_1, ind_0, ind_1) = f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1);

    symmetrize_single_particle_function::executeTimeOrFreq(f1, do_diff);

    for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1)
      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_1 = 0; b_1 < b::dmn_size(); ++b_1)
          for (int b_0 = 0; b_0 < b::dmn_size(); ++b_0)
            f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1) = f1(b_0, b_1, ind_0, ind_1);
  }
}

template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::symmetrize_over_electron_spin(
    func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0, f_dmn_1>>& f, bool /*do_diff*/) {
  for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ind_1++) {
    for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ind_0++) {
      // spin-symmetry ... --> G_(e_UP, e_DN) == G_(e_DN, e_UP) == 0 !!
      for (int i = 0; i < b::dmn_size(); i++) {
        for (int j = 0; j < b::dmn_size(); j++) {
          f(i, 0, j, 1, ind_0, ind_1) = 0;
          f(i, 1, j, 0, ind_0, ind_1) = 0;

          scalartype tmp = (f(i, 0, j, 0, ind_0, ind_1) + f(i, 1, j, 1, ind_0, ind_1)) / 2.;

          f(i, 0, j, 0, ind_0, ind_1) = tmp;
          f(i, 1, j, 1, ind_0, ind_1) = tmp;
        }
      }
    }
  }
}

template <typename scalartype>
void symmetrize_single_particle_function::difference(scalartype val, std::string function_name,
                                                     std::string dmn_name) {
  if (std::abs(val) > 1.e-6) {
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t"
              << std::abs(val) << "\n\n";

    difference_detected_ = true;
    // throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template <typename scalartype>
void symmetrize_single_particle_function::difference(scalartype val0, scalartype val1,
                                                     std::string function_name, std::string dmn_name) {
  if (abs(val0 - val1) > 1.e-6) {
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t"
              << abs(val0 - val1) << "\n\n";

    difference_detected_ = true;
    // throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template <>
void symmetrize_single_particle_function::difference(float val, std::string function_name,
                                                     std::string dmn_name);

template <>
void symmetrize_single_particle_function::difference(float val0, float val1,
                                                     std::string function_name, std::string dmn_name);

template <>
void symmetrize_single_particle_function::difference(std::complex<float> val,
                                                     std::string function_name, std::string dmn_name);

template <>
void symmetrize_single_particle_function::difference(std::complex<float> val0,
                                                     std::complex<float> val1,
                                                     std::string function_name, std::string dmn_name);

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, t>& f, bool do_diff) {
  int shift = t::dmn_size() / 2;

  double max = 0;
  for (int i = 0; i < t::dmn_size() / 2; i++) {
    max = std::max(max, abs((f(i) + f(i + shift)) / 2.));

    scalartype tmp = (f(i) - f(i + shift)) / 2.;

    f(i) = tmp;
    f(i + shift) = -tmp;
  }

  if (do_diff)
    difference(max, f.get_name(), "tau-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype, typename ClusterDmn>
void symmetrize_single_particle_function::executeTimeOrFreq(
    func::function<scalartype, func::dmn_variadic<b, b, ClusterDmn, t>>& f, bool do_diff) {
  func::function<scalartype, func::dmn_variadic<b, b, ClusterDmn, t>> f_new;
  // Antiperiodicity in time.

  int t_0 = t::dmn_size() / 2;

  for (int t_ind = 0; t_ind < t::dmn_size() / 2; ++t_ind) {
    for (int c_ind = 0; c_ind < ClusterDmn::dmn_size(); ++c_ind) {
      for (int b0 = 0; b0 < b::dmn_size(); ++b0) {
        for (int b1 = 0; b1 < b::dmn_size(); ++b1) {
          scalartype tmp = (f(b0, b1, c_ind, t_ind) - f(b0, b1, c_ind, t_ind + t_0)) / 2.;

          f_new(b0, b1, c_ind, t_ind) = tmp;
          f_new(b0, b1, c_ind, t_ind + t_0) = -tmp;
        }
      }
    }
  }

  double max = 0;
  for (int ind = 0; ind < f.size(); ++ind) {
    max = std::max(max, std::abs(f(ind) - f_new(ind)));
  }

  f = std::move(f_new);

  if (do_diff)
    difference(max, f.get_name(), "t-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, w>& f, bool do_diff) {
  double max = 0;
  for (int i = 0; i < w::dmn_size() / 2; i++) {
    max = std::max(max, abs((f(i) - std::conj(f(w::dmn_size() - i - 1))) / 2.));

    scalartype tmp = (f(i) + std::conj(f(w::dmn_size() - i - 1))) / 2.;

    f(i) = tmp;
    f(w::dmn_size() - 1 - i) = std::conj(tmp);
  }

  if (do_diff)
    difference(max, f.get_name(), "w-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype, typename ClusterDomain>
void symmetrize_single_particle_function::executeTimeOrFreq(
    func::function<scalartype, func::dmn_variadic<b, b, ClusterDomain, w>>& f, bool do_diff) {
  func::function<scalartype, func::dmn_variadic<b, b, ClusterDomain, w>> f_new;

  int w_0 = w::dmn_size() - 1;

  for (int w_ind = 0; w_ind < w::dmn_size() / 2; ++w_ind) {
    for (int c_ind = 0; c_ind < ClusterDomain::dmn_size(); ++c_ind) {
      const int opposite_idx = oppositeSite<ClusterDomain>(c_ind);

      for (int b0 = 0; b0 < b::dmn_size(); ++b0) {
        for (int b1 = 0; b1 < b::dmn_size(); ++b1) {
          scalartype tmp_0 = f(b0, b1, c_ind, w_ind);
          scalartype tmp_1 = f(b1, b0, opposite_idx, w_0 - w_ind);

          scalartype tmp = (tmp_0 + std::conj(tmp_1)) / 2.;

          f_new(b0, b1, c_ind, w_ind) = tmp;
          f_new(b1, b0, opposite_idx, w_0 - w_ind) = std::conj(tmp);
        }
      }
    }
  }

  double max = 0;
  for (int ind = 0; ind < f.size(); ++ind) {
    max = std::max(max, abs(f(ind) - f_new(ind)));
  }

  f = std::move(f_new);

  if (do_diff)
    difference(max, f.get_name(), "w-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, w_REAL>& /*f*/,
                                                  bool /*do_diff*/) {}

template <typename scalartype>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<b, b, w_REAL>>& /*f*/, bool /*do_diff*/) {}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, w_VERTEX>& f,
                                                  bool do_diff) {
  double max = 0;
  for (int i = 0; i < w_VERTEX::dmn_size() / 2; i++) {
    max = std::max(max, abs((f(i) - std::conj(f(w_VERTEX::dmn_size() - i - 1))) / 2.));

    scalartype tmp = (f(i) + std::conj(f(w_VERTEX::dmn_size() - i - 1))) / 2.;

    f(i) = tmp;
    f(w_VERTEX::dmn_size() - i - 1) = std::conj(tmp);
  }

  if (do_diff)
    difference(max, "w_VERTEX-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, w_VERTEX_EXTENDED>& f,
                                                  bool do_diff) {
  double max = 0;
  for (int i = 0; i < w_VERTEX_EXTENDED::dmn_size() / 2; i++) {
    max = std::max(max, abs((f(i) - std::conj(f(w_VERTEX_EXTENDED::dmn_size() - i - 1))) / 2.));

    scalartype tmp = (f(i) + std::conj(f(w_VERTEX_EXTENDED::dmn_size() - i - 1))) / 2.;

    f(i) = tmp;
    f(w_VERTEX_EXTENDED::dmn_size() - i - 1) = std::conj(tmp);
  }

  if (do_diff)
    difference(max, "w_VERTEX_EXTENDED-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
          domains::CLUSTER_SHAPE S>
void symmetrize_single_particle_function::executeCluster(
    func::function<scalartype,
                   func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S>>>& f,
    bool do_diff) {
  typedef domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S> r_cluster_type;
  typedef func::dmn_0<r_cluster_type> r_dmn_t;

  typedef
      typename domains::cluster_symmetry<r_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static func::function<std::pair<int, int>,
                        func::dmn_variadic<func::dmn_variadic<r_dmn_t, b>, sym_super_cell_dmn_t>>&
      r_symmetry_matrix = domains::cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, r_dmn_t> f_new;

  f_new = scalartype(0.);

  for (int S_ind = 0; S_ind < sym_super_cell_dmn_t::dmn_size(); ++S_ind) {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); ++r_ind) {
      int R_new_ind = r_symmetry_matrix(r_ind, 0, S_ind).first;

      f_new(r_ind) += f(R_new_ind);
    }
  }

  if (sym_super_cell_dmn_t::dmn_size() > 0)
    f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
    throw std::logic_error(__FUNCTION__);

  double max = 0;
  for (int ind = 0; ind < f.size(); ++ind) {
    max = std::max(max, std::abs(f(ind) - f_new(ind)));

    f(ind) = f_new(ind);
  }

  if (do_diff)
    difference(max, f.get_name(), "r-cluster-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
          domains::CLUSTER_SHAPE S>
void symmetrize_single_particle_function::executeCluster(
    func::function<
        scalartype,
        func::dmn_variadic<
            b, b, func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S>>>>& f,
    bool do_diff) {
  typedef domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S> r_cluster_type;
  typedef func::dmn_0<r_cluster_type> r_dmn_t;

  typedef
      typename domains::cluster_symmetry<r_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static func::function<std::pair<int, int>,
                        func::dmn_variadic<func::dmn_variadic<r_dmn_t, b>, sym_super_cell_dmn_t>>&
      r_symmetry_matrix = domains::cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, func::dmn_variadic<b, b, r_dmn_t>> f_new;

  f_new = scalartype(0.);
  const double norm_diag = 1. / sym_super_cell_dmn_t::dmn_size();

  for (int r_id = 0; r_id < r_dmn_t::dmn_size(); ++r_id)
    // TODO: remove hardcoding
    for (int b_id = 0; b_id < b::dmn_size(); ++b_id) {
      for (int l = 0; l < sym_super_cell_dmn_t::dmn_size(); ++l) {
        const int r_new = r_symmetry_matrix(r_id, b_id, l).first;
        const int b_new = r_symmetry_matrix(r_id, b_id, l).second;

        f_new(b_id, b_id, r_id) += f(b_new, b_new, r_new) * norm_diag;
      }
      const int b1 = b_id;
      const int b2 = 1 - b_id;

      for (auto r_new : getInterbandSym()[r_id]) {
        if (r_new >= 0)
          f_new(b1, b2, r_id) += f(b1, b2, r_new);
        else
          f_new(b1, b2, r_id) -= f(b1, b2, -r_new);
      }
      const int n_points = getInterbandSym()[r_id].size();
      if (n_points)
        f_new(b1, b2, r_id) /= n_points;
    }

  double max = 0;
  for (int r_id = 0; r_id < r_dmn_t::dmn_size(); ++r_id)
    for (int b2 = 0; b2 < b::dmn_size(); ++b2)
      for (int b1 = 0; b1 < b::dmn_size(); ++b1) {
        max = std::max(max, std::abs(f(b1, b2, r_id) - f_new(b1, b2, r_id)));

        f(b1, b2, r_id) = f_new(b1, b2, r_id);
      }

  if (do_diff)
    difference(max, f.get_name(), "r-cluster-domain of the function : " + f.get_name() + "\n");
}  // namespace phys

template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
          domains::CLUSTER_SHAPE S>
void symmetrize_single_particle_function::executeCluster(
    func::function<scalartype,
                   func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>& f,
    bool do_diff) {
  typedef domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S> k_cluster_type;
  typedef func::dmn_0<k_cluster_type> k_dmn_t;

  typedef
      typename domains::cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static func::function<std::pair<int, int>,
                        func::dmn_variadic<func::dmn_variadic<k_dmn_t, b>, sym_super_cell_dmn_t>>&
      k_symmetry_matrix = domains::cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, k_dmn_t> f_new;

  f_new = scalartype(0.);

  for (int S_ind = 0; S_ind < sym_super_cell_dmn_t::dmn_size(); ++S_ind) {
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
      int K_new_ind = k_symmetry_matrix(k_ind, 0, S_ind).first;

      f_new(k_ind) += f(K_new_ind);
    }
  }

  if (sym_super_cell_dmn_t::dmn_size() > 0)
    f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
    throw std::logic_error(__FUNCTION__);

  double max = 0;
  for (int ind = 0; ind < f.size(); ++ind) {
    max = std::max(max, abs(f(ind) - f_new(ind)));

    f(ind) = f_new(ind);
  }

  if (do_diff)
    difference(max, f.get_name(), "k-cluster-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
          domains::CLUSTER_SHAPE S>
void symmetrize_single_particle_function::executeCluster(
    func::function<
        scalartype,
        func::dmn_variadic<
            b, b, func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>>& f,
    bool do_diff) {
  typedef domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S> k_cluster_type;
  typedef func::dmn_0<k_cluster_type> k_dmn_t;

  typedef
      typename domains::cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static func::function<std::pair<int, int>,
                        func::dmn_variadic<func::dmn_variadic<k_dmn_t, b>, sym_super_cell_dmn_t>>&
      k_symmetry_matrix = domains::cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, func::dmn_variadic<b, b, k_dmn_t>> f_new;
  f_new = 0;

  const double norm_diag = 1. / sym_super_cell_dmn_t::dmn_size();

  for (int k_id = 0; k_id < k_dmn_t::dmn_size(); ++k_id)
    // TODO: remove hardcoding
    for (int b_id = 0; b_id < b::dmn_size(); ++b_id) {
      for (int l = 0; l < sym_super_cell_dmn_t::dmn_size(); ++l) {
        const int k_new = k_symmetry_matrix(k_id, b_id, l).first;
        const int b_new = k_symmetry_matrix(k_id, b_id, l).second;

        f_new(b_id, b_id, k_id) += f(b_new, b_new, k_new) * norm_diag;
      }
      const int b1 = b_id;
      const int b2 = 1 - b_id;

      for (auto k_new : getInterbandSym()[k_id]) {
        if (k_new >= 0)
          f_new(b1, b2, k_id) += f(b1, b2, k_new);
        else
          f_new(b1, b2, k_id) -= f(b1, b2, -k_new);
      }
      const int n_points = getInterbandSym()[k_id].size();
      if (n_points)
        f_new(b1, b2, k_id) /= n_points;
    }

  double max = 0;
  for (int k_id = 0; k_id < k_dmn_t::dmn_size(); ++k_id)
    for (int b2 = 0; b2 < b::dmn_size(); ++b2)
      for (int b1 = 0; b1 < b::dmn_size(); ++b1) {
        max = std::max(max, std::abs(f(b1, b2, k_id) - f_new(b1, b2, k_id)));

        f(b1, b2, k_id) = f_new(b1, b2, k_id);
      }

  if (do_diff)
    difference(max, f.get_name(), "k-clusterdomain of the function : " + f.get_name() + "\n");
}

template <typename ClusterDmn>
int symmetrize_single_particle_function::oppositeSite(const int idx) {
  using Cluster = typename ClusterDmn::parameter_type;
  const int origin = Cluster::origin_index();

  if (Cluster::REPRESENTATION == domains::MOMENTUM_SPACE) {
    return idx;
  }
  else if (Cluster::REPRESENTATION == domains::REAL_SPACE) {
    return Cluster::subtract(idx, origin);
  }
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_SINGLE_PARTICLE_FUNCTION_HPP
