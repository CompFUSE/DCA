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

class symmetrize_single_particle_function {
public:
  using TDmn = func::dmn_0<domains::time_domain>;
  using WDmn = func::dmn_0<domains::frequency_domain>;
  using WVertexDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using WVertexExtDmn = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;
  using WRealDmn = func::dmn_0<domains::frequency_domain_real_axis>;

  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;  // orbital-spin index

protected:
  template <class Lattice, typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
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
  static void execute(func::function<scalartype, func::dmn_variadic<BDmn, BDmn, f_dmn_0, f_dmn_1>>& f,
                      bool do_diff = false);

  template <typename scalartype, typename f_dmn_0>
  static void execute(func::function<scalartype, func::dmn_variadic<NuDmn, NuDmn, f_dmn_0>>& f,
                      bool do_diff = false);

  template <class Lattice, typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(func::function<scalartype, func::dmn_variadic<NuDmn, NuDmn, f_dmn_0, f_dmn_1>>& f,
                      bool do_diff = false);

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
      func::function<scalartype, func::dmn_variadic<NuDmn, NuDmn, f_dmn_0, f_dmn_1>>& f,
      bool do_diff);

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
  static void execute(func::function<scalartype, TDmn>& f, bool do_diff = false);

  template <typename scalartype, class ClusterDmn>
  static void executeTimeOrFreq(
      func::function<scalartype, func::dmn_variadic<BDmn, BDmn, ClusterDmn, TDmn>>& f,
      bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, WDmn>& f, bool do_diff = false);

  template <typename scalartype, typename ClusterDmn>
  static void executeTimeOrFreq(
      func::function<scalartype, func::dmn_variadic<BDmn, BDmn, ClusterDmn, WDmn>>& f,
      bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, WRealDmn>& f, bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, func::dmn_variadic<BDmn, BDmn, WRealDmn>>& f,
                      bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, WVertexDmn>& f, bool do_diff = false);

  template <typename scalartype>
  static void execute(func::function<scalartype, WVertexExtDmn>& f, bool do_diff = false);

  template <class Lattice, typename scalartype, typename scalar_type, int D,
            domains::CLUSTER_NAMES N, domains::CLUSTER_SHAPE S>
  static void executeCluster(
      func::function<scalartype, func::dmn_variadic<BDmn, BDmn,
                                                    func::dmn_0<domains::cluster_domain<
                                                        scalar_type, D, N, domains::REAL_SPACE, S>>>>& f,
      bool do_diff = false);

  template <class Lattice, typename scalartype, typename scalar_type, int D,
            domains::CLUSTER_NAMES N, domains::CLUSTER_SHAPE S>
  static void executeCluster(
      func::function<
          scalartype,
          func::dmn_variadic<
              BDmn, BDmn,
              func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>>& f,
      bool do_diff = false);

  template <typename ClusterDmn>
  static int oppositeSite(int idx);

  static bool difference_detected_;
};

bool symmetrize_single_particle_function::difference_detected_ = false;

template <class Lattice, typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1>>& f,
    func::function<int, func::dmn_variadic<nu_dmn_t, nu_dmn_t>>& /*H_symmetry*/, bool do_diff) {
  execute<Lattice>(f, do_diff);
}

template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<BDmn, BDmn, f_dmn_0, f_dmn_1>>& f, bool do_diff) {
  {
    func::function<scalartype, f_dmn_0> f0(f.get_name());

    for (int nu_0 = 0; nu_0 < BDmn::dmn_size(); ++nu_0) {
      for (int nu_1 = 0; nu_1 < BDmn::dmn_size(); ++nu_1) {
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

    for (int nu_0 = 0; nu_0 < BDmn::dmn_size(); ++nu_0) {
      for (int nu_1 = 0; nu_1 < BDmn::dmn_size(); ++nu_1) {
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
    func::function<scalartype, func::dmn_variadic<NuDmn, NuDmn, f_dmn_0>>& f, bool do_diff) {
  func::function<scalartype, func::dmn_variadic<BDmn, BDmn, f_dmn_0>> f0(f.get_name());

  for (int spin_ind = 0; spin_ind < SDmn::dmn_size(); ++spin_ind) {
    for (int b_0 = 0; b_0 < BDmn::dmn_size(); ++b_0)
      for (int b_1 = 0; b_1 < BDmn::dmn_size(); ++b_1)
        for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
          f0(b_0, b_1, ind_0) = f(b_0, spin_ind, b_1, spin_ind, ind_0);

    symmetrize_single_particle_function::execute(f0, do_diff);

    for (int b_0 = 0; b_0 < BDmn::dmn_size(); ++b_0)
      for (int b_1 = 0; b_1 < BDmn::dmn_size(); ++b_1)
        for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
          f(b_0, spin_ind, b_1, spin_ind, ind_0) = f0(b_0, b_1, ind_0);
  }
}

template <class Lattice, typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<NuDmn, NuDmn, f_dmn_0, f_dmn_1>>& f, bool do_diff) {
  symmetrize_over_electron_spin(f, do_diff);

  // Symmetrize over real space or momentum.
  func::function<scalartype, func::dmn_variadic<BDmn, BDmn, f_dmn_0>> f0(f.get_name());

  for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1) {
    for (int spin_ind = 0; spin_ind < SDmn::dmn_size(); ++spin_ind) {
      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_0 = 0; b_0 < BDmn::dmn_size(); ++b_0)
          for (int b_1 = 0; b_1 < BDmn::dmn_size(); ++b_1)

            f0(b_0, b_1, ind_0) = f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1);

      symmetrize_single_particle_function::executeCluster<Lattice>(f0, do_diff);

      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_0 = 0; b_0 < BDmn::dmn_size(); ++b_0)
          for (int b_1 = 0; b_1 < BDmn::dmn_size(); ++b_1)
            f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1) = f0(b_0, b_1, ind_0);
    }
  }

  // Symmetrize over time or frequency.
  func::function<scalartype, func::dmn_variadic<BDmn, BDmn, f_dmn_0, f_dmn_1>> f1(f.get_name());

  for (int spin_ind = 0; spin_ind < SDmn::dmn_size(); ++spin_ind) {
    for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1)
      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_1 = 0; b_1 < BDmn::dmn_size(); ++b_1)
          for (int b_0 = 0; b_0 < BDmn::dmn_size(); ++b_0)
            f1(b_0, b_1, ind_0, ind_1) = f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1);

    symmetrize_single_particle_function::executeTimeOrFreq(f1, do_diff);

    for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ++ind_1)
      for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ++ind_0)
        for (int b_1 = 0; b_1 < BDmn::dmn_size(); ++b_1)
          for (int b_0 = 0; b_0 < BDmn::dmn_size(); ++b_0)
            f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1) = f1(b_0, b_1, ind_0, ind_1);
  }
}

template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::symmetrize_over_electron_spin(
    func::function<scalartype, func::dmn_variadic<NuDmn, NuDmn, f_dmn_0, f_dmn_1>>& f,
    bool /*do_diff*/) {
  for (int ind_1 = 0; ind_1 < f_dmn_1::dmn_size(); ind_1++) {
    for (int ind_0 = 0; ind_0 < f_dmn_0::dmn_size(); ind_0++) {
      // spin-symmetry ... --> G_(e_UP, e_DN) == G_(e_DN, e_UP) == 0 !!
      for (int i = 0; i < BDmn::dmn_size(); i++) {
        for (int j = 0; j < BDmn::dmn_size(); j++) {
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
void symmetrize_single_particle_function::execute(func::function<scalartype, TDmn>& f, bool do_diff) {
  int shift = TDmn::dmn_size() / 2;

  double max = 0;
  for (int i = 0; i < TDmn::dmn_size() / 2; i++) {
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
    func::function<scalartype, func::dmn_variadic<BDmn, BDmn, ClusterDmn, TDmn>>& f, bool do_diff) {
  func::function<scalartype, func::dmn_variadic<BDmn, BDmn, ClusterDmn, TDmn>> f_new;
  // Antiperiodicity in time.

  int t_0 = TDmn::dmn_size() / 2;

  for (int t_ind = 0; t_ind < TDmn::dmn_size() / 2; ++t_ind) {
    for (int c_ind = 0; c_ind < ClusterDmn::dmn_size(); ++c_ind) {
      for (int b0 = 0; b0 < BDmn::dmn_size(); ++b0) {
        for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1) {
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
    difference(max, f.get_name(), "TDmn-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, WDmn>& f, bool do_diff) {
  double max = 0;
  for (int i = 0; i < WDmn::dmn_size() / 2; i++) {
    max = std::max(max, abs((f(i) - std::conj(f(WDmn::dmn_size() - i - 1))) / 2.));

    scalartype tmp = (f(i) + std::conj(f(WDmn::dmn_size() - i - 1))) / 2.;

    f(i) = tmp;
    f(WDmn::dmn_size() - 1 - i) = std::conj(tmp);
  }

  if (do_diff)
    difference(max, f.get_name(), "WDmn-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype, typename ClusterDomain>
void symmetrize_single_particle_function::executeTimeOrFreq(
    func::function<scalartype, func::dmn_variadic<BDmn, BDmn, ClusterDomain, WDmn>>& f, bool do_diff) {
  func::function<scalartype, func::dmn_variadic<BDmn, BDmn, ClusterDomain, WDmn>> f_new;

  int w_0 = WDmn::dmn_size() - 1;
  constexpr auto representation = ClusterDomain::parameter_type::REPRESENTATION;

  for (int w_ind = 0; w_ind < WDmn::dmn_size() / 2; ++w_ind) {
    for (int c_ind = 0; c_ind < ClusterDomain::dmn_size(); ++c_ind) {
      const int new_c_idx =
          representation == domains::REAL_SPACE ? oppositeSite<ClusterDomain>(c_ind) : c_ind;

      for (int b0 = 0; b0 < BDmn::dmn_size(); ++b0) {
        for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1) {
          scalartype tmp_0 = f(b0, b1, c_ind, w_ind);
          scalartype tmp_1 = f(b1, b0, new_c_idx, w_0 - w_ind);

          scalartype tmp = (tmp_0 + std::conj(tmp_1)) / 2.;

          f_new(b0, b1, c_ind, w_ind) = tmp;
          f_new(b1, b0, new_c_idx, w_0 - w_ind) = std::conj(tmp);
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
    difference(max, f.get_name(), "WDmn-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, WRealDmn>& /*f*/,
                                                  bool /*do_diff*/) {}

template <typename scalartype>
void symmetrize_single_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<BDmn, BDmn, WRealDmn>>& /*f*/, bool /*do_diff*/) {}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, WVertexDmn>& f,
                                                  bool do_diff) {
  double max = 0;
  for (int i = 0; i < WVertexDmn::dmn_size() / 2; i++) {
    max = std::max(max, abs((f(i) - std::conj(f(WVertexDmn::dmn_size() - i - 1))) / 2.));

    scalartype tmp = (f(i) + std::conj(f(WVertexDmn::dmn_size() - i - 1))) / 2.;

    f(i) = tmp;
    f(WVertexDmn::dmn_size() - i - 1) = std::conj(tmp);
  }

  if (do_diff)
    difference(max, "WVertexDmn-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype>
void symmetrize_single_particle_function::execute(func::function<scalartype, WVertexExtDmn>& f,
                                                  bool do_diff) {
  double max = 0;
  for (int i = 0; i < WVertexExtDmn::dmn_size() / 2; i++) {
    max = std::max(max, abs((f(i) - std::conj(f(WVertexExtDmn::dmn_size() - i - 1))) / 2.));

    scalartype tmp = (f(i) + std::conj(f(WVertexExtDmn::dmn_size() - i - 1))) / 2.;

    f(i) = tmp;
    f(WVertexExtDmn::dmn_size() - i - 1) = std::conj(tmp);
  }

  if (do_diff)
    difference(max, "WVertexExtDmn-domain of the function : " + f.get_name() + "\n");
}

template <typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
          domains::CLUSTER_SHAPE S>
void symmetrize_single_particle_function::executeCluster(
    func::function<scalartype,
                   func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S>>>& f,
    bool do_diff) {
  using RCluster = domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S>;
  using RDmn = func::dmn_0<RCluster>;

  typedef typename domains::cluster_symmetry<RCluster>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static func::function<std::pair<int, int>,
                        func::dmn_variadic<func::dmn_variadic<RDmn, BDmn>, sym_super_cell_dmn_t>>&
      r_symmetry_matrix = domains::cluster_symmetry<RCluster>::get_symmetry_matrix();

  static func::function<scalartype, RDmn> f_new;

  f_new = scalartype(0.);

  for (int S_ind = 0; S_ind < sym_super_cell_dmn_t::dmn_size(); ++S_ind) {
    for (int r_ind = 0; r_ind < RDmn::dmn_size(); ++r_ind) {
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

template <class Lattice, typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
          domains::CLUSTER_SHAPE S>
void symmetrize_single_particle_function::executeCluster(
    func::function<scalartype, func::dmn_variadic<BDmn, BDmn,
                                                  func::dmn_0<domains::cluster_domain<
                                                      scalar_type, D, N, domains::REAL_SPACE, S>>>>& f,
    bool do_diff) {
  typedef domains::cluster_domain<scalar_type, D, N, domains::REAL_SPACE, S> r_cluster_type;
  typedef func::dmn_0<r_cluster_type> RDmn;

  using SymDmn = typename domains::cluster_symmetry<r_cluster_type>::sym_super_cell_dmn_t;

  const auto& r_symmetry_matrix = domains::cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, func::dmn_variadic<BDmn, BDmn, RDmn>> f_new;

  f_new = scalartype(0.);

  for (int r_ind = 0; r_ind < RDmn::dmn_size(); ++r_ind) {
    for (int b0 = 0; b0 < BDmn::dmn_size(); ++b0) {
      for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1) {
        double norm = 0.;
        for (int s_ind = 0; s_ind < SymDmn::dmn_size(); ++s_ind) {
          int R_new_ind = r_symmetry_matrix(r_ind, 0, s_ind).first;

          int b0_new = r_symmetry_matrix(r_ind, b0, s_ind).second;
          int b1_new = r_symmetry_matrix(0, b1, s_ind).second;

          const double sign = Lattice::transformationSignOfR(b0, b1, s_ind);
          norm += std::abs(sign);

          f_new(b0, b1, r_ind) += sign * f(b0_new, b1_new, R_new_ind);
        }
        assert(std::abs(norm) > 0);
        f_new(b0, b1, r_ind) /= norm;
      }
    }
  }

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
    func::function<scalartype,
                   func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>& f,
    bool do_diff) {
  typedef domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S> k_cluster_type;
  typedef func::dmn_0<k_cluster_type> k_dmn_t;

  using SymDmn = typename domains::cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t;

  const auto& k_symmetry_matrix = domains::cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, k_dmn_t> f_new;

  f_new = scalartype(0.);

  for (int S_ind = 0; S_ind < SymDmn::dmn_size(); ++S_ind) {
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
      int K_new_ind = k_symmetry_matrix(k_ind, 0, S_ind).first;

      f_new(k_ind) += f(K_new_ind);
    }
  }

  if (SymDmn::dmn_size() > 0)
    f_new /= double(SymDmn::dmn_size());
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

template <class Lattice, typename scalartype, typename scalar_type, int D, domains::CLUSTER_NAMES N,
          domains::CLUSTER_SHAPE S>
void symmetrize_single_particle_function::executeCluster(
    func::function<
        scalartype,
        func::dmn_variadic<
            BDmn, BDmn,
            func::dmn_0<domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S>>>>& f,
    bool do_diff) {
  typedef domains::cluster_domain<scalar_type, D, N, domains::MOMENTUM_SPACE, S> k_cluster_type;
  typedef func::dmn_0<k_cluster_type> k_dmn_t;

  typedef
      typename domains::cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  const auto& k_symmetry_matrix = domains::cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, func::dmn_variadic<BDmn, BDmn, k_dmn_t>> f_new;

  f_new = scalartype(0.);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
    for (int b0 = 0; b0 < BDmn::dmn_size(); ++b0) {
      for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1) {
        double norm = 0.;
        for (int s_ind = 0; s_ind < sym_super_cell_dmn_t::dmn_size(); ++s_ind) {
          int k_new = k_symmetry_matrix(k_ind, b0, s_ind).first;  // FIXME: b0 -> b1

          int b0_new = k_symmetry_matrix(k_ind, b0, s_ind).second;
          int b1_new = k_symmetry_matrix(k_ind, b1, s_ind).second;

          const double sign = Lattice::transformationSignOfK(b0, b1, s_ind);
          norm += std::abs(sign);

          f_new(b0, b1, k_ind) += sign * f(b0_new, b1_new, k_new);
        }
        assert(std::abs(norm) > 0);
        f_new(b0, b1, k_ind) /= norm;
      }
    }
  }

  double max = 0;
  for (int ind = 0; ind < f.size(); ++ind) {
    max = std::max(max, std::abs(f(ind) - f_new(ind)));

    f(ind) = f_new(ind);
  }

  if (do_diff)
    difference(max, f.get_name(), "k-clusterdomain of the function : " + f.get_name() + "\n");
}

template <typename ClusterDmn>
int symmetrize_single_particle_function::oppositeSite(const int idx) {
  using Cluster = typename ClusterDmn::parameter_type;
  const int origin = Cluster::origin_index();
  return Cluster::subtract(idx, origin);
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_SINGLE_PARTICLE_FUNCTION_HPP
