// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class symmetrizes Greens functions according to cluster symmetries, matsubara frequencies
// and band-index symmetries.

#ifndef DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_HPP
#define DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_HPP

#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize_single_particle_function.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize_two_particle_function.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
// dca::phys::

class symmetrize : public symmetrize_single_particle_function,
                   public symmetrize_two_particle_function {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = func::dmn_variadic<nu, nu>;

  template <typename scalartype, typename f_dmn_0>
  static void execute(func::function<scalartype, f_dmn_0>& f, bool do_diff = false);

  template <typename scalartype, typename nu_dmn_t, typename f_dmn_0>
  static void execute(func::function<scalartype, func::dmn_variadic<nu_dmn_t, nu_dmn_t, f_dmn_0>>& f,
                      bool do_diff = false);

  template <typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
  static void execute(
      func::function<scalartype, func::dmn_variadic<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1>>& f,
      bool do_diff = false);

  template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0, f_dmn_1>>& f,
                      func::function<int, nu_nu>& H_symmetry, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                    func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                    func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>>>& f,
      func::function<int, nu_nu>& H_symmetry, std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
      func::function<int, nu_nu>& H_symmetry, std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype,
                     func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                        func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>, k_dmn_t>>& f,
      bool do_diff = false);
};

template <typename scalartype, typename f_dmn_0>
void symmetrize::execute(func::function<scalartype, f_dmn_0>& f, bool do_diff) {
  symmetrize_single_particle_function::execute(f, do_diff);
}

template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize::execute(func::function<scalartype, func::dmn_variadic<nu, nu, f_dmn_0, f_dmn_1>>& f,
                         func::function<int, nu_nu>& H_symmetry, bool do_diff) {
  symmetrize_single_particle_function::execute(f, H_symmetry, do_diff);
}

template <typename scalartype, typename nu_dmn_t, typename f_dmn_0>
void symmetrize::execute(func::function<scalartype, func::dmn_variadic<nu_dmn_t, nu_dmn_t, f_dmn_0>>& f,
                         bool do_diff) {
  symmetrize_single_particle_function::execute(f, do_diff);
}

template <typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
void symmetrize::execute(
    func::function<scalartype, func::dmn_variadic<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1>>& f,
    bool do_diff) {
  symmetrize_single_particle_function::execute(f, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    func::function<scalartype, func::dmn_variadic<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
    func::function<int, nu_nu>& /*H_symmetry*/, std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>>>& f,
    func::function<int, nu_nu>& /*H_symmetry*/, std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    func::function<scalartype, func::dmn_variadic<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
    std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>>>& f,
    std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>, k_dmn_t>>& f,
    bool do_diff) {
  symmetrize_two_particle_function::execute(f, do_diff);
}

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_HPP
