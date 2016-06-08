// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This class symmetrizes Greens functions according to cluster symmetries, matsubara frequencies
// and band-index symmetries.

#ifndef PHYS_LIBRARY_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_H
#define PHYS_LIBRARY_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_H

#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize_single_particle_function.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize_two_particle_function.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

class symmetrize : public symmetrize_single_particle_function,
                   public symmetrize_two_particle_function {
public:
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = dmn_variadic<nu, nu>;

public:
  template <typename scalartype, typename f_dmn_0>
  static void execute(FUNC_LIB::function<scalartype, f_dmn_0>& f, bool do_diff = false);

  template <typename scalartype, typename nu_dmn_t, typename f_dmn_0>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<nu_dmn_t, nu_dmn_t, f_dmn_0>>& f,
                      bool do_diff = false);

  template <typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
  static void execute(FUNC_LIB::function<scalartype, dmn_4<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1>>& f,
                      bool do_diff = false);

  template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(FUNC_LIB::function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1>>& f,
                      FUNC_LIB::function<int, nu_nu>& H_symmetry, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      FUNC_LIB::function<scalartype,
                         dmn_2<dmn_4<b, b, k_dmn_t, w_dmn_t>, dmn_4<b, b, k_dmn_t, w_dmn_t>>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      FUNC_LIB::function<scalartype,
                         dmn_2<dmn_4<b, b, k_dmn_t, w_dmn_t>, dmn_4<b, b, k_dmn_t, w_dmn_t>>>& f,
      FUNC_LIB::function<int, nu_nu>& H_symmetry, std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      FUNC_LIB::function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      FUNC_LIB::function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
      FUNC_LIB::function<int, nu_nu>& H_symmetry, std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<dmn_4<b, b, k_dmn_t, w_dmn_t>,
                                                           dmn_4<b, b, k_dmn_t, w_dmn_t>, k_dmn_t>>& f,
                      bool do_diff = false);
};

template <typename scalartype, typename f_dmn_0>
void symmetrize::execute(FUNC_LIB::function<scalartype, f_dmn_0>& f, bool do_diff) {
  symmetrize_single_particle_function::execute(f, do_diff);
}

template <typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize::execute(FUNC_LIB::function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1>>& f,
                         FUNC_LIB::function<int, nu_nu>& H_symmetry, bool do_diff) {
  symmetrize_single_particle_function::execute(f, H_symmetry, do_diff);
}

template <typename scalartype, typename nu_dmn_t, typename f_dmn_0>
void symmetrize::execute(FUNC_LIB::function<scalartype, dmn_3<nu_dmn_t, nu_dmn_t, f_dmn_0>>& f,
                         bool do_diff) {
  symmetrize_single_particle_function::execute(f, do_diff);
}

template <typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
void symmetrize::execute(FUNC_LIB::function<scalartype, dmn_4<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1>>& f,
                         bool do_diff) {
  symmetrize_single_particle_function::execute(f, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    FUNC_LIB::function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
    FUNC_LIB::function<int, nu_nu>& /*H_symmetry*/, std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    FUNC_LIB::function<scalartype, dmn_2<dmn_4<b, b, k_dmn_t, w_dmn_t>, dmn_4<b, b, k_dmn_t, w_dmn_t>>>& f,
    FUNC_LIB::function<int, nu_nu>& /*H_symmetry*/, std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    FUNC_LIB::function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
    std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    FUNC_LIB::function<scalartype, dmn_2<dmn_4<b, b, k_dmn_t, w_dmn_t>, dmn_4<b, b, k_dmn_t, w_dmn_t>>>& f,
    std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(
    FUNC_LIB::function<scalartype,
                       dmn_3<dmn_4<b, b, k_dmn_t, w_dmn_t>, dmn_4<b, b, k_dmn_t, w_dmn_t>, k_dmn_t>>& f,
    bool do_diff) {
  symmetrize_two_particle_function::execute(f, do_diff);
}

#endif  // PHYS_LIBRARY_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_H
