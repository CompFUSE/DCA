// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// Description

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_UPDATE_CHEMICAL_POTENTIAL_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_UPDATE_CHEMICAL_POTENTIAL_HPP

#include <complex>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/util/print_time.hpp"

#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"

namespace dca {
namespace phys {
// dca::phys::

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
class update_chemical_potential {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using t = func::dmn_0<time_domain>;
  using w = func::dmn_0<frequency_domain>;

  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using DCA_r_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                            CLUSTER, REAL_SPACE, BRILLOUIN_ZONE>;
  using r_DCA = func::dmn_0<DCA_r_cluster_type>;
  using DCA_k_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                            CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;

public:
  update_chemical_potential(parameters_type& parameters_ref, MOMS_type& MOMS_ref,
                            coarsegraining_type& coarsegraining_ref);

  void execute();

  double compute_density();

private:
  double get_new_chemical_potential(double d_0,

                                    double mu_lb, double mu_ub,

                                    double n_lb, double n_ub);

  void compute_density_correction(func::function<double, nu>& result);

  void compute_density_coefficients(
      func::function<double, func::dmn_variadic<nu, k_DCA>>& A,
      func::function<double, func::dmn_variadic<nu, k_DCA>>& B,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& G);

  void search_bounds(double dens);

  void print_bounds();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  MOMS_type& MOMS;
  coarsegraining_type& coarsegraining;

  std::pair<double, double> lower_bound;
  std::pair<double, double> upper_bound;
};

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::update_chemical_potential(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref, coarsegraining_type& coarsegraining_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      MOMS(MOMS_ref),
      coarsegraining(coarsegraining_ref) {}

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::execute() {
  double dens = compute_density();

  if (concurrency.id() == 0) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::cout << "\n\t\t initial chemical potential : " << parameters.get_chemical_potential()
              << " (" << dens << ")\n\n";
  }

  if (std::abs(dens - parameters.get_density()) < 1.e-3)
    return;

  search_bounds(dens);

  while (true) {
    double d_0 = parameters.get_density();

    double mu_lb = lower_bound.first;
    double mu_ub = upper_bound.first;

    double n_lb = lower_bound.second;
    double n_ub = upper_bound.second;

    parameters.get_chemical_potential() = get_new_chemical_potential(d_0, mu_lb, mu_ub, n_lb, n_ub);

    dens = compute_density();

    if (std::abs(dens - parameters.get_density()) < 1.e-3) {
      if (concurrency.id() == 0) {
        std::cout.precision(6);
        std::cout << std::scientific;

        std::cout << "\n\t\t final chemical potential : " << parameters.get_chemical_potential()
                  << " (" << dens << ")\n";
      }

      break;
    }
    else {
      if (dens < parameters.get_density()) {
        lower_bound.first = parameters.get_chemical_potential();
        lower_bound.second = dens;
      }
      else {
        upper_bound.first = parameters.get_chemical_potential();
        upper_bound.second = dens;
      }
    }

    print_bounds();
  }
}

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
double update_chemical_potential<parameters_type, MOMS_type,
                                 coarsegraining_type>::get_new_chemical_potential(double d_0,

                                                                                  double mu_lb,
                                                                                  double mu_ub,

                                                                                  double n_lb,
                                                                                  double n_ub) {
  double new_mu = (mu_ub - mu_lb) / (n_ub - n_lb) * (d_0 - n_lb) + mu_lb;
  return new_mu;
}

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::search_bounds(
    double dens) {
  double factor = 2;
  double delta = 0.1;

  if (dens < parameters.get_density()) {
    lower_bound.first = parameters.get_chemical_potential();
    lower_bound.second = dens;

    while (true) {
      parameters.get_chemical_potential() += delta;

      dens = compute_density();

      upper_bound.first = parameters.get_chemical_potential();
      upper_bound.second = dens;

      print_bounds();

      if (parameters.get_density() < dens)
        break;
      else
        lower_bound = upper_bound;

      delta *= factor;
    }
  }
  else {
    upper_bound.first = parameters.get_chemical_potential();
    upper_bound.second = dens;

    while (true) {
      parameters.get_chemical_potential() -= delta;

      dens = compute_density();

      lower_bound.first = parameters.get_chemical_potential();
      lower_bound.second = dens;

      print_bounds();

      if (dens < parameters.get_density())
        break;
      else
        upper_bound = lower_bound;

      delta *= factor;
    }
  }
}

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
double update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::compute_density() {
  if (parameters.do_DCA_plus())
    coarsegraining.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma_lattice, MOMS.G_k_w);
  else
    coarsegraining.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma, MOMS.G_k_w);

  MOMS.G_k_w -= MOMS.G0_k_w;

  math::transform::FunctionTransform<w, t>::execute(MOMS.G_k_w, MOMS.G_k_t);

  MOMS.G_k_t += MOMS.G0_k_t;

  math::transform::FunctionTransform<k_DCA, r_DCA>::execute(MOMS.G_k_t, MOMS.G_r_t);

  MOMS.G_k_w += MOMS.G0_k_w;

  func::function<double, nu> result;
  result = 0.0;
  compute_density_correction(result);
  double result_total = 0.0;
  for (int i = 0; i < nu::dmn_size(); i++) {
    result(i) += 1. - MOMS.G_r_t(i, i, r_DCA::parameter_type::origin_index(), 0);
    MOMS.orbital_occupancy(i) = result(i);
    result_total += result(i);
  }
  return result_total;
}

/*!
 *  We assume that G_ii(w>>0) ~ 1/(i w_m + A + B/(i w_m))
 */
template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::compute_density_correction(
    func::function<double, nu>& result) {
  std::complex<double> I(0, 1);

  double N_k = k_DCA::dmn_size();
  double beta = parameters.get_beta();

  func::function<double, func::dmn_variadic<nu, k_DCA>> A;
  func::function<double, func::dmn_variadic<nu, k_DCA>> B;

  func::function<double, func::dmn_variadic<nu, k_DCA>> A0;
  func::function<double, func::dmn_variadic<nu, k_DCA>> B0;

  compute_density_coefficients(A, B, MOMS.G_k_w);
  compute_density_coefficients(A0, B0, MOMS.G0_k_w);

  for (int k_i = 0; k_i < k_DCA::dmn_size(); k_i++) {
    for (int nu_i = 0; nu_i < nu::dmn_size(); nu_i++) {
      double tmp = 0.0;
      double sum = 1.e-16;

      int l = w::dmn_size() / 2;

      do {
        std::complex<double> I_wn = (M_PI / beta) * (1 + 2 * l) * I;

        std::complex<double> G = 1. / (I_wn + A(nu_i, k_i) + B(nu_i, k_i) / I_wn);
        std::complex<double> G0 = 1. / (I_wn + A0(nu_i, k_i) + B0(nu_i, k_i) / I_wn);

        tmp = std::real(G - G0);
        sum += tmp;

        l += 1;
      } while (std::abs(tmp / sum) > 1.e-6 and l < 1.e6);

      result(nu_i) += sum;
    }
  }
  for (int nu_i = 0; nu_i < nu::dmn_size(); nu_i++)
    result(nu_i) *= (2. / (beta * N_k));
}

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::compute_density_coefficients(
    func::function<double, func::dmn_variadic<nu, k_DCA>>& A,
    func::function<double, func::dmn_variadic<nu, k_DCA>>& B,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& G) {
  A = 0;
  B = 0;

  int nb_wm = parameters.get_number_of_tail_frequencies();

  if (nb_wm > 0) {
    for (int k_i = 0; k_i < k_DCA::dmn_size(); k_i++) {
      for (int nu_i = 0; nu_i < nu::dmn_size(); nu_i++) {
        for (int w_i = w::dmn_size() - nb_wm; w_i < w::dmn_size(); w_i++) {
          double wm = w::get_elements()[w_i];

          A(nu_i, k_i) += std::real(1. / G(nu_i, nu_i, k_i, w_i));
          B(nu_i, k_i) += wm * (wm - std::imag(1. / G(nu_i, nu_i, k_i, w_i)));
        }
      }
    }

    A /= nb_wm;
    B /= nb_wm;

    if (nb_wm == 1) {
      std::complex<double> I(0, 1);

      for (int k_i = 0; k_i < k_DCA::dmn_size(); k_i++) {
        for (int nu_i = 0; nu_i < nu::dmn_size(); nu_i++) {
          int w_i = w::dmn_size() - 1;
          double wm = w::get_elements()[w_i];

          if (abs((G(nu_i, nu_i, k_i, w_i)) -
                  1. / (I * wm + A(nu_i, k_i) + B(nu_i, k_i) / (I * wm))) > 1.e-12)
            throw std::logic_error(__FUNCTION__);
        }
      }
    }
  }
}

template <typename parameters_type, typename MOMS_type, typename coarsegraining_type>
void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::print_bounds() {
  if (concurrency.id() == 0) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::cout << "\t";
    std::cout << "\t mu : " << lower_bound.first << " (n = " << lower_bound.second << ")";
    std::cout << "\t mu : " << upper_bound.first << " (n = " << upper_bound.second << ")\t";
    std::cout << dca::util::print_time() << "\n";
  }
}

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_UPDATE_CHEMICAL_POTENTIAL_HPP
