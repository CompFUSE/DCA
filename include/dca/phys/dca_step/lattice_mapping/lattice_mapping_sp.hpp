// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the lattice mapping for single-particle functions.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_SP_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_SP_HPP

#include <complex>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/lattice_mapping/deconvolution/deconvolution_sp.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/interpolation_sp.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class lattice_mapping_sp {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using DCA_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;
  using host_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::LATTICE_SP,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_HOST = func::dmn_0<host_k_cluster_type>;

public:
  lattice_mapping_sp(parameters_type& parameters_ref);

  void execute(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn_t, w>>& f_source,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& Sigma_interp,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& Sigma_deconv,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target);

  template <typename MOMS_type, typename HTS_solver_type, typename coarsegraining_sp_type>
  void execute_with_HTS_approximation(
      MOMS_type& HTS_MOMS, HTS_solver_type& HTS_solver, coarsegraining_sp_type& cluster_mapping_obj,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn_t, w>>& f_source,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& Sigma_interp,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& Sigma_deconv,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target);

private:
  template <typename k_dmn_t>
  void plot_function(func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& f);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  interpolation_sp<parameters_type, source_k_dmn_t, target_k_dmn_t> interpolation_obj;
  deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t> deconvolution_obj;

  func::function<double, nu> Sigma_shift;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> Sigma_cluster;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST, w>> Sigma_lattice;
};

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
lattice_mapping_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::lattice_mapping_sp(
    parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      interpolation_obj(parameters),
      deconvolution_obj(parameters),

      Sigma_shift("Sigma_shift"),

      Sigma_cluster("Sigma_cluster_HTS"),
      Sigma_lattice("Sigma_lattice_HTS") {}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void lattice_mapping_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn_t, w>>& f_source,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_interp,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_approx,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target) {
  symmetrize::execute(f_source);

  // plot_function(f_source);

  interpolation_obj.execute_with_alpha_transformation(f_source, f_interp);

  // plot_function(f_interp);

  symmetrize::execute(f_interp);

  deconvolution_obj.execute(f_source, f_interp, f_approx, f_target);

  // plot_function(f_target);

  symmetrize::execute(f_target);
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename MOMS_type, typename HTS_solver_type, typename coarsegraining_sp_type>
void lattice_mapping_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute_with_HTS_approximation(
    MOMS_type& HTS_MOMS, HTS_solver_type& HTS_solver, coarsegraining_sp_type& cluster_mapping_obj,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn_t, w>>& f_source,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_interp,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_approx,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target) {
  {
    HTS_solver.initialize(0);

    HTS_solver.execute();

    HTS_solver.finalize();

    HTS_MOMS.Sigma = f_source;

    cluster_mapping_obj.compute_S_K_w(HTS_MOMS.Sigma_lattice, HTS_MOMS.Sigma_cluster);

    Sigma_cluster = HTS_MOMS.Sigma_cluster;
    Sigma_lattice = HTS_MOMS.Sigma_lattice;

    //       HTS_MOMS.write("data_HTS.json");
  }

  plot_function(Sigma_lattice);

  {
    f_source -= Sigma_cluster;

    execute(f_source, f_interp, f_approx, f_target);

    f_source += Sigma_cluster;
    f_interp += Sigma_lattice;
    f_approx += Sigma_lattice;
    f_target += Sigma_lattice;
  }

  plot_function(f_interp);
  plot_function(f_target);
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename k_dmn_t>
void lattice_mapping_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::plot_function(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& f) {
  std::vector<double> x(0);
  std::vector<double> y(0);

  std::vector<double> z_re(0);
  std::vector<double> z_im(0);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
    x.push_back(k_dmn_t::get_elements()[k_ind][0]);
    y.push_back(k_dmn_t::get_elements()[k_ind][1]);

    z_re.push_back(std::real(f(0, 0, k_ind, w::dmn_size() / 2)));
    z_im.push_back(std::imag(f(0, 0, k_ind, w::dmn_size() / 2)));
  }

  util::Plot::heatMap(x, y, z_re, f.get_name());
  util::Plot::heatMap(x, y, z_im, f.get_name());
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_SP_HPP
