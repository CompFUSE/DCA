// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the lattice mapping for two-particle functions.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_TP_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_TP_HPP

#include <complex>
#include <iostream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/lattice_mapping/deconvolution/deconvolution_tp.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/interpolation_tp.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/util/plot.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class lattice_mapping_tp {
  struct tmp_cluster_domain {
    typedef typename target_k_dmn_t::parameter_type::element_type element_type;
    typedef typename target_k_dmn_t::parameter_type::dmn_specifications_type dmn_specifications_type;

    static int get_size() {
      return target_k_dmn_t::dmn_size();
    }
    static std::vector<element_type>& get_elements() {
      return target_k_dmn_t::get_elements();
    }
  };

public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using w = func::dmn_0<domains::frequency_domain>;
  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;

  using tmp_k_dmn_t = func::dmn_0<tmp_cluster_domain>;
  using tmp_vector_dmn_t = func::dmn_variadic<b, b, tmp_k_dmn_t, w_VERTEX>;
  using source_vector_dmn_t = func::dmn_variadic<b, b, source_k_dmn_t, w_VERTEX>;
  using target_vector_dmn_t = func::dmn_variadic<b, b, target_k_dmn_t, w_VERTEX>;

public:
  lattice_mapping_tp(parameters_type& parameters_ref);

  template <typename scalartype>
  void execute(func::function<std::complex<scalartype>,
                              func::dmn_variadic<source_vector_dmn_t, source_vector_dmn_t>>& f_source,
               func::function<std::complex<scalartype>,
                              func::dmn_variadic<target_vector_dmn_t, target_vector_dmn_t>>& f_target);

private:
  template <typename k_dmn_t>
  void plot_function(func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& f);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  interpolation_tp<parameters_type, source_k_dmn_t, target_k_dmn_t> interpolation_obj;
  deconvolution_tp<parameters_type, source_k_dmn_t, target_k_dmn_t> deconvolution_obj;
};

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
lattice_mapping_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::lattice_mapping_tp(
    parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      interpolation_obj(parameters),
      deconvolution_obj(parameters) {}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename scalartype>
void lattice_mapping_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute(
    func::function<std::complex<scalartype>,
                   func::dmn_variadic<source_vector_dmn_t, source_vector_dmn_t>>& f_source,
    func::function<std::complex<scalartype>,
                   func::dmn_variadic<target_vector_dmn_t, target_vector_dmn_t>>& f_target) {
  func::function<std::complex<scalartype>, func::dmn_variadic<tmp_vector_dmn_t, tmp_vector_dmn_t>> f_interp(
      "f_interp");

  {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\n start tp-interpolation of Gamma \n\n";

    interpolation_obj.execute(f_source, f_target);
  }

  {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\n start tp-deconvolution of Gamma \n\n";

    for (int i = 0; i < f_target.size(); i++)
      f_interp(i) = f_target(i);

    deconvolution_obj.execute(f_interp, f_target);
  }
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename k_dmn_t>
void lattice_mapping_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::plot_function(
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

  util::Plot::heatMap(x, y, z_re);
  util::Plot::heatMap(x, y, z_im);
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_TP_HPP
