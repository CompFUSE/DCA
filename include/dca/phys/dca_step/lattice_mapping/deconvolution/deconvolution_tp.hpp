// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the deconvolution step of the lattice mapping for two-particle functions.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_TP_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_TP_HPP

#include <complex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/dca_step/lattice_mapping/deconvolution/deconvolution_routines.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class deconvolution_tp
    : public deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t> {
  using concurrency_type = typename parameters_type::concurrency_type;

  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using host_vertex_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::LATTICE_TP,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_HOST_VERTEX = func::dmn_0<host_vertex_k_cluster_type>;

public:
  deconvolution_tp(parameters_type& parameters_ref);

  template <typename k_dmn_t, typename scalartype>
  void execute(func::function<std::complex<scalartype>,
                              func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_VERTEX>,
                                                 func::dmn_variadic<b, b, k_dmn_t, w_VERTEX>>>&
                   Gamma_lattice_interp,
               func::function<std::complex<scalartype>,
                              func::dmn_variadic<func::dmn_variadic<b, b, target_k_dmn_t, w_VERTEX>,
                                                 func::dmn_variadic<b, b, target_k_dmn_t, w_VERTEX>>>&
                   Gamma_lattice_deconv);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;
};

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
deconvolution_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::deconvolution_tp(
    parameters_type& parameters_ref)
    : deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>(parameters_ref),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()) {}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename k_dmn_t, typename scalartype>
void deconvolution_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute(
    func::function<std::complex<scalartype>,
                   func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_VERTEX>,
                                      func::dmn_variadic<b, b, k_dmn_t, w_VERTEX>>>& Gamma_lattice_interp,
    func::function<std::complex<scalartype>,
                   func::dmn_variadic<func::dmn_variadic<b, b, target_k_dmn_t, w_VERTEX>,
                                      func::dmn_variadic<b, b, target_k_dmn_t, w_VERTEX>>>&
        Gamma_lattice_deconv) {
  int N = k_HOST_VERTEX::dmn_size();

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> phi_inv(
      "phi_inv", std::pair<int, int>(N, N));

  this->compute_T_inv_matrix(parameters.get_Gamma_deconvolution_cut_off(), phi_inv);

  math::transform::FunctionTransform<k_dmn_t, target_k_dmn_t>::execute_on_all(
      Gamma_lattice_interp, Gamma_lattice_deconv, phi_inv);
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_TP_HPP
