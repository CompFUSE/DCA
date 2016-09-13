// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the deconvolution step of the lattice mapping for two-particle functions.

#ifndef PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_TP_H
#define PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_TP_H

#include <complex>
#include <utility>

#include "comp_library/function_library/include_function_library.h"
#include "math_library/functional_transforms/function_transforms/function_transforms.hpp"
#include "phys_library/DCA+_step/lattice_mapping/deconvolution/deconvolution_routines.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"

namespace DCA {

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class deconvolution_tp
    : public deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t> {
  using concurrency_type = typename parameters_type::concurrency_type;

  using compact_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::COMPACT>;
  using w_VERTEX = dmn_0<compact_vertex_frequency_domain_type>;
  using b = dmn_0<electron_band_domain>;
  using host_vertex_k_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                                    LATTICE_TP, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
  using k_HOST_VERTEX = dmn_0<host_vertex_k_cluster_type>;

public:
  deconvolution_tp(parameters_type& parameters_ref);

  template <typename k_dmn_t, typename scalartype>
  void execute(
      FUNC_LIB::function<std::complex<scalartype>,
                         dmn_2<dmn_4<b, b, k_dmn_t, w_VERTEX>, dmn_4<b, b, k_dmn_t, w_VERTEX>>>&
          Gamma_lattice_interp,
      FUNC_LIB::function<std::complex<scalartype>, dmn_2<dmn_4<b, b, target_k_dmn_t, w_VERTEX>,
                                                         dmn_4<b, b, target_k_dmn_t, w_VERTEX>>>&
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
    FUNC_LIB::function<std::complex<scalartype>,
                       dmn_2<dmn_4<b, b, k_dmn_t, w_VERTEX>, dmn_4<b, b, k_dmn_t, w_VERTEX>>>&
        Gamma_lattice_interp,
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<dmn_4<b, b, target_k_dmn_t, w_VERTEX>,
                                                       dmn_4<b, b, target_k_dmn_t, w_VERTEX>>>&
        Gamma_lattice_deconv) {
  int N = k_HOST_VERTEX::dmn_size();

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> phi_inv(
      "phi_inv", std::pair<int, int>(N, N));

  this->compute_T_inv_matrix(parameters.get_singular_value_cut_off(), phi_inv);

  math_algorithms::functional_transforms::TRANSFORM<k_dmn_t, target_k_dmn_t>::execute_on_all(
      Gamma_lattice_interp, Gamma_lattice_deconv, phi_inv);
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_TP_H
