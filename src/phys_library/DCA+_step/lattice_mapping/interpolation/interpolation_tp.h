// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the interpolated cluster vertex.

#ifndef PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_TP_H
#define PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_TP_H

#include <complex>
#include <utility>

#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg.hpp"
#include "math_library/functional_transforms/function_transforms/function_transforms.hpp"
#include "phys_library/DCA+_step/lattice_mapping/interpolation/interpolation_routines.h"
#include "phys_library/domains/cluster/centered_cluster_domain.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"

namespace DCA {

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
class interpolation_tp : public interpolation_routines<parameters_type, source_k_dmn, target_k_dmn> {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using source_r_cluster_type = typename source_k_dmn::parameter_type::dual_type;
  using r_centered_dmn = dmn_0<centered_cluster_domain<source_r_cluster_type>>;

  using basis_function_type = math_algorithms::functional_transforms::basis_function<
      typename target_k_dmn::parameter_type, math_algorithms::KRONECKER_DELTA,
      typename source_k_dmn::parameter_type, math_algorithms::HERMITE_CUBIC_SPLINE>;

  using b = dmn_0<electron_band_domain>;

  using compact_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::COMPACT>;
  using w_VERTEX = dmn_0<compact_vertex_frequency_domain_type>;

  using DCA_k_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                            CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
  using k_DCA = dmn_0<DCA_k_cluster_type>;
  using host_vertex_k_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                                    LATTICE_TP, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
  using k_HOST_VERTEX = dmn_0<host_vertex_k_cluster_type>;

public:
  interpolation_tp(parameters_type& parameters_ref);

  template <typename scalartype>
  void initialize_T_K_to_k(LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& T_K_to_k);

  template <typename scalartype>
  void execute(FUNC_LIB::function<std::complex<scalartype>,
                                  dmn_2<dmn_4<b, b, k_DCA, w_VERTEX>, dmn_4<b, b, k_DCA, w_VERTEX>>>&
                   Gamma_cluster,
               FUNC_LIB::function<std::complex<scalartype>,
                                  dmn_2<dmn_4<b, b, k_HOST_VERTEX, w_VERTEX>,
                                        dmn_4<b, b, k_HOST_VERTEX, w_VERTEX>>>& Gamma_lattice);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;
};

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::interpolation_tp(
    parameters_type& parameters_ref)
    : interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>(parameters_ref),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()) {}

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
template <typename scalartype>
void interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::initialize_T_K_to_k(
    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& T_K_to_k) {
  int Nr = target_k_dmn::dmn_size();
  int Nc = source_k_dmn::dmn_size();

  T_K_to_k.resize_no_copy(std::pair<int, int>(Nr, Nc));

  for (int j = 0; j < Nc; j++)
    for (int i = 0; i < Nr; i++)
      T_K_to_k(i, j) = basis_function_type::execute(i, j);
}

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
template <typename scalartype>
void interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::execute(
    FUNC_LIB::function<std::complex<scalartype>,
                       dmn_2<dmn_4<b, b, k_DCA, w_VERTEX>, dmn_4<b, b, k_DCA, w_VERTEX>>>& Gamma_cluster,
    FUNC_LIB::function<std::complex<scalartype>,
                       dmn_2<dmn_4<b, b, k_HOST_VERTEX, w_VERTEX>, dmn_4<b, b, k_HOST_VERTEX, w_VERTEX>>>&
        Gamma_lattice) {
  LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> T_K_to_k("T_K_to_k");

  initialize_T_K_to_k(T_K_to_k);

  math_algorithms::functional_transforms::TRANSFORM<k_DCA, k_HOST_VERTEX>::execute_on_all(
      Gamma_cluster, Gamma_lattice, T_K_to_k);
}

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_TP_H
