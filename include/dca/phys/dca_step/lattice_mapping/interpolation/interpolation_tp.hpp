// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the interpolated cluster vertex.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_TP_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_TP_HPP

#include <complex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/interpolation_routines.hpp"
#include "dca/phys/domains/cluster/centered_cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
class interpolation_tp : public interpolation_routines<parameters_type, source_k_dmn, target_k_dmn> {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using source_r_cluster_type = typename source_k_dmn::parameter_type::dual_type;
  using r_centered_dmn = func::dmn_0<domains::centered_cluster_domain<source_r_cluster_type>>;

  using basis_function_type = math::transform::basis_function<
      typename target_k_dmn::parameter_type, math::transform::KRONECKER_DELTA,
      typename source_k_dmn::parameter_type, math::transform::HERMITE_CUBIC_SPLINE>;

  using b = func::dmn_0<domains::electron_band_domain>;

  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;

  using DCA_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;
  using host_vertex_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::LATTICE_TP,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_HOST_VERTEX = func::dmn_0<host_vertex_k_cluster_type>;

public:
  interpolation_tp(parameters_type& parameters_ref);

  template <typename scalartype>
  void initialize_T_K_to_k(dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& T_K_to_k);

  template <typename scalartype>
  void execute(
      func::function<std::complex<scalartype>,
                     func::dmn_variadic<func::dmn_variadic<b, b, k_DCA, w_VERTEX>,
                                        func::dmn_variadic<b, b, k_DCA, w_VERTEX>>>& Gamma_cluster,
      func::function<std::complex<scalartype>,
                     func::dmn_variadic<func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>,
                                        func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>>>& Gamma_lattice);

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
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& T_K_to_k) {
  int Nr = target_k_dmn::dmn_size();
  int Nc = source_k_dmn::dmn_size();

  T_K_to_k.resizeNoCopy(std::pair<int, int>(Nr, Nc));

  for (int j = 0; j < Nc; j++)
    for (int i = 0; i < Nr; i++)
      T_K_to_k(i, j) = basis_function_type::execute(i, j);
}

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
template <typename scalartype>
void interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::execute(
    func::function<std::complex<scalartype>,
                   func::dmn_variadic<func::dmn_variadic<b, b, k_DCA, w_VERTEX>,
                                      func::dmn_variadic<b, b, k_DCA, w_VERTEX>>>& Gamma_cluster,
    func::function<std::complex<scalartype>,
                   func::dmn_variadic<func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>,
                                      func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>>>& Gamma_lattice) {
  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> T_K_to_k("T_K_to_k");

  initialize_T_K_to_k(T_K_to_k);

  math::transform::FunctionTransform<k_DCA, k_HOST_VERTEX>::execute_on_all(Gamma_cluster,
                                                                           Gamma_lattice, T_K_to_k);
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_TP_HPP
