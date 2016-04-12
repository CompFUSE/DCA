//-*-C++-*-

#ifndef DCA_DECONVOLUTION_TP_H
#define DCA_DECONVOLUTION_TP_H

#include "phys_library/DCA+_step/lattice_mapping/deconvolution/deconvolution_routines.h"

namespace DCA {
/*! \ingroup LATTICE-MAPPING
 *
 *  \author Peter Staar
 *  \brief  This class implements the deconvolution in the lattice-mapping.
 *
 */
template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class deconvolution_tp
    : public deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t> {
  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef typename source_k_dmn_t::parameter_type source_k_cluster_type;
  typedef typename target_k_dmn_t::parameter_type target_k_cluster_type;

  typedef typename source_k_cluster_type::dual_type source_r_cluster_type;
  typedef typename target_k_cluster_type::dual_type target_r_cluster_type;

  typedef dmn_0<source_r_cluster_type> source_r_dmn_t;
  typedef dmn_0<target_r_cluster_type> target_r_dmn_t;

  typedef math_algorithms::functional_transforms::basis_transform<target_k_dmn_t, target_r_dmn_t>
      trafo_k_to_r_type;
  typedef math_algorithms::functional_transforms::basis_transform<target_r_dmn_t, target_k_dmn_t>
      trafo_r_to_k_type;

public:
  deconvolution_tp(parameters_type& parameters_ref);
  ~deconvolution_tp();

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
deconvolution_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::~deconvolution_tp() {}

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

  LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> phi_inv("phi_inv",
                                                                  std::pair<int, int>(N, N));

  this->compute_T_inv_matrix(parameters.get_singular_value_cut_off(), phi_inv);

  math_algorithms::functional_transforms::TRANSFORM<k_dmn_t, target_k_dmn_t>::execute_on_all(
      Gamma_lattice_interp, Gamma_lattice_deconv, phi_inv);
}
}

#endif
