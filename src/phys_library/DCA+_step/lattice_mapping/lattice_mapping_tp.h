// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the lattice mapping for two-particle functions.

#ifndef PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_TP_H
#define PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_TP_H

#include <complex>
#include <iostream>
#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "comp_library/function_plotting/include_plotting.h"
#include "phys_library/DCA+_step/lattice_mapping/deconvolution/deconvolution_tp.h"
#include "phys_library/DCA+_step/lattice_mapping/interpolation/interpolation_tp.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"

namespace DCA {

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

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index

  using w = dmn_0<frequency_domain>;
  using compact_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::COMPACT>;
  using w_VERTEX = dmn_0<compact_vertex_frequency_domain_type>;

  using tmp_k_dmn_t = dmn_0<tmp_cluster_domain>;
  using tmp_vector_dmn_t = dmn_4<b, b, tmp_k_dmn_t, w_VERTEX>;
  using source_vector_dmn_t = dmn_4<b, b, source_k_dmn_t, w_VERTEX>;
  using target_vector_dmn_t = dmn_4<b, b, target_k_dmn_t, w_VERTEX>;

public:
  lattice_mapping_tp(parameters_type& parameters_ref);

  template <typename scalartype>
  void execute(FUNC_LIB::function<std::complex<scalartype>,
                                  dmn_2<source_vector_dmn_t, source_vector_dmn_t>>& f_source,
               FUNC_LIB::function<std::complex<scalartype>,
                                  dmn_2<target_vector_dmn_t, target_vector_dmn_t>>& f_target);

private:
  void initialize();

  template <typename k_dmn_t>
  void plot_function(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_dmn_t, w>>& f);

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
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<source_vector_dmn_t, source_vector_dmn_t>>& f_source,
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<target_vector_dmn_t, target_vector_dmn_t>>&
        f_target) {
  FUNC_LIB::function<std::complex<scalartype>, dmn_2<tmp_vector_dmn_t, tmp_vector_dmn_t>> f_interp(
      "f_interp");

  {
    if (concurrency.id() == 0)
      std::cout << "\n\n start tp-interpolation of Gamma \n\n";

    interpolation_obj.execute(f_source, f_target);
  }

  {
    if (concurrency.id() == 0)
      std::cout << "\n\n start tp-deconvolution of Gamma \n\n";

    for (int i = 0; i < f_target.size(); i++)
      f_interp(i) = f_target(i);

    deconvolution_obj.execute(f_interp, f_target);
  }
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename k_dmn_t>
void lattice_mapping_tp<parameters_type, source_k_dmn_t, target_k_dmn_t>::plot_function(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_dmn_t, w>>& f) {
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

  SHOW::heatmap(x, y, z_re);
  SHOW::heatmap(x, y, z_im);
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_LATTICE_MAPPING_TP_H
