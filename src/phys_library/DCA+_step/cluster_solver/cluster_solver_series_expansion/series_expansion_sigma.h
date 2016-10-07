// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class computes the self-energy using a series expansion.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_SERIES_EXPANSION_SIGMA_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_SERIES_EXPANSION_SIGMA_H

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/1st_order_perturbation_sigma.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/2nd_order_perturbation_sigma.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/3rd_order_perturbation_sigma.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/4th_order_perturbation_sigma.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/compute_bare_bubble.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/compute_interaction.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/compute_lattice_Greens_function.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"

namespace DCA {
namespace SERIES_EXPANSION {

template <class parameters_type, class MOMS_type>
class series_expansion {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using k_HOST = func::dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                            LATTICE_SP, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;
  using k_dmn_t = k_HOST;

  using w = func::dmn_0<frequency_domain>;
  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using sigma_function_t =
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>;

public:
  series_expansion(parameters_type& parameter_ref, MOMS_type& MOMS_ref);

  template <class stream_type>
  void to_JSON(stream_type& ss);

  void execute(bool do_not_adjust_mu = true);

  sigma_function_t& get_Sigma() {
    return Sigma;
  }

  template <typename Writer>
  void write(Writer& writer);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;
  MOMS_type& MOMS;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>> Sigma;

  compute_interaction interaction_obj;

  compute_bubble<ph, parameters_type, k_dmn_t, w> ph_bubble;
  compute_bubble<pp, parameters_type, k_dmn_t, w> pp_bubble;

  sigma_perturbation<1, parameters_type, k_dmn_t> sigma_perturbation_1_obj;
  sigma_perturbation<2, parameters_type, k_dmn_t> sigma_perturbation_2_obj;
  sigma_perturbation<3, parameters_type, k_dmn_t> sigma_perturbation_3_obj;
  sigma_perturbation<4, parameters_type, k_dmn_t> sigma_perturbation_4_obj;
};

template <class parameters_type, class MOMS_type>
series_expansion<parameters_type, MOMS_type>::series_expansion(parameters_type& parameters_ref,
                                                               MOMS_type& MOMS_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),
      MOMS(MOMS_ref),

      Sigma("perturbation-Sigma"),

      interaction_obj(),

      ph_bubble(parameters),
      pp_bubble(parameters),

      sigma_perturbation_1_obj(parameters, interaction_obj),
      sigma_perturbation_2_obj(parameters, interaction_obj, ph_bubble, pp_bubble),
      sigma_perturbation_3_obj(parameters, interaction_obj, ph_bubble, pp_bubble),
      sigma_perturbation_4_obj(parameters, interaction_obj, ph_bubble, pp_bubble) {
  interaction_obj.execute(MOMS.H_interactions);
}

template <class parameters_type, class MOMS_type>
void series_expansion<parameters_type, MOMS_type>::execute(bool /*do_not_adjust_mu*/) {
  compute_lattice_Greens_function<parameters_type, MOMS_type, k_dmn_t, w>
      compute_lattice_Greens_function_obj(parameters, MOMS);

  compute_lattice_Greens_function_obj.execute();

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST, w>>& G_k_w =
      compute_lattice_Greens_function_obj.get_G_k_w();

  ph_bubble.threaded_execute_on_cluster(G_k_w);
  pp_bubble.threaded_execute_on_cluster(G_k_w);

  sigma_perturbation_1_obj.execute_on_cluster(MOMS.orbital_occupancy);
  sigma_perturbation_2_obj.threaded_execute_on_cluster(G_k_w);

  Sigma = 0.;

  Sigma += sigma_perturbation_1_obj.get_function();
  Sigma += sigma_perturbation_2_obj.get_function();

  if (true) {
    std::complex<double> I(0, 1);
    for (int b_ind = 0; b_ind < 2 * b::dmn_size(); ++b_ind) {
      for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
        int wc_ind = w::dmn_size() / 8;

        double wc = w::get_elements()[wc_ind];

        std::complex<double> Sigma_wc = Sigma(b_ind, b_ind, k_ind, wc_ind);

        double alpha = real(Sigma_wc);
        double beta = imag(Sigma_wc * wc);

        for (int w_ind = 0; w_ind < wc_ind; ++w_ind) {
          Sigma(b_ind, b_ind, k_ind, w_ind) = alpha + beta * I / w::get_elements()[w_ind];
          Sigma(b_ind, b_ind, k_ind, w::dmn_size() - 1 - w_ind) =
              alpha - beta * I / w::get_elements()[w_ind];
        }
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
template <typename Writer>
void series_expansion<parameters_type, MOMS_type>::write(Writer& writer) {
  writer.execute(Sigma);

  ph_bubble.write(writer);
  pp_bubble.write(writer);

  sigma_perturbation_1_obj.write(writer);
  sigma_perturbation_2_obj.write(writer);
  sigma_perturbation_3_obj.write(writer);
  sigma_perturbation_4_obj.write(writer);
}
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_SERIES_EXPANSION_SIGMA_H
