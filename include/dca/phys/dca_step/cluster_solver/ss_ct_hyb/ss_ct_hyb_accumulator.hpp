// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Bart Ydens
//
// This class organizes the measurements in the SS CT-HYB QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_ACCUMULATOR_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/feynman_expansion_order_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/mc_accumulator_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/accumulator/sp/sp_accumulator_nfft.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_walker.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_hybridization_solver_routines.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class SsCtHybAccumulator : public MC_accumulator_data,
                           public ss_hybridization_solver_routines<parameters_type, MOMS_type> {
public:
  using this_type = SsCtHybAccumulator<device_t, parameters_type, MOMS_type>;

  typedef parameters_type my_parameters_type;
  typedef MOMS_type my_MOMS_type;

  typedef SsCtHybWalker<device_t, parameters_type, MOMS_type> walker_type;

  typedef ss_hybridization_solver_routines<parameters_type, MOMS_type> ss_hybridization_solver_routines_type;

  typedef
      typename walker_type::ss_hybridization_walker_routines_type ss_hybridization_walker_routines_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = func::dmn_variadic<nu, nu>;

  using r_DCA =
      func::dmn_0<domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;
  typedef r_DCA r_dmn_t;

  typedef func::dmn_variadic<nu, nu, r_dmn_t> p_dmn_t;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef double scalar_type;

  typedef typename SsCtHybTypedefs<parameters_type, MOMS_type>::vertex_vertex_matrix_type
      vertex_vertex_matrix_type;
  typedef typename SsCtHybTypedefs<parameters_type, MOMS_type>::orbital_configuration_type
      orbital_configuration_type;

  typedef typename SsCtHybTypedefs<parameters_type, MOMS_type>::configuration_type configuration_type;

  typedef func::function<vertex_vertex_matrix_type, nu> M_matrix_type;

public:
  SsCtHybAccumulator(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id = 0);

  void initialize(int dca_iteration);

  void finalize();  // func::function<double, nu> mu_DC);

  void update_from(walker_type& walker);
  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sum_to(this_type& other);

  configuration_type& get_configuration() {
    return configuration;
  }

  func::function<double, func::dmn_0<ctaux::Feynman_expansion_order_domain>>& get_visited_expansion_order_k() {
    return visited_expansion_order_k;
  }

  const auto& get_G_r_w() const {
    return G_r_w;
  }

  // TODO: Remove getter methods that return a non-const reference.
  auto& get_G_r_w() {
    return G_r_w;
  }
  const auto& get_GS_r_w() const {
    return GS_r_w;
  }
  auto& get_GS_r_w() {
    return GS_r_w;
  }

  void accumulate_length(walker_type& walker);
  void accumulate_overlap(walker_type& walker);

  func::function<double, nu>& get_length() {
    return length;
  }

  func::function<double, nu_nu>& get_overlap() {
    return overlap;
  }

  /*!
   *  \brief Print the functions G_r_w and G_k_w.
   */
  template <typename Writer>
  void write(Writer& writer);

protected:
  using MC_accumulator_data::DCA_iteration;
  using MC_accumulator_data::number_of_measurements;

  using MC_accumulator_data::current_sign;
  using MC_accumulator_data::accumulated_sign;

  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

  int thread_id;

  configuration_type configuration;
  func::function<vertex_vertex_matrix_type, nu> M_matrices;

  func::function<double, func::dmn_0<ctaux::Feynman_expansion_order_domain>> visited_expansion_order_k;

  func::function<double, nu> length;
  func::function<double, func::dmn_variadic<nu, nu>> overlap;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>> G_r_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>> GS_r_w;

  SpAccumulatorNfft<parameters_type, MOMS_type> single_particle_accumulator_obj;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::SsCtHybAccumulator(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id)
    : ss_hybridization_solver_routines<parameters_type, MOMS_type>(parameters_ref, MOMS_ref),

      parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thread_id(id),

      configuration(),
      M_matrices("accumulator-M-matrices"),

      visited_expansion_order_k("visited-expansion-order-k"),

      length("length"),
      overlap("overlap"),

      G_r_w("G-r-w-measured"),
      GS_r_w("GS-r-w-measured"),

      single_particle_accumulator_obj(parameters) {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::initialize(int dca_iteration) {
  MC_accumulator_data::initialize(dca_iteration);

  visited_expansion_order_k = 0;

  single_particle_accumulator_obj.initialize(G_r_w, GS_r_w);

  length = 0;
  overlap = 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybAccumulator<device_t, parameters_type,
                        MOMS_type>::finalize()  // func::function<double, nu> mu_DC)
{
  single_particle_accumulator_obj.finalize(G_r_w, GS_r_w);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename Writer>
void SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::write(Writer& writer) {
  writer.execute(G_r_w);
  writer.execute(GS_r_w);
}

/*************************************************************
 **                                                         **
 **                    G2 - MEASUREMENTS                    **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::update_from(walker_type& walker) {
  current_sign = walker.get_sign();

  configuration.copy_from(walker.get_configuration());

  for (int l = 0; l < nu::dmn_size(); l++)
    M_matrices(l) = walker.get_M_matrices()(l);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::measure() {
  number_of_measurements += 1;
  accumulated_sign += current_sign;

  int k = configuration.size();
  if (k < visited_expansion_order_k.size())
    visited_expansion_order_k(k) += 1;

  single_particle_accumulator_obj.accumulate(current_sign, configuration, M_matrices,
                                             MOMS.H_interactions);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::accumulate_length(walker_type& walker) {
  ss_hybridization_walker_routines_type& hybridization_routines =
      walker.get_ss_hybridization_walker_routines();

  Hybridization_vertex full_segment(0, parameters.get_beta());

  for (int ind = 0; ind < b::dmn_size() * s::dmn_size(); ind++) {
    length(ind) += hybridization_routines.compute_overlap(
        full_segment, walker.get_configuration().get_vertices(ind),
        walker.get_configuration().get_full_line(ind), parameters.get_beta());
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::accumulate_overlap(walker_type& walker) {
  ss_hybridization_walker_routines_type& hybridization_routines =
      walker.get_ss_hybridization_walker_routines();

  Hybridization_vertex full_segment(0, parameters.get_beta());

  for (int ind_1 = 0; ind_1 < b::dmn_size() * s::dmn_size(); ind_1++) {
    for (int ind_2 = 0; ind_2 < b::dmn_size() * s::dmn_size(); ind_2++) {
      if (walker.get_configuration().get_full_line(ind_1)) {
        overlap(ind_1, ind_2) += hybridization_routines.compute_overlap(
            full_segment, walker.get_configuration().get_vertices(ind_2),
            walker.get_configuration().get_full_line(ind_2), parameters.get_beta());
      }
      else {
        for (typename orbital_configuration_type::iterator it =
                 walker.get_configuration().get_vertices(ind_1).begin();
             it != walker.get_configuration().get_vertices(ind_1).end(); it++) {
          overlap(ind_1, ind_2) += hybridization_routines.compute_overlap(
              *it, walker.get_configuration().get_vertices(ind_2),
              walker.get_configuration().get_full_line(ind_2), parameters.get_beta());
        }
      }
    }
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybAccumulator<device_t, parameters_type, MOMS_type>::sum_to(this_type& other) {
  finalize();

  other.get_sign() += get_sign();
  other.get_number_of_measurements() += get_number_of_measurements();

  for (int i = 0; i < visited_expansion_order_k.size(); i++)
    other.get_visited_expansion_order_k()(i) += visited_expansion_order_k(i);

  {  // sp-measurements
    for (int i = 0; i < G_r_w.size(); i++)
      other.get_G_r_w()(i) += G_r_w(i);

    for (int i = 0; i < GS_r_w.size(); i++)
      other.get_GS_r_w()(i) += GS_r_w(i);
  }
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_ACCUMULATOR_HPP
