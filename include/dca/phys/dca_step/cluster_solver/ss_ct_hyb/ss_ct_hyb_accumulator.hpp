// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// This class organizes the measurements in the SS CT-HYB QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_ACCUMULATOR_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/feynman_expansion_order_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/mc_accumulator_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/accumulator/sp/sp_accumulator_nfft.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_walker.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_hybridization_solver_routines.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
class SsCtHybAccumulator : public MC_accumulator_data,
                           public ss_hybridization_solver_routines<parameters_type, Data> {
public:
  using this_type = SsCtHybAccumulator<device_t, parameters_type, Data>;

  using ParametersType = parameters_type;
  using DataType = Data;

  typedef SsCtHybWalker<device_t, parameters_type, Data> walker_type;

  typedef ss_hybridization_solver_routines<parameters_type, Data> ss_hybridization_solver_routines_type;

  typedef
      typename walker_type::ss_hybridization_walker_routines_type ss_hybridization_walker_routines_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = func::dmn_variadic<nu, nu>;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  typedef RClusterDmn r_dmn_t;

  typedef func::dmn_variadic<nu, nu, r_dmn_t> p_dmn_t;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef double scalar_type;

  typedef
      typename SsCtHybTypedefs<parameters_type, Data>::vertex_vertex_matrix_type vertex_vertex_matrix_type;
  typedef
      typename SsCtHybTypedefs<parameters_type, Data>::orbital_configuration_type orbital_configuration_type;

  typedef typename SsCtHybTypedefs<parameters_type, Data>::configuration_type configuration_type;

  typedef func::function<vertex_vertex_matrix_type, nu> M_matrix_type;

public:
  SsCtHybAccumulator(const parameters_type& parameters_ref, Data& data_ref, int id = 0);

  void initialize(int dca_iteration);

  void finalize();  // func::function<double, nu> mu_DC);

  void updateFrom(walker_type& walker);
  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sumTo(this_type& other);

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

  //  TODO: implement.
  std::size_t deviceFingerprint() const {
    return 0;
  }
  static std::size_t staticDeviceFingerprint() {
    return 0;
  }

protected:
  using MC_accumulator_data::DCA_iteration;
  using MC_accumulator_data::number_of_measurements;

  using MC_accumulator_data::current_sign;
  using MC_accumulator_data::accumulated_sign;

  const parameters_type& parameters_;
  Data& data_;
  const concurrency_type& concurrency;

  int thread_id;

  configuration_type configuration;
  func::function<vertex_vertex_matrix_type, nu> M_matrices;

  func::function<double, func::dmn_0<ctaux::Feynman_expansion_order_domain>> visited_expansion_order_k;

  func::function<double, nu> length;
  func::function<double, func::dmn_variadic<nu, nu>> overlap;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>> G_r_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>> GS_r_w;

  SpAccumulatorNfft<parameters_type, Data> single_particle_accumulator_obj;

  bool finalized_;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
SsCtHybAccumulator<device_t, parameters_type, Data>::SsCtHybAccumulator(const parameters_type& parameters_ref,
                                                                        Data& data_ref, int id)
    : ss_hybridization_solver_routines<parameters_type, Data>(parameters_ref, data_ref),

      parameters_(parameters_ref),
      data_(data_ref),
      concurrency(parameters_.get_concurrency()),

      thread_id(id),

      configuration(),
      M_matrices("accumulator-M-matrices"),

      visited_expansion_order_k("visited-expansion-order-k"),

      length("length"),
      overlap("overlap"),

      G_r_w("G-r-w-measured"),
      GS_r_w("GS-r-w-measured"),

      single_particle_accumulator_obj(parameters_),
      finalized_(false) {}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybAccumulator<device_t, parameters_type, Data>::initialize(int dca_iteration) {
  MC_accumulator_data::initialize(dca_iteration);

  visited_expansion_order_k = 0;

  single_particle_accumulator_obj.initialize(G_r_w, GS_r_w);

  length = 0;
  overlap = 0;

  finalized_ = false;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybAccumulator<device_t, parameters_type,
                        Data>::finalize()  // func::function<double, nu> mu_DC)
{
  if (finalized_)
    return;
  single_particle_accumulator_obj.finalize(G_r_w, GS_r_w);
  finalized_ = true;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
template <typename Writer>
void SsCtHybAccumulator<device_t, parameters_type, Data>::write(Writer& writer) {
  writer.execute(G_r_w);
  writer.execute(GS_r_w);
}

/*************************************************************
 **                                                         **
 **                    G2 - MEASUREMENTS                    **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybAccumulator<device_t, parameters_type, Data>::updateFrom(walker_type& walker) {
  current_sign = walker.get_sign();

  configuration.copy_from(walker.get_configuration());

  for (int l = 0; l < nu::dmn_size(); l++)
    M_matrices(l) = walker.get_M_matrices()(l);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybAccumulator<device_t, parameters_type, Data>::measure() {
  number_of_measurements += 1;
  accumulated_sign += current_sign;

  int k = configuration.size();
  if (k < visited_expansion_order_k.size())
    visited_expansion_order_k(k) += 1;

  single_particle_accumulator_obj.accumulate(current_sign, configuration, M_matrices,
                                             data_.H_interactions);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybAccumulator<device_t, parameters_type, Data>::accumulate_length(walker_type& walker) {
  ss_hybridization_walker_routines_type& hybridization_routines =
      walker.get_ss_hybridization_walker_routines();

  Hybridization_vertex full_segment(0, parameters_.get_beta());

  for (int ind = 0; ind < b::dmn_size() * s::dmn_size(); ind++) {
    length(ind) += hybridization_routines.compute_overlap(
        full_segment, walker.get_configuration().get_vertices(ind),
        walker.get_configuration().get_full_line(ind), parameters_.get_beta());
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybAccumulator<device_t, parameters_type, Data>::accumulate_overlap(walker_type& walker) {
  ss_hybridization_walker_routines_type& hybridization_routines =
      walker.get_ss_hybridization_walker_routines();

  Hybridization_vertex full_segment(0, parameters_.get_beta());

  for (int ind_1 = 0; ind_1 < b::dmn_size() * s::dmn_size(); ind_1++) {
    for (int ind_2 = 0; ind_2 < b::dmn_size() * s::dmn_size(); ind_2++) {
      if (walker.get_configuration().get_full_line(ind_1)) {
        overlap(ind_1, ind_2) += hybridization_routines.compute_overlap(
            full_segment, walker.get_configuration().get_vertices(ind_2),
            walker.get_configuration().get_full_line(ind_2), parameters_.get_beta());
      }
      else {
        for (typename orbital_configuration_type::iterator it =
                 walker.get_configuration().get_vertices(ind_1).begin();
             it != walker.get_configuration().get_vertices(ind_1).end(); it++) {
          overlap(ind_1, ind_2) += hybridization_routines.compute_overlap(
              *it, walker.get_configuration().get_vertices(ind_2),
              walker.get_configuration().get_full_line(ind_2), parameters_.get_beta());
        }
      }
    }
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybAccumulator<device_t, parameters_type, Data>::sumTo(this_type& other) {
  other.accumulated_sign += accumulated_sign;
  other.number_of_measurements += number_of_measurements;

  other.get_visited_expansion_order_k() += visited_expansion_order_k;

  single_particle_accumulator_obj.sumTo(other.single_particle_accumulator_obj);
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_ACCUMULATOR_HPP
