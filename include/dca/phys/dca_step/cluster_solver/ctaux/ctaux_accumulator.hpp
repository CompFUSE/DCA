// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class organizes the measurements in the CT-AUX QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_ACCUMULATOR_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/tp_equal_time_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/feynman_expansion_order_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ct_aux_hs_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_pair.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/mc_accumulator_data.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/quantum/numerical_error_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/four_point_type.hpp"
#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real = double>
class CtauxAccumulator : public MC_accumulator_data {
public:
  using this_type = CtauxAccumulator<device_t, Parameters, Data, Real>;

  using ParametersType = Parameters;
  using DataType = Data;

  typedef vertex_pair<Parameters> vertex_pair_type;
  typedef vertex_singleton vertex_singleton_type;

  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;
  using WVertexDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  typedef RClusterDmn r_dmn_t;
  typedef KClusterDmn k_dmn_t;

  typedef typename Parameters::profiler_type profiler_type;
  typedef typename Parameters::concurrency_type concurrency_type;

  typedef CT_AUX_HS_configuration<Parameters> configuration_type;

  CtauxAccumulator(Parameters& parameters_ref, Data& data_ref, int id);

  template <typename Writer>
  void write(Writer& writer);

  void initialize(int dca_iteration);

  template <typename walker_type>
  void updateFrom(walker_type& walker);

  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects
  // of the 'other' accumulator.
  void sumTo(this_type& other);

  void finalize();

  std::vector<vertex_singleton_type>& get_configuration(e_spin_states_type e_spin = e_UP);

  func::function<Real, func::dmn_0<domains::numerical_error_domain>>& get_error_distribution() {
    return error;
  }

  func::function<Real, func::dmn_0<Feynman_expansion_order_domain>>& get_visited_expansion_order_k() {
    return visited_expansion_order_k;
  }

  // equal time-measurements
  // TODO: Make equal time getters const.
  func::function<Real, func::dmn_variadic<nu, nu, r_dmn_t, t>>& get_G_r_t() {
    return equal_time_accumulator_ptr_->get_G_r_t();
  }
  func::function<Real, func::dmn_variadic<nu, nu, r_dmn_t, t>>& get_G_r_t_stddev() {
    return equal_time_accumulator_ptr_->get_G_r_t_stddev();
  }

  func::function<Real, func::dmn_variadic<b, r_dmn_t>>& get_charge_cluster_moment() {
    return equal_time_accumulator_ptr_->get_charge_cluster_moment();
  }
  func::function<Real, func::dmn_variadic<b, r_dmn_t>>& get_magnetic_cluster_moment() {
    return equal_time_accumulator_ptr_->get_magnetic_cluster_moment();
  }
  func::function<Real, func::dmn_variadic<b, r_dmn_t>>& get_dwave_pp_correlator() {
    return equal_time_accumulator_ptr_->get_dwave_pp_correlator();
  }

  // sp-measurements
  const auto& get_sign_times_M_r_w() const {
    return single_particle_accumulator_obj.get_sign_times_M_r_w();
  }

  const auto& get_sign_times_M_r_w_sqr() const {
    return single_particle_accumulator_obj.get_sign_times_M_r_w_sqr();
  }

  // tp-measurements
  const auto& get_sign_times_G4() {
    return two_particle_accumulator_.get_sign_times_G4();
  }

  bool compute_std_deviation() const {
    return compute_std_deviation_;
  }

#ifdef MEASURE_ERROR_BARS
  void store_standard_deviation(int nr_measurements, std::ofstream& points_file,
                                std::ofstream& norm_file);
  void update_sum_squares();
#endif

  std::size_t deviceFingerprint() const {
    return M_[0].deviceFingerprint() + M_[1].deviceFingerprint() +
           single_particle_accumulator_obj.deviceFingerprint() +
           two_particle_accumulator_.deviceFingerprint();
  }

  static std::size_t staticDeviceFingerprint() {
    return accumulator::TpAccumulator<Parameters, device_t>::staticDeviceFingerprint();
  }

private:
  void accumulate_single_particle_quantities();

  void accumulate_equal_time_quantities();
  void accumulate_equal_time_quantities(const std::array<linalg::Matrix<Real, linalg::GPU>, 2>& M);
  void accumulate_equal_time_quantities(const std::array<linalg::Matrix<Real, linalg::CPU>, 2>& M);

  void accumulate_two_particle_quantities();

protected:
  Parameters& parameters_;
  Data& data_;
  concurrency_type& concurrency;

  int thread_id;

  using MC_accumulator_data::GFLOP;

  using MC_accumulator_data::DCA_iteration;
  using MC_accumulator_data::number_of_measurements;

  using MC_accumulator_data::accumulated_sign;
  using MC_accumulator_data::current_sign;

  const bool compute_std_deviation_;

  std::array<std::vector<vertex_singleton_type>, 2> hs_configuration_;

  std::array<dca::linalg::Matrix<Real, device_t>, 2> M_;
  std::array<dca::linalg::Matrix<Real, linalg::CPU>, 2> M_host_;

  func::function<Real, func::dmn_0<domains::numerical_error_domain>> error;
  func::function<Real, func::dmn_0<Feynman_expansion_order_domain>> visited_expansion_order_k;

  func::function<std::complex<Real>, func::dmn_variadic<nu, nu, r_dmn_t, w>> M_r_w_stddev;

  accumulator::SpAccumulator<Parameters, device_t, Real> single_particle_accumulator_obj;

  std::unique_ptr<ctaux::TpEqualTimeAccumulator<Parameters, Data, Real>> equal_time_accumulator_ptr_;

  accumulator::TpAccumulator<Parameters, device_t> two_particle_accumulator_;

  bool perform_tp_accumulation_ = false;
};

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
CtauxAccumulator<device_t, Parameters, Data, Real>::CtauxAccumulator(Parameters& parameters_ref,
                                                                     Data& data_ref, int id)
    : MC_accumulator_data(),

      parameters_(parameters_ref),
      data_(data_ref),
      concurrency(parameters_.get_concurrency()),

      thread_id(id),

      compute_std_deviation_(parameters_.get_error_computation_type() ==
                             ErrorComputationType::STANDARD_DEVIATION),

      error("numerical-error-distribution-of-N-matrices"),
      visited_expansion_order_k("<k>"),

      M_r_w_stddev("M_r_w_stddev"),

      single_particle_accumulator_obj(parameters_, compute_std_deviation_),

      two_particle_accumulator_(data_.G0_k_w_cluster_excluded, parameters_) {}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::initialize(int dca_iteration) {
  // Note: profiling this function breaks the PAPI profiler as both the master
  // thread and the first
  // worker call this with the same thread_id.
  // TODO: fix thread id assignment.
  //  profiler_type profiler(__FUNCTION__, "CT-AUX accumulator", __LINE__,
  //  thread_id);

  MC_accumulator_data::initialize(dca_iteration);

  if (dca_iteration == parameters_.get_dca_iterations() - 1 && parameters_.isAccumulatingG4())
    perform_tp_accumulation_ = true;

  for (int i = 0; i < visited_expansion_order_k.size(); i++)
    visited_expansion_order_k(i) = 0;

  single_particle_accumulator_obj.resetAccumulation();

  if (perform_tp_accumulation_)
    two_particle_accumulator_.resetAccumulation(dca_iteration);

  if (dca_iteration == parameters_.get_dca_iterations() - 1 &&
      parameters_.additional_time_measurements()) {
    equal_time_accumulator_ptr_ =
        std::make_unique<ctaux::TpEqualTimeAccumulator<Parameters, Data, Real>>(parameters_, data_,
                                                                                thread_id);
    equal_time_accumulator_ptr_->resetAccumulation();
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::finalize() {
  // Note: only one thread calls this function.
  profiler_type profiler(__FUNCTION__, "CT-AUX accumulator", __LINE__);

  single_particle_accumulator_obj.finalize();

  if (compute_std_deviation_) {
    const auto& M_r_w = single_particle_accumulator_obj.get_sign_times_M_r_w();
    const auto& M_r_w_squared = single_particle_accumulator_obj.get_sign_times_M_r_w_sqr();
    for (int l = 0; l < M_r_w_stddev.size(); l++)
      M_r_w_stddev(l) = std::sqrt(abs(M_r_w_squared(l)) - std::pow(abs(M_r_w(l)), 2));

    Real factor = 1. / std::sqrt(parameters_.get_measurements() - 1);

    M_r_w_stddev *= factor;
  }

  if (parameters_.additional_time_measurements())
    equal_time_accumulator_ptr_->finalize();

  if (perform_tp_accumulation_)
    two_particle_accumulator_.finalize();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
std::vector<vertex_singleton>& CtauxAccumulator<device_t, Parameters, Data, Real>::get_configuration(
    e_spin_states_type e_spin) {
  if (e_spin == e_UP)
    return hs_configuration_[0];
  else
    return hs_configuration_[1];
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
template <typename Writer>
void CtauxAccumulator<device_t, Parameters, Data, Real>::write(Writer& writer) {
  //       writer.open_group("CT-AUX-SOLVER-functions");

#ifdef DCA_WITH_QMC_BIT
  writer.execute(error);
#endif  // DCA_WITH_QMC_BIT

  writer.execute(visited_expansion_order_k);

  //       writer.execute(M_r_w);
  //       writer.execute(M_r_w_stddev);

  if (parameters_.additional_time_measurements()) {
    writer.execute(get_charge_cluster_moment());
    writer.execute(get_magnetic_cluster_moment());
    writer.execute(get_dwave_pp_correlator());

    writer.execute(get_G_r_t());
    writer.execute(get_G_r_t_stddev());
  }

  //       writer.close_group();
}

/*!
 *  \brief Get all the information from the walker in order to start a
 * measurement.
 *
 *   \f{eqnarray}{
 *    M_{i,j} &=& (e^{V_i}-1) N_{i,j}
 *   \f}
 */
template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
template <typename walker_type>
void CtauxAccumulator<device_t, Parameters, Data, Real>::updateFrom(walker_type& walker) {
  profiler_type profiler("update from", "CT-AUX accumulator", __LINE__, thread_id);

  GFLOP += walker.get_Gflop();

  current_sign = walker.get_sign();

  const linalg::util::CudaEvent* event = walker.computeM(M_);

  single_particle_accumulator_obj.synchronizeCopy();
  two_particle_accumulator_.synchronizeCopy();

  configuration_type& full_configuration = walker.get_configuration();
  hs_configuration_[0] = full_configuration.get(e_UP);
  hs_configuration_[1] = full_configuration.get(e_DN);

  const int k = full_configuration.get_number_of_interacting_HS_spins();
  if (k < visited_expansion_order_k.size())
    visited_expansion_order_k(k) += 1.;

#ifdef DCA_WITH_QMC_BIT
  error += walker.get_error_distribution();
  walker.get_error_distribution() = 0;
#endif  // DCA_WITH_QMC_BIT

  single_particle_accumulator_obj.syncStreams(*event);
  two_particle_accumulator_.syncStreams(*event);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::measure() {
  number_of_measurements += 1;
  accumulated_sign += current_sign;

  if (perform_tp_accumulation_)
    accumulate_two_particle_quantities();

  accumulate_single_particle_quantities();

  if (DCA_iteration == parameters_.get_dca_iterations() - 1 &&
      parameters_.additional_time_measurements())
    accumulate_equal_time_quantities();
}

#ifdef MEASURE_ERROR_BARS

/*!
 *  \brief Output and store standard deviation and error.
 *
 *  It computes and write to the given files the standard deviation of the
 * measurements of the one
 * particle accumulator.
 *  It outputs the L1-Norm, i.e. \f$\sum_{i=1}^N \left|x_i\right|/N\f$, the
 * L2-Norm, i.e.
 * \f$\sqrt{\sum_{i=1}^N \left|x_i\right|^2/N}\f$,
 *  and the Linf-Norm, i.e. \f$\max_{i=1}^N \left|x_i\right|\f$ of the standard
 * deviation and of the
 * error.
 */
template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::store_standard_deviation(
    int nr_measurements, std::ofstream& points_file, std::ofstream& norm_file) {
  single_particle_accumulator_obj.store_standard_deviation(nr_measurements, points_file, norm_file);
}

/*!
 *  \brief Update the sum of the squares of the measurements of the single
 * particle accumulator.
 *         It has to be called after each measurement.
 */
template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::update_sum_squares() {
  single_particle_accumulator_obj.update_sum_squares();
}
#endif

/*************************************************************
 **                                                         **
 **                    G2 - MEASUREMENTS                    **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::accumulate_single_particle_quantities() {
  profiler_type profiler("sp-accumulation", "CT-AUX accumulator", __LINE__, thread_id);

  single_particle_accumulator_obj.accumulate(M_, hs_configuration_, current_sign);

  GFLOP += 2. * 8. * M_[1].size().first * M_[1].size().first * (1.e-9);
  GFLOP += 2. * 8. * M_[0].size().first * M_[0].size().first * (1.e-9);
}

/*************************************************************
 **                                                         **
 **                 equal-time - MEASUREMENTS               **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::accumulate_equal_time_quantities() {
  profiler_type profiler("equal-time-measurements", "CT-AUX accumulator", __LINE__, thread_id);

  return accumulate_equal_time_quantities(M_);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::accumulate_equal_time_quantities(
    const std::array<linalg::Matrix<Real, linalg::GPU>, 2>& M) {
  for (int s = 0; s < 2; ++s)
    M_host_[s].setAsync(M[s], thread_id, s);
  for (int s = 0; s < 2; ++s)
    linalg::util::syncStream(thread_id, s);

  return accumulate_equal_time_quantities(M_host_);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::accumulate_equal_time_quantities(
    const std::array<linalg::Matrix<Real, linalg::CPU>, 2>& M) {
  equal_time_accumulator_ptr_->accumulateAll(hs_configuration_[0], M[0], hs_configuration_[1], M[1],
                                             current_sign);

  GFLOP += equal_time_accumulator_ptr_->get_GFLOP();
}

/*************************************************************
 **                                                         **
 **                 nonlocal \chi - MEASUREMENTS            **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::accumulate_two_particle_quantities() {
  profiler_type profiler("tp-accumulation", "CT-AUX accumulator", __LINE__, thread_id);
  GFLOP += 1e-9 * two_particle_accumulator_.accumulate(M_, hs_configuration_, current_sign);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, typename Real>
void CtauxAccumulator<device_t, Parameters, Data, Real>::sumTo(this_type& other) {
  other.GFLOP += GFLOP;

  other.accumulated_sign += accumulated_sign;
  other.number_of_measurements += number_of_measurements;

  other.get_visited_expansion_order_k() += visited_expansion_order_k;
  other.get_error_distribution() += error;

  // sp-measurements
  single_particle_accumulator_obj.sumTo(other.single_particle_accumulator_obj);

  // equal time measurements
  if (DCA_iteration == parameters_.get_dca_iterations() - 1 &&
      parameters_.additional_time_measurements())
    equal_time_accumulator_ptr_->sumTo(*other.equal_time_accumulator_ptr_);

  // tp-measurements
  if (perform_tp_accumulation_)
    two_particle_accumulator_.sumTo(other.two_particle_accumulator_);
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_ACCUMULATOR_HPP
