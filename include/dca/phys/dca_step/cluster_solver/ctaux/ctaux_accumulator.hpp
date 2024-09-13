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
#include "dca/linalg/util/gpu_event.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_id.hpp"
#ifdef DCA_HAVE_GPU
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"
#else
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator.hpp"
#endif  // DCA_HAVE_GPU
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_cpu.hpp"

#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/tp_equal_time_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/feynman_expansion_order_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ct_aux_hs_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_pair.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/mc_accumulator_data.hpp"

#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/quantum/numerical_error_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/four_point_type.hpp"
#ifdef DCA_HAVE_MPI
// #include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_mpi_blocked_gpu.hpp"
// #include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_mpi_gpu.hpp"
#endif  // DCA_HAVE_MPI

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
class CtauxAccumulator : public MC_accumulator_data<typename Parameters::Scalar> {
public:
  static constexpr ClusterSolverId solver_id{ClusterSolverId::CT_AUX};
  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;
  using this_type = CtauxAccumulator<device_t, Parameters, Data, DIST>;
  using TpAccumulator = accumulator::TpAccumulator<Parameters, DIST, device_t>;
  using ParametersType = Parameters;
  using DataType = Data;
  using BaseClass = MC_accumulator_data<Scalar>;
  
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
  using SpAccumulator = typename accumulator::SpAccumulator<Parameters, device_t>;
  using MFunction = typename SpAccumulator::MFunction;
  using MFunctionTime = typename SpAccumulator::MFunctionTime;
  using MFunctionTimePair = typename SpAccumulator::MFunctionTimePair;
  using FTauPair = typename SpAccumulator::FTauPair;
  using PaddedTimeDmn = typename SpAccumulator::PaddedTimeDmn;

  typedef CT_AUX_HS_configuration<Parameters> configuration_type;

  CtauxAccumulator(const Parameters& parameters_ref, Data& data_ref, int id);

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

  auto& get_visited_expansion_order_k() {
    return visited_expansion_order_k;
  }

  // equal time-measurements
  // TODO: Make equal time getters const.
  auto& get_G_r_t() {
    return equal_time_accumulator_ptr_->get_G_r_t();
  }
  func::function<Scalar, func::dmn_variadic<nu, nu, r_dmn_t, t>>& get_G_r_t_stddev() {
    return equal_time_accumulator_ptr_->get_G_r_t_stddev();
  }

  auto& get_charge_cluster_moment() {
    return equal_time_accumulator_ptr_->get_charge_cluster_moment();
  }
  auto& get_magnetic_cluster_moment() {
    return equal_time_accumulator_ptr_->get_magnetic_cluster_moment();
  }
  auto& get_dwave_pp_correlator() {
    return equal_time_accumulator_ptr_->get_dwave_pp_correlator();
  }

  // sp-measurements
  const auto& get_sign() const {
    return current_phase_;
  }

  const auto& get_sign_times_M_r_w() const {
    return single_particle_accumulator_obj.get_sign_times_M_r_w();
  }

  const MFunction& get_single_measurement_sign_times_MFunction() {
    return single_particle_accumulator_obj.get_single_measurement_sign_times_MFunction();
  }

  const FTauPair& get_single_measurement_sign_times_MFunction_time() {
    return single_particle_accumulator_obj.get_single_measurement_sign_times_MFunction_time();
  }

  void clearSingleMeasurement();

  const auto& get_sign_times_M_r_w_sqr() const {
    return single_particle_accumulator_obj.get_sign_times_M_r_w_sqr();
  }

  // tp-measurements
  const auto& get_sign_times_G4() {
    return two_particle_accumulator_.get_G4();
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
    return accumulator::TpAccumulator<Parameters, DIST, device_t>::staticDeviceFingerprint();
  }

  bool perform_tp_accumulation() const {
    return perform_tp_accumulation_;
  }
  bool perform_equal_time_accumulation() const {
    return perform_equal_time_accumulation_;
  }

private:
  void accumulate_single_particle_quantities();

  void accumulate_equal_time_quantities();
  void accumulate_equal_time_quantities(const std::array<linalg::Matrix<Scalar, linalg::GPU>, 2>& M);
  void accumulate_equal_time_quantities(const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M);

  void accumulate_two_particle_quantities();

protected:
  const Parameters& parameters_;
  Data& data_;
  const concurrency_type& concurrency;

  int thread_id;

  using MC_accumulator_data<Scalar>::gflop_;

  using MC_accumulator_data<Scalar>::dca_iteration_;
  using MC_accumulator_data<Scalar>::number_of_measurements_;

  using MC_accumulator_data<Scalar>::accumulated_phase_;
  using MC_accumulator_data<Scalar>::current_phase_;

  const bool compute_std_deviation_;

  std::array<std::vector<vertex_singleton_type>, 2> hs_configuration_;

  std::array<dca::linalg::Matrix<Scalar, device_t>, 2> M_;
  std::array<dca::linalg::Matrix<Scalar, linalg::CPU>, 2> M_host_;

  func::function<Real, func::dmn_0<domains::numerical_error_domain>> error;
  func::function<Real, func::dmn_0<Feynman_expansion_order_domain>> visited_expansion_order_k;

  func::function<std::complex<Real>, func::dmn_variadic<nu, nu, r_dmn_t, w>> M_r_w_stddev;

  accumulator::SpAccumulator<Parameters, device_t> single_particle_accumulator_obj;

  accumulator::TpAccumulator<Parameters, DIST, device_t> two_particle_accumulator_;

  std::unique_ptr<ctaux::TpEqualTimeAccumulator<Parameters, Data>> equal_time_accumulator_ptr_;

  bool perform_tp_accumulation_ = false;
  bool perform_equal_time_accumulation_ = false;
};

/** This constructor takes a number of references for later convenience.
 *
 *  There is good possibility of initialization order issues here.
 *  So its better the use the passed in parameters_ref and data_ref than their local references
 *  in the initializers.
 */
template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
CtauxAccumulator<device_t, Parameters, Data, DIST>::CtauxAccumulator(
    const Parameters& parameters_ref, Data& data_ref, int id)
    : MC_accumulator_data<Scalar>(),

      parameters_(parameters_ref),
      data_(data_ref),
      concurrency(parameters_ref.get_concurrency()),

      thread_id(id),

      compute_std_deviation_(parameters_ref.get_error_computation_type() ==
                             ErrorComputationType::STANDARD_DEVIATION),

      error("numerical-error-distribution-of-N-matrices"),
      visited_expansion_order_k("<k>"),

      M_r_w_stddev("M_r_w_stddev"),

      single_particle_accumulator_obj(parameters_ref, compute_std_deviation_),

      two_particle_accumulator_(data_ref.G0_k_w_cluster_excluded, parameters_ref, id) {}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::initialize(int dca_iteration) {
  // Note: profiling this function breaks the PAPI profiler as both the master
  // thread and the first
  // worker call this with the same thread_id.
  // TODO: fix thread id assignment.
  //  profiler_type profiler(__FUNCTION__, "CT-AUX accumulator", __LINE__,
  //  thread_id);

  dca_iteration_ = dca_iteration;
  MC_accumulator_data<Scalar>::initialize(dca_iteration);

  perform_tp_accumulation_ =
      parameters_.isAccumulatingG4() && ((dca_iteration_ == parameters_.get_dca_iterations() - 1) ||
                                         parameters_.dump_every_iteration());

  for (int i = 0; i < visited_expansion_order_k.size(); i++)
    visited_expansion_order_k(i) = 0;

  single_particle_accumulator_obj.resetAccumulation();
  single_particle_accumulator_obj.clearSingleMeasurement();

  if (perform_tp_accumulation_)
    two_particle_accumulator_.resetAccumulation(dca_iteration_);

  perform_equal_time_accumulation_ = parameters_.additional_time_measurements() &&
                                     ((dca_iteration_ == parameters_.get_dca_iterations() - 1) ||
                                      parameters_.dump_every_iteration());

  if (perform_equal_time_accumulation_) {
    equal_time_accumulator_ptr_ =
        std::make_unique<ctaux::TpEqualTimeAccumulator<Parameters, Data>>(parameters_, data_,
                                                                                thread_id);
    equal_time_accumulator_ptr_->resetAccumulation();
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::finalize() {
  // Note: only one thread calls this function.
  profiler_type profiler(__FUNCTION__, "CT-AUX accumulator", __LINE__);

  single_particle_accumulator_obj.finalize();

  if (compute_std_deviation_) {
    const auto& M_r_w = single_particle_accumulator_obj.get_sign_times_M_r_w();
    const auto& M_r_w_squared = single_particle_accumulator_obj.get_sign_times_M_r_w_sqr();
    for (int l = 0; l < M_r_w_stddev.size(); l++)
      M_r_w_stddev(l) = std::sqrt(abs(M_r_w_squared(l)) - std::pow(abs(M_r_w(l)), 2));

    Real factor = 1. / std::sqrt(parameters_.get_measurements().at(dca_iteration_) - 1);

    M_r_w_stddev *= factor;
  }

  if (perform_equal_time_accumulation_)
    equal_time_accumulator_ptr_->finalize();

  if (perform_tp_accumulation_)
    two_particle_accumulator_.finalize();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
std::vector<vertex_singleton>& CtauxAccumulator<device_t, Parameters, Data, DIST>::get_configuration(
    e_spin_states_type e_spin) {
  if (e_spin == e_UP)
    return hs_configuration_[0];
  else
    return hs_configuration_[1];
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
template <typename Writer>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::write(Writer& writer) {
  // it is assumed this this is called from CtauxClusterSolver and so we are in a
  // CT-AUX-SOLVER-functions group.

#ifdef DCA_WITH_QMC_BIT
  writer.execute(error);
#endif  // DCA_WITH_QMC_BIT

  writer.execute(visited_expansion_order_k);

  // equal time should just have a write
  if (perform_equal_time_accumulation_) {
    writer.execute(get_charge_cluster_moment());
    writer.execute(get_magnetic_cluster_moment());
    writer.execute(get_dwave_pp_correlator());

    writer.execute(get_G_r_t());
    writer.execute(get_G_r_t_stddev());
  }
}

/*!
 *  \brief Get all the information from the walker in order to start a
 * measurement.
 *
 *   \f{eqnarray}{
 *    M_{i,j} &=& (e^{V_i}-1) N_{i,j}
 *   \f}
 */
template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
template <typename walker_type>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::updateFrom(walker_type& walker) {
  profiler_type profiler("update from", "CT-AUX accumulator", __LINE__, thread_id);

  gflop_ += walker.get_Gflop();

  current_phase_ = walker.get_sign();

  const linalg::util::GpuEvent* event = walker.computeM(M_);

  single_particle_accumulator_obj.synchronizeCopy();
  two_particle_accumulator_.synchronizeCopy();

  configuration_type& full_configuration = walker.get_configuration();
  hs_configuration_[0] = full_configuration.get(e_DN);
  hs_configuration_[1] = full_configuration.get(e_UP);

  const int k = full_configuration.get_number_of_interacting_HS_spins();
  if (k < visited_expansion_order_k.size())
    visited_expansion_order_k(k) += 1.;

#ifdef DCA_WITH_QMC_BIT
  error += walker.get_error_distribution();
  walker.get_error_distribution() = 0;
#endif  // DCA_WITH_QMC_BIT

  //single_particle_accumulator_obj.syncStreams(*event);
  //two_particle_accumulator_.syncStreams(*event);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::measure() {
  number_of_measurements_ += 1;
  accumulated_phase_.addSample(current_phase_.getSign());

  if (perform_tp_accumulation_)
    accumulate_two_particle_quantities();

  accumulate_single_particle_quantities();

  if (perform_equal_time_accumulation_)
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
template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::store_standard_deviation(
    int nr_measurements, std::ofstream& points_file, std::ofstream& norm_file) {
  single_particle_accumulator_obj.store_standard_deviation(nr_measurements, points_file, norm_file);
}

/*!
 *  \brief Update the sum of the squares of the measurements of the single
 * particle accumulator.
 *         It has to be called after each measurement.
 */
template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::update_sum_squares() {
  single_particle_accumulator_obj.update_sum_squares();
}
#endif

/*************************************************************
 **                                                         **
 **                    G2 - MEASUREMENTS                    **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::accumulate_single_particle_quantities() {
  profiler_type profiler("sp-accumulation", "CT-AUX accumulator", __LINE__, thread_id);

  single_particle_accumulator_obj.accumulate(M_, hs_configuration_, current_phase_.getSign());

  gflop_ += 2. * 8. * M_[1].size().first * M_[1].size().first * (1.e-9);
  gflop_ += 2. * 8. * M_[0].size().first * M_[0].size().first * (1.e-9);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::clearSingleMeasurement() {
  single_particle_accumulator_obj.clearSingleMeasurement();
}

/*************************************************************
 **                                                         **
 **                 equal-time - MEASUREMENTS               **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::accumulate_equal_time_quantities() {
  profiler_type profiler("equal-time-measurements", "CT-AUX accumulator", __LINE__, thread_id);

  return accumulate_equal_time_quantities(M_);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::accumulate_equal_time_quantities(
    const std::array<linalg::Matrix<Scalar, linalg::GPU>, 2>& M) {
  for (int s = 0; s < 2; ++s)
    M_host_[s].setAsync(M[s], thread_id, s);
  for (int s = 0; s < 2; ++s)
    linalg::util::syncStream(thread_id, s);

  return accumulate_equal_time_quantities(M_host_);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::accumulate_equal_time_quantities(
    const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M) {
  equal_time_accumulator_ptr_->accumulateAll(hs_configuration_[0], M[0], hs_configuration_[1], M[1],
                                             current_phase_.getSign());

  gflop_ += equal_time_accumulator_ptr_->get_gflop();
}

/*************************************************************
 **                                                         **
 **                 nonlocal \chi - MEASUREMENTS            **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::accumulate_two_particle_quantities() {
  profiler_type profiler("tp-accumulation", "CT-AUX accumulator", __LINE__, thread_id);
  gflop_ += 1e-9 * two_particle_accumulator_.accumulate(M_, hs_configuration_, current_phase_.getSign());
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data, DistType DIST>
void CtauxAccumulator<device_t, Parameters, Data, DIST>::sumTo(this_type& other) {
  other.gflop_ += gflop_;

  other.accumulated_phase_ += accumulated_phase_;
  other.number_of_measurements_ += number_of_measurements_;

  other.get_visited_expansion_order_k() += visited_expansion_order_k;
  other.get_error_distribution() += error;

  // sp-measurements
  single_particle_accumulator_obj.sumTo(other.single_particle_accumulator_obj);

  // equal time measurements
  if (perform_equal_time_accumulation_)
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
