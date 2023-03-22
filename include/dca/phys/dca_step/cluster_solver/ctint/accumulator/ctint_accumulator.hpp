// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Accumulator for the CT-INT solver.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_CTINT_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_CTINT_ACCUMULATOR_HPP

#include <cassert>
#include <memory>
#include <stdexcept>

#include "dca/linalg/util/gpu_event.hpp"
#include "dca/linalg/util/gpu_stream.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_cpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/mc_accumulator_data.hpp"
#include "dca/util/type_utils.hpp"
#ifdef DCA_HAVE_GPU
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"
#endif  // DCA_HAVE_GPU

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <class Parameters, linalg::DeviceType device,
          DistType DIST = dca::DistType::NONE>
class CtintAccumulator : public MC_accumulator_data<typename dca::util::ScalarSelect<typename dca::config::McOptions::MC_REAL,Parameters::complex_g0>::type> {
public:
  constexpr static ClusterSolverId solver_id{ClusterSolverId::CT_INT};
  using Real = typename dca::config::McOptions::MC_REAL;
  using Scalar = typename dca::util::ScalarSelect<Real,Parameters::complex_g0>::type;
  using Base = MC_accumulator_data<Scalar>;
  using this_type = CtintAccumulator<Parameters, device, DIST>;
  using Base::accumulated_phase_;
  using Base::current_phase_;
  using Base::number_of_measurements_;
  
  using ParametersType = Parameters;
  using DataType = phys::DcaData<Parameters, DIST>;
  using SpAccumulator = accumulator::SpAccumulator<Parameters, device>;
  // Same form as SpGreensFunction
  using MFunction = typename SpAccumulator::MFunction;
  using MFunctionTime = typename SpAccumulator::MFunctionTime;
  using MFunctionTimePair = typename SpAccumulator::MFunctionTimePair;
  using FTauPair = typename SpAccumulator::FTauPair;
  using PaddedTimeDmn = typename SpAccumulator::PaddedTimeDmn;

  template <class Data>
  CtintAccumulator(const Parameters& pars, const Data& data, int id = 0);
  CtintAccumulator(this_type&& other_one) = default;

  void initialize(int dca_iteration_);

  void measure();

  template <class Walker>
  void updateFrom(Walker& walker);
  template <class Walker>
  void accumulate(Walker& walker);

  void finalize();

  void sumTo(this_type& other_acc);

  /** write runtime parameters used by ctint_accumulator and its important owned objects */
  template <class Writer>
  void write(Writer& writer) {
    writer.open_group("accumulator");
    sp_accumulator_.write(writer);
    writer.close_group();
  }

  const auto& get_sign() const {
    return current_phase_;
  }

  const auto& get_sign_times_M_r_w() const;

  const MFunction& get_single_measurement_sign_times_MFunction() {
    return sp_accumulator_.get_single_measurement_sign_times_MFunction();
  }

  const FTauPair& get_single_measurement_sign_times_MFunction_time() {
    return sp_accumulator_.get_single_measurement_sign_times_MFunction_time();
  }

  void clearSingleMeasurement();

  const auto& get_sign_times_G4() const;

  double avgSign() const {
    return accumulated_phase_.mean();
  }

  auto get_accumulated_phase() const {
    return accumulated_phase_.sum();
  }

  int get_number_of_measurements() const {
    assert(accumulated_phase_.count() == number_of_measurements_);
    return number_of_measurements_;
  }

  int order() const {
    return (configuration_.size(0) + configuration_.size(1)) / 2;
  }
  double avgOrder() const {
    return accumulated_order_.mean();
  }

  std::size_t deviceFingerprint() const {
    return sp_accumulator_.deviceFingerprint() + tp_accumulator_.deviceFingerprint();
  }

  static std::size_t staticDeviceFingerprint() {
    return accumulator::TpAccumulator<Parameters, DIST, device>::staticDeviceFingerprint();
  }

  bool perform_tp_accumulation() const {
    return perform_tp_accumulation_;
  }

  double getFLOPs() const {
    return flop_;
  }

private:
  const Parameters& parameters_;

  // Internal instantaneous configuration.
  std::array<linalg::Matrix<Scalar, device>, 2> M_;
  MatrixConfiguration configuration_;

  std::vector<const linalg::util::GpuStream*> streams_;
  linalg::util::GpuEvent event_;

  util::Accumulator<unsigned long> accumulated_order_;

  const int thread_id_;

  SpAccumulator sp_accumulator_;
  accumulator::TpAccumulator<Parameters, DIST, device> tp_accumulator_;

  bool perform_tp_accumulation_ = false;
  bool ready_ = false;
  bool finalized_ = false;

  double flop_ = 0.;

  // cache flop calculations to avoid doing them in hot functions
  // Flops for a measurement
  double measure_flops_ = 0.0;
};

template <class Parameters, linalg::DeviceType device, DistType DIST>
template <class Data>
CtintAccumulator<Parameters, device, DIST>::CtintAccumulator(const Parameters& pars,
                                                                   const Data& data, int id)
    : parameters_(pars),
      thread_id_(id),
      sp_accumulator_(pars),
      tp_accumulator_(data.G0_k_w_cluster_excluded, pars, thread_id_) {
  streams_.insert(streams_.end(), sp_accumulator_.get_streams().begin(),
                  sp_accumulator_.get_streams().end());
  streams_.push_back(tp_accumulator_.get_stream());
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
void CtintAccumulator<Parameters, device, DIST>::initialize(const int dca_iteration) {
  perform_tp_accumulation_ =
      parameters_.isAccumulatingG4() && ((dca_iteration == parameters_.get_dca_iterations() - 1) ||
                                         parameters_.dump_every_iteration());
  accumulated_order_.reset();
  
  Base::initialize(dca_iteration);
  sp_accumulator_.resetAccumulation();
  sp_accumulator_.clearSingleMeasurement();
  if (perform_tp_accumulation_)
    tp_accumulator_.resetAccumulation(dca_iteration);

  flop_ = 0.0;
  measure_flops_ = M_[0].nrCols() * M_[0].nrCols() * 2 * 2 * 8 * 19;

  finalized_ = false;
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
template <class Walker>
void CtintAccumulator<Parameters, device, DIST>::updateFrom(Walker& walker) {
  // Compute M.
  auto m_size = M_[0].nrCols();
  walker.computeM(M_);

  // if accumulating the M_ squared too + another term with 22 instead of 19
  // first factor of 2 is because there are 2 M_
  // second 2 is a factor from the device accumulation
  // 8 is oversampling
  if (m_size != M_[0].nrCols())
    measure_flops_ = M_[0].nrCols() * M_[0].nrCols() * 2 * 2 * 8 * 19;

  if constexpr (device == linalg::GPU) {
    for (int s = 0; s < 2; ++s) {
      event_.record(walker.get_stream(s));
      //  Synchronize sp accumulator streams with walker.
      event_.block(*sp_accumulator_.get_streams()[s]);
      //  Synchronize both walker streams with tp accumulator.
      event_.block(*tp_accumulator_.get_stream());
    }
  }

  configuration_ = walker.getConfiguration();
  current_phase_ = walker.get_sign();
  flop_ += walker.stealFLOPs();

  ready_ = true;
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
template <class Walker>
void CtintAccumulator<Parameters, device, DIST>::accumulate(Walker& walker) {
  updateFrom(walker);
  measure();
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
void CtintAccumulator<Parameters, device, DIST>::measure() {
  if (!ready_ || current_phase_.isNull())
    throw(std::logic_error("No or invalid configuration to accumulate."));
  accumulated_phase_.addSample(current_phase_.getSign());
  accumulated_order_.addSample(order());
  number_of_measurements_ += 1;

  sp_accumulator_.accumulate(M_, configuration_.get_sectors(), current_phase_.getSign());
  flop_ += measure_flops_;

  if (perform_tp_accumulation_)
    tp_accumulator_.accumulate(M_, configuration_.get_sectors(), current_phase_.getSign());

  ready_ = false;
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
void CtintAccumulator<Parameters, device, DIST>::sumTo(this_type& other_one) {
  other_one.accumulated_order_ += accumulated_order_;
  other_one.accumulated_phase_ += accumulated_phase_;

  sp_accumulator_.sumTo(other_one.sp_accumulator_);
  if (perform_tp_accumulation_) {
    assert(other_one.perform_tp_accumulation_);
    tp_accumulator_.sumTo(other_one.tp_accumulator_);
  }
  other_one.flop_ += flop_;
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
void CtintAccumulator<Parameters, device, DIST>::finalize() {
  if (finalized_)
    return;

  sp_accumulator_.finalize();
  if (perform_tp_accumulation_)
    tp_accumulator_.finalize();

  finalized_ = true;
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
const auto& CtintAccumulator<Parameters, device, DIST>::get_sign_times_M_r_w() const {
  return sp_accumulator_.get_sign_times_M_r_w();
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
void CtintAccumulator<Parameters, device, DIST>::clearSingleMeasurement() {
  sp_accumulator_.clearSingleMeasurement();
}

template <class Parameters, linalg::DeviceType device, DistType DIST>
const auto& CtintAccumulator<Parameters, device, DIST>::get_sign_times_G4() const {
  if (!perform_tp_accumulation_)
    throw(std::logic_error("G4 was not accumulated."));
  return tp_accumulator_.get_G4();
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_CTINT_ACCUMULATOR_HPP
