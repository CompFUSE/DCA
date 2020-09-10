// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
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

#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/cuda_stream.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <class Parameters, linalg::DeviceType device, typename Real>
class CtintAccumulator {
public:
  using this_type = CtintAccumulator<Parameters, device, Real>;
  using ParametersType = Parameters;
  using DataType = phys::DcaData<Parameters>;

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

  const auto& get_sign_times_M_r_w() const;

  const auto& get_sign_times_G4() const;

  double avgSign() const {
    return accumulated_sign_.mean();
  }

  int get_accumulated_sign() const {
    return accumulated_sign_.sum();
  }
  int get_number_of_measurements() const {
    return accumulated_sign_.count();
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
    return accumulator::TpAccumulator<Parameters, device>::staticDeviceFingerprint();
  }

  float getFLOPs() const {
    return flop_;
  }

private:
  const Parameters& parameters_;

  // Internal instantaneous configuration.
  std::array<linalg::Matrix<Real, device>, 2> M_;
  MatrixConfiguration configuration_;
  int sign_ = 0;

  std::vector<const linalg::util::CudaStream*> streams_;
  linalg::util::CudaEvent event_;

  util::Accumulator<int> accumulated_sign_;
  util::Accumulator<unsigned long> accumulated_order_;

  const int thread_id_;

  accumulator::SpAccumulator<Parameters, device, Real> sp_accumulator_;
  accumulator::TpAccumulator<Parameters, device> tp_accumulator_;

  bool perform_tp_accumulation_ = false;
  bool ready_ = false;
  bool finalized_ = false;

  float flop_ = 0.;
};

template <class Parameters, linalg::DeviceType device, typename Real>
template <class Data>
CtintAccumulator<Parameters, device, Real>::CtintAccumulator(const Parameters& pars, const Data& data,
                                                       int id)
    : parameters_(pars),
      thread_id_(id),
      sp_accumulator_(pars),
      tp_accumulator_(data.G0_k_w_cluster_excluded, pars, thread_id_) {
  streams_.insert(streams_.end(), sp_accumulator_.get_streams().begin(),
                  sp_accumulator_.get_streams().end());
  streams_.push_back(tp_accumulator_.get_stream());
}

template <class Parameters, linalg::DeviceType device, typename Real>
void CtintAccumulator<Parameters, device, Real>::initialize(const int dca_iteration) {
  perform_tp_accumulation_ =
      parameters_.isAccumulatingG4() && dca_iteration == parameters_.get_dca_iterations() - 1;
  accumulated_order_.reset();
  accumulated_sign_.reset();

  sp_accumulator_.resetAccumulation();
  if (perform_tp_accumulation_)
    tp_accumulator_.resetAccumulation(dca_iteration);

  finalized_ = false;
}

template <class Parameters, linalg::DeviceType device, typename Real>
template <class Walker>
void CtintAccumulator<Parameters, device, Real>::updateFrom(Walker& walker) {
  // Compute M.
  walker.computeM(M_);

  if (device == linalg::GPU) {
    for (int s = 0; s < 2; ++s) {
      event_.record(walker.get_stream(s));
      //  Synchronize sp accumulator streams with walker.
      event_.block(*sp_accumulator_.get_streams()[s]);
      //  Synchronize both walker streams with tp accumulator.
      event_.block(*tp_accumulator_.get_stream());
    }
  }

  configuration_ = walker.getConfiguration();
  sign_ = walker.get_sign();
  flop_ += walker.stealFLOPs();

  ready_ = true;
}

template <class Parameters, linalg::DeviceType device, typename Real>
template <class Walker>
void CtintAccumulator<Parameters, device, Real>::accumulate(Walker& walker) {
  updateFrom(walker);
  measure();
}

template <class Parameters, linalg::DeviceType device, typename Real>
void CtintAccumulator<Parameters, device, Real>::measure() {
  if (!ready_ || sign_ == 0)
    throw(std::logic_error("No or invalid configuration to accumulate."));

  accumulated_sign_.addSample(sign_);
  accumulated_order_.addSample(order());

  sp_accumulator_.accumulate(M_, configuration_.get_sectors(), sign_);
  if (perform_tp_accumulation_)
    tp_accumulator_.accumulate(M_, configuration_.get_sectors(), sign_);

  ready_ = false;
}

template <class Parameters, linalg::DeviceType device, typename Real>
void CtintAccumulator<Parameters, device, Real>::sumTo(this_type& other_one) {
  other_one.accumulated_order_ += accumulated_order_;
  other_one.accumulated_sign_ += accumulated_sign_;

  sp_accumulator_.sumTo(other_one.sp_accumulator_);
  if (perform_tp_accumulation_) {
    assert(other_one.perform_tp_accumulation_);
    tp_accumulator_.sumTo(other_one.tp_accumulator_);
  }
  other_one.flop_ += flop_;
}

template <class Parameters, linalg::DeviceType device, typename Real>
void CtintAccumulator<Parameters, device, Real>::finalize() {
  if (finalized_)
    return;

  sp_accumulator_.finalize();
  if (perform_tp_accumulation_)
    tp_accumulator_.finalize();

  finalized_ = true;
}

template <class Parameters, linalg::DeviceType device, typename Real>
const auto& CtintAccumulator<Parameters, device, Real>::get_sign_times_M_r_w() const {
  return sp_accumulator_.get_sign_times_M_r_w();
}

template <class Parameters, linalg::DeviceType device, typename Real>
const auto& CtintAccumulator<Parameters, device, Real>::get_sign_times_G4() const {
  if (!perform_tp_accumulation_)
    throw(std::logic_error("G4 was not accumulated."));
  return tp_accumulator_.get_sign_times_G4();
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_CTINT_ACCUMULATOR_HPP
