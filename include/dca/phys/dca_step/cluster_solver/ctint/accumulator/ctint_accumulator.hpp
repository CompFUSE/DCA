// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_CTINT_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_CTINT_ACCUMULATOR_HPP

#include <cassert>
#include <memory>
#include <stdexcept>

#include "dca/linalg/util/cuda_stream.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <class Parameters, linalg::DeviceType device>
class CtintAccumulator {
public:
  using this_type = CtintAccumulator<Parameters, device>;
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
    return double(total_sign_) / double(total_meas_);
  }

  int get_accumulated_sign() const {
    return total_sign_;
  }
  int get_number_of_measurements() const {
    return total_meas_;
  }

  int order() const {
    return (configuration_.size(0) + configuration_.size(1)) / 2;
  }
  double avgOrder() const {
    return (double)order_sum_ / (double)total_meas_;
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
  std::array<linalg::Matrix<double, device>, 2> M_;
  MatrixConfiguration configuration_;
  int sign_ = 0;

  std::vector<linalg::util::CudaStream*> streams_;

  int total_sign_ = 0;
  uint total_meas_ = 0;
  ulong order_sum_ = 0;

  const int thread_id_;

  accumulator::SpAccumulator<Parameters, device> sp_accumulator_;
  accumulator::TpAccumulator<Parameters, device> tp_accumulator_;

  bool perform_tp_accumulation_ = false;
  bool ready_ = false;
  bool finalized_ = false;

  float flop_ = 0.;
};

template <class Parameters, linalg::DeviceType device>
template <class Data>
CtintAccumulator<Parameters, device>::CtintAccumulator(const Parameters& pars, const Data& data,
                                                       int id)
    : parameters_(pars),
      thread_id_(id),
      sp_accumulator_(pars),
      tp_accumulator_(data.G0_k_w_cluster_excluded, pars, thread_id_) {
  streams_.insert(streams_.end(), sp_accumulator_.get_streams().begin(),
                  sp_accumulator_.get_streams().end());
  streams_.push_back(tp_accumulator_.get_stream());
}

template <class Parameters, linalg::DeviceType device>
void CtintAccumulator<Parameters, device>::initialize(const int dca_iteration) {
  perform_tp_accumulation_ =
      parameters_.accumulateG4() && dca_iteration == parameters_.get_dca_iterations() - 1;
  total_sign_ = 0;
  total_meas_ = 0;

  sp_accumulator_.resetAccumulation();
  if (perform_tp_accumulation_)
    tp_accumulator_.resetAccumulation(dca_iteration);

  finalized_ = false;
}

template <class Parameters, linalg::DeviceType device>
template <class Walker>
void CtintAccumulator<Parameters, device>::updateFrom(Walker& walker) {
  walker.computeM(M_, streams_);
  configuration_ = walker.getConfiguration();
  sign_ = walker.getSign();
  flop_ += walker.stealFLOPs();

  ready_ = true;
}

template <class Parameters, linalg::DeviceType device>
template <class Walker>
void CtintAccumulator<Parameters, device>::accumulate(Walker& walker) {
  updateFrom(walker);
  measure();
}

template <class Parameters, linalg::DeviceType device>
void CtintAccumulator<Parameters, device>::measure() {
  if (!ready_ || sign_ == 0)
    throw(std::logic_error("No or invalid configuration to accumulate."));

  total_sign_ += sign_;
  ++total_meas_;
  order_sum_ += order();
  sp_accumulator_.accumulate(M_, configuration_.get_sectors(), sign_);
  if (perform_tp_accumulation_)
    tp_accumulator_.accumulate(M_, configuration_.get_sectors(), sign_);

  ready_ = false;
}

template <class Parameters, linalg::DeviceType device>
void CtintAccumulator<Parameters, device>::sumTo(this_type& other_one) {
  other_one.total_meas_ += total_meas_;
  other_one.total_sign_ += total_sign_;
  sp_accumulator_.sumTo(other_one.sp_accumulator_);
  if (perform_tp_accumulation_) {
    assert(other_one.perform_tp_accumulation_);
    tp_accumulator_.sumTo(other_one.tp_accumulator_);
  }
  other_one.flop_ += flop_;
  other_one.order_sum_ += order_sum_;
}

template <class Parameters, linalg::DeviceType device>
void CtintAccumulator<Parameters, device>::finalize() {
  if (finalized_)
    return;

  sp_accumulator_.finalize();
  if (perform_tp_accumulation_)
    tp_accumulator_.finalize();

  finalized_ = true;
}

template <class Parameters, linalg::DeviceType device>
const auto& CtintAccumulator<Parameters, device>::get_sign_times_M_r_w() const {
  return sp_accumulator_.get_sign_times_M_r_w();
}

template <class Parameters, linalg::DeviceType device>
const auto& CtintAccumulator<Parameters, device>::get_sign_times_G4() const {
  if (!perform_tp_accumulation_)
    throw(std::logic_error("G4 was not accumulated."));
  return tp_accumulator_.get_sign_times_G4();
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_CTINT_ACCUMULATOR_HPP
