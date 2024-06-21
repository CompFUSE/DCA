// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// A std::thread jacket that measures correlations in the MC walker independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_WALKER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_WALKER_HPP

#include <mutex>

#include "dca/math/statistics/autocorrelation.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/time_correlator.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_walker.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/qmci_autocorrelation_data.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {
// dca::phys::solver::stdthreadqmci::

template <class QmciWalker, class DATA>
class StdThreadQmciWalker final
    : public QmciWalker {  // , public QmciAutocorrelationData<QmciWalker> {
  using ThisType = StdThreadQmciWalker<QmciWalker, DATA>;
  using Parameters = typename QmciWalker::parameters_type;
  using Data = DATA;
  using Concurrency = typename Parameters::concurrency_type;
  using Rng = typename Parameters::random_number_generator;
  using Real = typename QmciWalker::Scalar;

  constexpr static auto device = QmciWalker::device;
  constexpr static int bands = Parameters::bands;

public:
  StdThreadQmciWalker(Parameters& parameters_ref, DATA& data_ref, Rng& rng, int concurrency_id,
                      int id, const std::shared_ptr<io::Writer<Concurrency>>& writer,
                      G0Interpolation<device, Real>& g0);

  void initialize(int iteration_);
  void doSweep();

  void markThermalized() {
    //    QmciAutocorrelationData<QmciWalker>::markThermalized();
    QmciWalker::markThermalized();
  }

  int get_thread_id() const {
    return thread_id_;
  }
  bool get_last_iteration() const {
    return last_iteration_;
  }

  std::size_t get_meas_id() const {
    return meas_id_;
  }

private:
  void logConfiguration() const;

  const unsigned stamping_period_;
  int concurrency_id_;
  int thread_id_;

  std::shared_ptr<io::Writer<Concurrency>> writer_;
  std::size_t meas_id_ = 0;

  const int total_iterations_;

  bool last_iteration_ = false;
};  // namespace stdthreadqmci

template <class QmciWalker, class DATA>
StdThreadQmciWalker<QmciWalker, DATA>::StdThreadQmciWalker(
    Parameters& parameters, DATA& data_ref, Rng& rng, int concurrency_id, int id,
    const std::shared_ptr<io::Writer<Concurrency>>& writer, [[maybe_unused]]G0Interpolation<device, Real>& g0)
    : QmciWalker(parameters, data_ref, rng, id),
      // QmciAutocorrelationData<QmciWalker>(parameters, id, g0),
      stamping_period_(parameters.stamping_period()),
      concurrency_id_(concurrency_id),
      thread_id_(id),
      writer_(writer),
      total_iterations_(parameters.get_dca_iterations()) {}

template <class QmciWalker, class DATA>
void StdThreadQmciWalker<QmciWalker, DATA>::initialize(int iteration) {
  QmciWalker::initialize(iteration);

  meas_id_ = 0;
  last_iteration_ = iteration == total_iterations_ - 1;
}

template <class QmciWalker, class DATA>
void StdThreadQmciWalker<QmciWalker, DATA>::doSweep() {
  QmciWalker::doSweep();
  // QmciAutocorrelationData<QmciWalker>::accumulateAutocorrelation(*this);

  if (QmciWalker::is_thermalized()) {
    // This must be before or the G_k_w and configuration meas_id will not match
    ++meas_id_;
    if (last_iteration_)
      logConfiguration();
  }
}

template <class QmciWalker, class DATA>
void StdThreadQmciWalker<QmciWalker, DATA>::logConfiguration() const {
  const bool print_to_log = writer_ && static_cast<bool>(*writer_);  // File exists and it is open. \todo possibly this should always be true
  if (print_to_log &&
      (writer_->isADIOS2() ||
       (writer_->isHDF5() && writer_->get_concurrency().id() == writer_->get_concurrency().first()))) {
    if (stamping_period_ && (meas_id_ % stamping_period_) == 0) {
      const std::string stamp_name = "r_" + std::to_string(concurrency_id_) + "_meas_" +
                                     std::to_string(meas_id_) + "_w_" + std::to_string(thread_id_);

      writer_->lock();
      writer_->open_group("STQW_Configurations");

      const auto& config = QmciWalker::get_configuration();
      config.write(*writer_, stamp_name);
      writer_->open_group(stamp_name);
      writer_->execute("log-weight", QmciWalker::get_MC_log_weight());
      writer_->close_group();
      writer_->close_group();
      writer_->unlock();
    }
  }
}

// No additional ops specialization for SS-HYB
template <linalg::DeviceType device, class Parameters, class DATA>
class StdThreadQmciWalker<cthyb::SsCtHybWalker<device, Parameters, DATA>, DATA>
    : public cthyb::SsCtHybWalker<device, Parameters, DATA>,
      public QmciAutocorrelationData<cthyb::SsCtHybWalker<device, Parameters, DATA>> {
  using QmciWalker = cthyb::SsCtHybWalker<device, Parameters, DATA>;
  using ThisType = StdThreadQmciWalker<QmciWalker, DATA>;
  using Concurrency = typename Parameters::concurrency_type;

  using Rng = typename Parameters::random_number_generator;
  using Real = typename QmciWalker::Scalar;

public:
  StdThreadQmciWalker(Parameters& parameters_ref, DATA& data_ref, Rng& rng, int concurrency_id,
                      int id, const std::shared_ptr<io::Writer<Concurrency>>& writer,
                      G0Interpolation<device, Real>& g0)
      : QmciWalker(parameters_ref, data_ref, rng, id),
        QmciAutocorrelationData<cthyb::SsCtHybWalker<device, Parameters, DATA>>(parameters_ref, id,
                                                                                g0),
        thread_id_(id),
        concurrency_id_(concurrency_id),
        writer_(writer) {}

  static void write([[maybe_unused]] io::Writer<Concurrency>& writer) {}

  int get_thread_id() const {
    return thread_id_;
  }

  std::size_t get_meas_id() const {
    return 0;
  }

private:
  int thread_id_;
  int concurrency_id_;
  std::shared_ptr<io::Writer<Concurrency>> writer_;
};

}  // namespace stdthreadqmci
}  // namespace solver
}  // namespace phys
}  // namespace dca


#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_WALKER_HPP
