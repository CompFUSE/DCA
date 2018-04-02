// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class provides the common interface between a walker on the CPU and one on the GPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_BASE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_BASE_HPP

#include <cassert>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <vector>

#include "dca/linalg/linalg.hpp"
#include "dca/linalg/make_constant_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/accumulator/ctint_accumulator_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/walker_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/domains/common_domains.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <dca::linalg::DeviceType device_t, class Parameters>
class CtintWalker;

template <dca::linalg::DeviceType device_t, class Parameters>
class CtintWalkerSubmatrix;

template <linalg::DeviceType device_t, class Parameters>
class CtintWalkerBase {
public:
  using parameters_type = Parameters;
  using Data = DcaData<Parameters>;
  using Rng = typename Parameters::random_number_generator;
  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename Parameters::concurrency_type;

protected:
  using Matrix = linalg::Matrix<double, device_t>;
  using MatrixView = linalg::MatrixView<double, device_t>;

protected:  // The class is not instantiable.
  CtintWalkerBase(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
                  const DMatrixBuilder<device_t>& builder_ref, int id = 0);

public:
  void markThermalized();

  bool is_thermalized() const {
    return thermalized_;
  }
  void doSweep();
  void do_step();

  int order() const {
    return configuration_.size();
  }
  double avgOrder() const {
    return total_steps_ ? double(order_sum_) / double(total_steps_) : order();
  }
  int sign() const {
    return sign_;
  }

  double acceptanceRatio() {
    return double(number_of_annihilations_ + number_of_creations_) / double(total_steps_);
  }

  void initialize() {}

  const auto& get_configuration() const {
    return configuration_;
  }

  void updateShell(int meas_id, int meas_to_do) const;

  void printSummary() const {}

protected:  // typedefs
  using TPosDmn = func::dmn_0<ctint::PositiveTimeDomain>;

protected:  // Auxiliary methods.
  virtual bool tryVertexInsert(bool forced = false) = 0;
  virtual bool tryVertexRemoval() = 0;

  template <linalg::DeviceType matrix_device>
  void applyInsertion(const MatrixPair<matrix_device>& S, const MatrixPair<matrix_device>& Q,
                      const MatrixPair<matrix_device>& R);
  double removalProbability();
  template <linalg::DeviceType matrix_device = linalg::CPU>
  void applyRemoval();

  void pushToEnd(const std::array<std::vector<ushort>, 2>& matrix_indices,
                 const std::pair<short, short>& vertex_indices);

  // For testing purposes
  void setMFromConfig();

  virtual void smallInverse(const MatrixView& in, MatrixView& out, int s) = 0;
  virtual void smallInverse(MatrixView& in_out, int s) = 0;
  virtual double separateIndexDeterminant(Matrix& m, const std::vector<ushort>& indices, int s) = 0;

protected:  // Members.
  Parameters& parameters_;
  Concurrency& concurrency_;

  const int thread_id_;

  Rng& rng_;
  SolverConfiguration<device_t> configuration_;

  MatrixPair<device_t> M_;
  MatrixPair<linalg::CPU> M_host_;

  const double beta_;
  const DMatrixBuilder<device_t>& d_builder_;

  const double total_interaction_;  // Space integrated interaction Hamiltonian.

  ulong order_sum_ = 0;
  ulong total_steps_ = 0;
  ulong total_sweeps_ = 0;
  ulong partial_order_sum_ = 0;
  ulong partial_num_steps_ = 0;
  ulong number_of_creations_ = 0;
  ulong number_of_annihilations_ = 0;
  int nb_steps_per_sweep_ = -1;

  bool thermalized_ = false;

  int sign_ = 1;

  // Store for testing purposes:
  std::array<double, 2> det_ratio_;

protected:
  // work sapces

  MatrixPair<device_t> M_Q_;
  MatrixPair<device_t> ws_dn_;
  std::array<std::vector<ushort>, 2> removal_matrix_indices_;
  std::pair<short, short> removal_candidates_;
  linalg::Vector<double, linalg::CPU> v_work_;
  linalg::Vector<int, linalg::CPU> ipiv_;
};

template <linalg::DeviceType device_t, class Parameters>
CtintWalkerBase<device_t, Parameters>::CtintWalkerBase(Parameters& parameters_ref, Rng& rng_ref,
                                                       const InteractionVertices& vertices,
                                                       const DMatrixBuilder<device_t>& builder_ref,
                                                       int id)
    : parameters_(parameters_ref),
      concurrency_(parameters_.get_concurrency()),

      thread_id_(id),

      rng_(rng_ref),

      configuration_(parameters_.get_beta(), Bdmn::dmn_size(), vertices,
                     parameters_.getDoubleUpdateProb()),

      beta_(parameters_.get_beta()),
      d_builder_(builder_ref),
      total_interaction_(vertices.integratedInteraction()),
      det_ratio_{1, 1} {}

template <linalg::DeviceType device_t, class Parameters>
void CtintWalkerBase<device_t, Parameters>::setMFromConfig() {
  // compute Mij = g0(t_i,t_j) - I* alpha(s_i)
  sign_ = 1;
  for (int s = 0; s < 2; ++s) {
    const auto& sector = configuration_.getSector(s);
    auto& M = M_host_[s];
    const int n = sector.size();
    M.resize(n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        M(i, j) = d_builder_.computeD(i, j, sector);
    linalg::matrixop::inverse(M);
    // Set the initial sign
    const double det = linalg::matrixop::determinant(M);
    if (det < 0)
      sign_ *= -1;
    if (device_t == linalg::GPU)
      M_[s].setAsync(M_host_[s], thread_id_, s);
  }

  if (device_t == linalg::CPU)
    for (int s = 0; s < 2; ++s)
      M_[s] = std::move(M_host_[s]);
  else
    for (int s = 0; s < 2; ++s)
      linalg::util::syncStream(thread_id_, s);
}

template <linalg::DeviceType device_t, class Parameters>
void CtintWalkerBase<device_t, Parameters>::doSweep() {
  int nb_of_steps;
  if (not thermalized_)
    nb_of_steps = avgOrder() + 1;
  else
    nb_of_steps = nb_steps_per_sweep_ * parameters_.get_sweeps_per_measurement();

  for (int i = 0; i < nb_of_steps; i++) {
    do_step();
    ++total_steps_;
    order_sum_ += order();
  }
  ++total_sweeps_;

  // Keep tha average after half tantalization for deciding the order.
  if ((not thermalized_) and (total_sweeps_ == parameters_.get_warm_up_sweeps() / 2)) {
    partial_order_sum_ = order_sum_;
    partial_num_steps_ = total_steps_;
  }
}

template <linalg::DeviceType device_t, class Parameters>
void CtintWalkerBase<device_t, Parameters>::markThermalized() {
  thermalized_ = true;
  nb_steps_per_sweep_ =
      dca::util::ceilDiv(order_sum_ - partial_order_sum_, total_steps_ - partial_num_steps_);

  order_sum_ = 0;
  total_steps_ = 0;
  total_sweeps_ = 0;
  number_of_creations_ = 0;
  number_of_annihilations_ = 0;
}

template <linalg::DeviceType device_t, class Parameters>
void CtintWalkerBase<device_t, Parameters>::do_step() {
  if (int(rng_() * 2)) {
    number_of_creations_ += tryVertexInsert();
  }
  else {
    if (configuration_.size())
      number_of_annihilations_ += tryVertexRemoval();
  }
}

template <linalg::DeviceType device_t, class Parameters>
void CtintWalkerBase<device_t, Parameters>::pushToEnd(
    const std::array<std::vector<ushort>, 2>& matrix_indices,
    const std::pair<short, short>& vertex_indices) {
  for (int s = 0; s < 2; ++s) {
    auto& M = M_[s];
    const auto indices = matrix_indices[s];
    ushort destination = M.nrCols() - 1;
    assert(M.nrCols() >= indices.size());
    // TODO check
    for (int idx = indices.size() - 1; idx >= 0; --idx) {
      const ushort source = indices[idx];
      linalg::matrixop::swapRowAndCol(M, source, destination, thread_id_, s);
      configuration_.swapSectorLabels(source, destination, s);
      --destination;
    }
  }

  if (vertex_indices.second == -1)
    configuration_.swapVertices(vertex_indices.first, configuration_.size() - 1);
  else {
    const ushort a = std::min(vertex_indices.first, vertex_indices.second);
    const ushort b = std::max(vertex_indices.first, vertex_indices.second);
    configuration_.swapVertices(b, configuration_.size() - 1);
    configuration_.swapVertices(a, configuration_.size() - 2);
  }
}

template <linalg::DeviceType device_t, class Parameters>
template <linalg::DeviceType matrix_device>
void CtintWalkerBase<device_t, Parameters>::applyInsertion(const MatrixPair<matrix_device>& Sp,
                                                           const MatrixPair<matrix_device>& Qp,
                                                           const MatrixPair<matrix_device>& Rp) {
  for (int s = 0; s < 2; ++s) {
    const int delta = Qp[s].nrCols();
    if (not delta)
      continue;
    // update M matrix.
    const auto& R = Rp[s];
    const auto S = linalg::makeConstantView(Sp[s]);
    auto& M = M_[s];
    const auto& M_Q = M_Q_[s];
    const int m_size = M.nrCols();

    if (not m_size) {
      M.resizeNoCopy(delta);
      auto M_view = MatrixView(M);
      smallInverse(*S, M_view, s);
      continue;
    }

    auto& R_M = ws_dn_[device_t == linalg::GPU ? s : 0];
    R_M.resizeNoCopy(R.size());
    linalg::matrixop::gemm(R, M, R_M, thread_id_, s);

    M.resize(m_size + delta);
    using MatrixView = linalg::MatrixView<double, device_t>;
    //  S_tilde = S^-1.
    MatrixView S_tilde(M, m_size, m_size, delta, delta);
    smallInverse(*S, S_tilde, s);

    // R_tilde = - S * R * M
    MatrixView R_tilde(M, m_size, 0, delta, m_size);
    linalg::matrixop::gemm(-1., S_tilde, R_M, double(0.), R_tilde, thread_id_, s);

    // Q_tilde = -M * Q * S
    MatrixView Q_tilde(M, 0, m_size, m_size, delta);
    linalg::matrixop::gemm(-1., M_Q, S_tilde, 0., Q_tilde, thread_id_, s);

    // update bulk: M += M*Q*S*R*M
    MatrixView M_bulk(M, 0, 0, m_size, m_size);
    linalg::matrixop::gemm(-1., Q_tilde, R_M, 1., M_bulk, thread_id_, s);
  }
}

template <linalg::DeviceType device_t, class Parameters>
double CtintWalkerBase<device_t, Parameters>::removalProbability() {
  removal_candidates_ = configuration_.randomRemovalCandidate(rng_);

  const int n = configuration_.size();
  const int n_removed = removal_candidates_.second == -1 ? 1 : 2;

  for (int s = 0; s < 2; ++s) {
    removal_matrix_indices_[s].clear();
    for (int i = 0; i < n_removed; ++i) {
      const int index = i ? removal_candidates_.second : removal_candidates_.first;
      const std::vector<ushort> vertex_indices = configuration_.findIndices(index, s);
      removal_matrix_indices_[s].insert(removal_matrix_indices_[s].end(), vertex_indices.begin(),
                                        vertex_indices.end());
    }
    if (removal_matrix_indices_[s].size() == 0) {
      det_ratio_[s] = 1.;
      continue;
    }
    std::sort(removal_matrix_indices_[s].begin(), removal_matrix_indices_[s].end());
    det_ratio_[s] = separateIndexDeterminant(M_[s], removal_matrix_indices_[s], s);
  }
  assert((removal_matrix_indices_[0].size() + removal_matrix_indices_[1].size()) % 2 == 0);

  const int combinatorial_factor =
      n_removed == 1 ? n : n * configuration_.occupationNumber(removal_candidates_.second);

  return details::computeAcceptanceProbability(n_removed, det_ratio_[0] * det_ratio_[1],
                                               total_interaction_, beta_,
                                               configuration_.getStrength(removal_candidates_.first),
                                               combinatorial_factor, details::VERTEX_REMOVAL);
}

template <linalg::DeviceType device_t, class Parameters>
template <linalg::DeviceType matrix_device>
void CtintWalkerBase<device_t, Parameters>::applyRemoval() {
  const int n_removed = removal_candidates_.second == -1 ? 1 : 2;
  // TODO maybe: don't copy values to be popped.
  pushToEnd(removal_matrix_indices_, removal_candidates_);
  configuration_.pop(n_removed);

  for (int s = 0; s < 2; ++s) {
    auto& M = M_[s];
    const int delta = removal_matrix_indices_[s].size();
    if (delta == 0)  // Nothing to update.
      continue;
    const int m_size = M.nrCols() - delta;
    if (m_size == 0) {  // removing last vertex
      M.resize(0);
      continue;
    }

    MatrixView Q(M, 0, m_size, m_size, delta);
    MatrixView R(M, m_size, 0, delta, m_size);

    MatrixView S(M, m_size, m_size, delta, delta);
    smallInverse(S, s);

    auto& Q_S = M_Q_[s];
    Q_S.resizeNoCopy(Q.size());
    linalg::matrixop::gemm(Q, S, Q_S, thread_id_, s);

    // M -= Q*S^-1*R
    MatrixView M_bulk(M, 0, 0, m_size, m_size);
    linalg::matrixop::gemm(-1., Q_S, R, 1., M_bulk, thread_id_, s);
    M.resize(m_size);
  }
}

template <linalg::DeviceType device_t, class Parameters>
void CtintWalkerBase<device_t, Parameters>::updateShell(int meas_id, int meas_to_do) const {
  if (concurrency_.id() == concurrency_.first() && meas_id > 1 && (meas_id % (meas_to_do / 10)) == 0) {
    std::cout << "\t\t\t" << int(double(meas_id) / double(meas_to_do) * 100) << " % completed \t ";
    std::cout << "\t k :" << order();
    const double avg_order = avgOrder();
    if (avg_order != -1)
      std::cout << "\t <k> :" << std::setprecision(1) << std::fixed << avg_order;
    std::cout << "\t\t" << dca::util::print_time() << "\n";
  }
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_BASE_HPP
