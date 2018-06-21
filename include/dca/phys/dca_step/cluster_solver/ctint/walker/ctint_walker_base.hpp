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
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/util/print_time.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"
#endif

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <linalg::DeviceType device_type, class Parameters>
class CtintWalker;

template <linalg::DeviceType device_type, class Parameters>
class CtintWalkerSubmatrix;

template <class Parameters>
class CtintWalkerBase {
public:
  using parameters_type = Parameters;
  using Data = DcaData<Parameters>;
  using Rng = typename Parameters::random_number_generator;
  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename Parameters::concurrency_type;

protected:  // The class is not instantiable.
  CtintWalkerBase(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
                  const DMatrixBuilder<linalg::CPU>& builder_ref, int id = 0);

public:
  AccumulatorConfiguration getConfiguration() const;
  AccumulatorConfiguration moveConfiguration();
  void setConfiguration(AccumulatorConfiguration&& config);

  void markThermalized();

  bool is_thermalized() const {
    return thermalized_;
  }

  int order() const {
    return configuration_.size();
  }
  double avgOrder() const {
    return order_avg_.count() ? order_avg_.mean() : order();
  }
  int sign() const {
    return sign_;
  }

  double acceptanceRatio() {
    return double(n_accepted_) / double(n_steps_);
  }

  void initialize() {}

  const auto& get_configuration() const {
    return configuration_;
  }

  void updateShell(int meas_id, int meas_to_do) const;

  void printSummary() const {}

  virtual void synchronize() const {}

protected:  // typedefs
  using Matrix = linalg::Matrix<double, linalg::CPU>;
  using MatrixPair = std::array<linalg::Matrix<double, linalg::CPU>, 2>;

  using TPosDmn = func::dmn_0<ctint::PositiveTimeDomain>;

protected:  // Auxiliary methods.
  void updateSweepAverages();

  void pushToEnd(const std::array<std::vector<ushort>, 2>& matrix_indices,
                 const std::pair<short, short>& vertex_indices);

  void setMFromConfig();

protected:  // Members.
  Parameters& parameters_;
  Concurrency& concurrency_;

  const int thread_id_;

  Rng& rng_;
  SolverConfiguration configuration_;

  MatrixPair M_;

  const double beta_;
  static constexpr int n_bands_ = Parameters::bands;

  const DMatrixBuilder<linalg::CPU>& d_builder_;

  const double total_interaction_;  // Space integrated interaction Hamiltonian.

  util::Accumulator<uint> partial_order_avg_;
  util::Accumulator<uint> order_avg_;
  ulong n_steps_ = 0;
  ulong n_accepted_ = 0;
  int nb_steps_per_sweep_ = -1;

  bool thermalized_ = false;

  int sign_ = 1;

  // Store for testing purposes:
  double acceptance_prob_;

  std::array<std::vector<ushort>, 2> removal_matrix_indices_;
  std::pair<short, short> removal_candidates_;

private:
  linalg::Vector<int, linalg::CPU> ipiv_;
  linalg::Vector<double, linalg::CPU> work_;
};

template <class Parameters>
CtintWalkerBase<Parameters>::CtintWalkerBase(Parameters& parameters_ref, Rng& rng_ref,
                                                       const InteractionVertices& vertices,
                                                       const DMatrixBuilder<linalg::CPU>& builder_ref,
                                                       int id)
    : parameters_(parameters_ref),
      concurrency_(parameters_.get_concurrency()),

      thread_id_(id),

      rng_(rng_ref),

      configuration_(parameters_.get_beta(), Bdmn::dmn_size(), vertices,
                     parameters_.getDoubleUpdateProb()),

      beta_(parameters_.get_beta()),
      d_builder_(builder_ref),
      total_interaction_(vertices.integratedInteraction())
      {
  while (parameters_.getInitialConfigurationSize() > configuration_.size())
    configuration_.insertRandom(rng_);

  setMFromConfig();
}

template <class Parameters>
void CtintWalkerBase<Parameters>::setMFromConfig() {
  // compute Mij = g0(t_i,t_j) - I* alpha(s_i)
  sign_ = 1;
  for (int s = 0; s < 2; ++s) {
    const auto& sector = configuration_.getSector(s);
    auto& M = M_[s];
    const int n = sector.size();
    M.resize(n);
    if (!n)
      continue;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        M(i, j) = d_builder_.computeD(i, j, sector);

    const double det = linalg::matrixop::inverseAndDeterminant(M);

    // Set the initial sign
    if (det < 0)
      sign_ *= -1;
  }
}

template <class Parameters>
AccumulatorConfiguration CtintWalkerBase<Parameters>::getConfiguration() const {
  synchronize();
  return AccumulatorConfiguration{sign_, M_, configuration_};
}

template <class Parameters>
AccumulatorConfiguration CtintWalkerBase<Parameters>::moveConfiguration() {
  return AccumulatorConfiguration{sign_, std::move(M_), std::move(configuration_)};
}

template <class Parameters>
void CtintWalkerBase<Parameters>::setConfiguration(AccumulatorConfiguration&& config) {
  sign_ = config.sign;
  M_ = std::move(config.M);
  static_cast<MatrixConfiguration&>(configuration_) = std::move(config.matrix_configuration);
}

template <class Parameters>
void CtintWalkerBase<Parameters>::updateSweepAverages() {
    order_avg_.addSample(order());
  // Track avg order for the final number of steps / sweep.
  if (order_avg_.count() >= parameters_.get_warm_up_sweeps() / 2)
    partial_order_avg_.addSample(order());
}

template <class Parameters>
void CtintWalkerBase<Parameters>::markThermalized() {
  thermalized_ = true;
  nb_steps_per_sweep_ = std::ceil(partial_order_avg_.mean());

  order_avg_.reset();
  n_accepted_ = 0;
  n_steps_ = 0;
}

template <class Parameters>
void CtintWalkerBase<Parameters>::pushToEnd(
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

template <class Parameters>
void CtintWalkerBase<Parameters>::updateShell(int meas_id, int meas_to_do) const {
  if (concurrency_.id() == concurrency_.first() && meas_id > 1 && (meas_id % dca::util::ceilDiv(meas_to_do, 10)) == 0) {
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
