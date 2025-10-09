// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This class provides the common interface between a walker on the CPU and one on the GPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_BASE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_BASE_HPP

#include <cassert>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <vector>

#include "dca/distribution/dist_types.hpp"
#include "dca/io/buffer.hpp"
#include "dca/math/util/phase.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/linalg/util/gpu_stream.hpp"
#include "dca/linalg/util/stream_container.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/walker_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/domains/common_domains.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/util/print_time.hpp"

#ifdef DCA_HAVE_GPU
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"
#endif

//#define DEBUG_SUBMATRIX

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <linalg::DeviceType device_type, class Parameters, DistType DIST>
class CtintWalker;

template <class Parameters, DistType DIST = DistType::NONE>
class CtintWalkerBase {
public:
  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;
  using parameters_type = Parameters;
  using Data = DcaData<Parameters, DIST>;
  using Rng = typename Parameters::random_number_generator;
  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename Parameters::concurrency_type;

  using Matrix = linalg::Matrix<Scalar, linalg::CPU>;
  using MatrixPair = std::array<Matrix, 2>;
  using GpuStream = linalg::util::GpuStream;

  constexpr static linalg::DeviceType device = linalg::CPU;

protected:  // The class is not instantiable.
  CtintWalkerBase(const Parameters& pars_ref, Rng& rng_ref, int id = 0);

  virtual ~CtintWalkerBase() = default;

public:
  const auto& getConfiguration() const {
    return configuration_;
  }

  void computeM(MatrixPair& m_accum) const;

  // Reset the counters and recompute the configuration sign and weight.
  virtual void markThermalized() = 0;

  /** Recompute the matrix M from the configuration in O(expansion_order^3) time.
   *  Postcondition: sign_ and mc_log_weight_ are recomputed.
   *  mc_log_weight is the negative sum of the log det of both sectors of M
   *  + log of the abs of each vertices interaction strength.
   */
  virtual void setMFromConfig() = 0;

  virtual void doStep() = 0;

  virtual void doSweep() = 0;

  template <linalg::DeviceType DEVICE>
  void setMFromConfigImpl(DMatrixBuilder<DEVICE, Scalar>& d_matrix_builder);

  bool is_thermalized() const {
    return thermalized_;
  }

  unsigned long get_steps() const {
    return n_steps_;
  }

  int order() const {
    return configuration_.size();
  }
  double avgOrder() const {
    return order_avg_.count() ? order_avg_.mean() : order();
  }
  auto get_sign() const {
    return phase_;
  }

  Real get_MC_log_weight() const {
    return mc_log_weight_;
  }

  double acceptanceRatio() const {
    return Real(n_accepted_) / Real(n_steps_ - thermalization_steps_);
  }

  void initialize(int iter);

  const auto& get_configuration() const {
    return configuration_;
  }

  const auto& get_matrix_configuration() const {
    return configuration_.get_sectors();
  }

  void updateShell(int meas_id, int meas_to_do) const;

  void printSummary() const;

  virtual void synchronize() const {}

  // For testing purposes.
  // Fixes the numbers of proposed steps per sweep.
  void fixStepsPerSweep(const int nb_steps_per_sweep) {
    assert(nb_steps_per_sweep > 0);
    nb_steps_per_sweep_ = nb_steps_per_sweep;
  }

  virtual std::size_t deviceFingerprint() const {
    return 0;
  };

  io::Buffer dumpConfig() const {
    io::Buffer buff;
    buff << configuration_;
    return buff;
  }

  void readConfig(io::Buffer& buff) {
    buff >> configuration_;
  }

  double stealFLOPs() {
    auto flop = flop_;
    flop_ = 0.;
    return flop;
  }

  const auto& get_stream() const {
    return stream_;
  }

  static void sumConcurrency(const Concurrency&) {}

  void writeAlphas() const;

  static void setInteractionVertices(const Data& data, const Parameters& parameters);

protected:
  // typedefs
  using RDmn = typename Parameters::RClusterDmn;

  // Auxiliary methods.
  void updateSweepAverages();

protected:  // Members.
  static inline InteractionVertices vertices_;

  const Parameters& parameters_;
  const Concurrency& concurrency_;

  const int thread_id_;
  linalg::util::GpuStream* stream_;

  Rng& rng_;
  SolverConfiguration configuration_;

  MatrixPair M_;

  const Real beta_;
  static inline constexpr int n_bands_ = Parameters::bands;

  const Real total_interaction_;  // Space integrated interaction Hamiltonian.

  util::Accumulator<uint> partial_order_avg_;
  util::Accumulator<uint> order_avg_;
  util::Accumulator<math::Phase<Scalar>> sign_avg_;
  unsigned long n_steps_ = 0;
  unsigned long thermalization_steps_ = 0;
  unsigned long n_accepted_ = 0;
  int nb_steps_per_sweep_ = -1;

  bool thermalized_ = false;

  math::Phase<Scalar> phase_;

  // Store for testing purposes:
  Scalar acceptance_prob_;

  double flop_ = 0.;

  double sweeps_per_meas_ = 1.;

  double mc_log_weight_ = 0;

private:
  linalg::Vector<int, linalg::CPU> ipiv_;
  linalg::Vector<Scalar, linalg::CPU> work_;
};

template <class Parameters, DistType DIST>
CtintWalkerBase<Parameters, DIST>::CtintWalkerBase(const Parameters& parameters_ref, Rng& rng_ref,
                                                   int id)
    : parameters_(parameters_ref),
      concurrency_(parameters_.get_concurrency()),
      thread_id_(id),
      stream_{&linalg::util::getStreamContainer()(thread_id_, 0)},
      rng_(rng_ref),

      configuration_(parameters_.get_beta(), Bdmn::dmn_size(), vertices_,
                     parameters_.getDoubleUpdateProbability()),

      beta_(parameters_.get_beta()),
      total_interaction_(vertices_.integratedInteraction()) {}

template <class Parameters, DistType DIST>
void CtintWalkerBase<Parameters, DIST>::initialize(int iteration) {
  assert(total_interaction_);
  phase_.reset();

  mc_log_weight_ = 0.;

  sweeps_per_meas_ = parameters_.get_sweeps_per_measurement().at(iteration);

  if (!configuration_.size()) {  // Do not initialize config if it was read.
    while (parameters_.getInitialConfigurationSize() > configuration_.size()) {
      configuration_.insertRandom(rng_);
      for (int i = configuration_.lastInsertionSize(); i > 0; --i)
        configuration_.commitInsertion(configuration_.size() - i);
    }
  }

  setMFromConfig();
}

template <class Parameters, DistType DIST>
void CtintWalkerBase<Parameters, DIST>::updateSweepAverages() {
  order_avg_.addSample(order());
  sign_avg_.addSample(phase_.getSign());
  // Track avg order for the final number of steps / sweep.
  if (!thermalized_ && order_avg_.count() >= parameters_.get_warm_up_sweeps() / 2)
    partial_order_avg_.addSample(order());
}

// template <class Parameters, DistType DIST>
// void CtintWalkerBase<Parameters,DIST>::writeAlphas() const {
//   std::cout << "For initial configuration integration:\n";
//   for (int isec = 0; isec < 2; ++isec) {
//     std::cout << "Sector: " << isec << '\n';
//     for (int ic = 0; ic < order(); ++ic) {
//       auto aux_spin = configuration_.getSector(isec).getAuxFieldType(ic);
//       auto left_b = configuration_.getSector(isec).getLeftB(ic);
//       auto alpha_left = get_dmatrix_builder().computeAlpha(aux_spin, left_b);
//       std::cout << "vertex: " << std::setw(6) << ic;
//       std::cout << " | aux spin: " << aux_spin << " | left B: " << left_b
//                 << " | alpha left = " << alpha_left << '\n';
//     }
//   }
// }

template <class Parameters, DistType DIST>
void CtintWalkerBase<Parameters, DIST>::updateShell(int meas_id, int meas_to_do) const {
  if (concurrency_.id() == concurrency_.first() && meas_id > 1 &&
      (meas_id % dca::util::ceilDiv(meas_to_do, 20)) == 0) {
    std::cout << "\t\t\t" << int(double(meas_id) / double(meas_to_do) * 100) << " % completed \t ";
    std::cout << "\t k :" << order();
    const double avg_order = avgOrder();
    if (avg_order != -1) {
      auto precision = std::cout.precision();
      std::cout.precision(1);
      std::cout << "\t <k> :" << std::fixed << avg_order;
      std::cout.precision(precision);
    }
    std::cout << "\t\t" << dca::util::print_time() << "\n";
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerBase<Parameters, DIST>::printSummary() const {
  std::cout << "\n"
            << "Walker: process ID = " << concurrency_.id() << ", thread ID = " << thread_id_ << "\n"
            << "-------------------------------------------\n";

  if (partial_order_avg_.count())
    std::cout << "Estimate for sweep size: " << partial_order_avg_.mean() << "\n";
  if (order_avg_.count())
    std::cout << "Average expansion order: " << order_avg_.mean() << "\n";
  if (sign_avg_.count())
    std::cout << "Average sign: " << sign_avg_.mean() << "\n";

  std::cout << "Acceptance ratio: " << acceptanceRatio() << "\n";
  std::cout << std::endl;
}

template <class Parameters, DistType DIST>
void CtintWalkerBase<Parameters, DIST>::setInteractionVertices(const Data& data,
                                                               const Parameters& parameters) {
  vertices_.reset();
  vertices_.initialize(parameters.getDoubleUpdateProbability(), parameters.getAllSitesPartnership());
  vertices_.initializeFromHamiltonian(data.H_interactions);
  if (data.has_non_density_interactions()) {
    vertices_.checkForInterbandPropagators(data.G0_r_t_cluster_excluded);
    vertices_.initializeFromNonDensityHamiltonian(data.get_non_density_interactions());
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerBase<Parameters, DIST>::computeM(MatrixPair& m_accum) const {
  m_accum = M_;
}

// template<class WALKER, linalg::DeviceType DEVICE>
// void setMFromConfigHelper(WALKER& walker, DMatrixBuilder<DEVICE, Scalar>& d_matrix_builder) {
//   walker.mc_log_weight_ = 0.;
//   walker.phase_.reset();

//   for (int s = 0; s < 2; ++s) {
//     // compute Mij = g0(t_i,t_j) - I* alpha(s_i)

//     const auto& sector = walker.configuration_.getSector(s);
//     auto& M = walker.M_[s];
//     const int n = sector.size();
//     M.resize(n);
//     if (!n)
//       continue;
//     for (int j = 0; j < n; ++j)
//       for (int i = 0; i < n; ++i)
//         M(i, j) = d_matrix_builder.computeD(i, j, sector);

//     if (M.nrRows()) {
//       const auto [log_det, phase] = linalg::matrixop::inverseAndLogDeterminant(M);

//       walker.mc_log_weight_ -= log_det;  // Weight proportional to det(M^{-1})
//       walker.phase_.divide(phase);
//     }
//   }

//   // So what is going on here.
//   for (int i = 0; i < walker.configuration_.size(); ++i) {
//     // This is actual interaction strength of the vertex i.e H_int(nu1, nu2, delta_r)
//     const typename decltype(walker)::Real term = -walker.configuration_.getStrength(i);
//     walker.mc_log_weight_ += std::log(std::abs(term));
//     walker.phase_.multiply(term);
//   }
// }

template <class Parameters, DistType DIST>
template <linalg::DeviceType DEVICE>
void CtintWalkerBase<Parameters, DIST>::setMFromConfigImpl(
    DMatrixBuilder<DEVICE, Scalar>& d_matrix_builder) {
  mc_log_weight_ = 0.;
  phase_.reset();

  for (int s = 0; s < 2; ++s) {
    // compute Mij = g0(t_i,t_j) - I* alpha(s_i)

    const auto& sector = configuration_.getSector(s);
    auto& M = M_[s];
    const int n = sector.size();
    M.resize(n);
    if (!n)
      continue;
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        M(i, j) = d_matrix_builder.computeD(i, j, sector);

    if (M.nrRows()) {
      const auto [log_det, phase] = linalg::matrixop::inverseAndLogDeterminant(M);

      mc_log_weight_ -= log_det;  // Weight proportional to det(M^{-1})
      phase_.divide(phase);
    }
  }

  // So what is going on here.
  for (int i = 0; i < configuration_.size(); ++i) {
    // This is actual interaction strength of the vertex i.e H_int(nu1, nu2, delta_r)
    const Real term = -configuration_.getStrength(i);
    mc_log_weight_ += std::log(std::abs(term));
    phase_.multiply(term);
  }
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_BASE_HPP
