// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the MC walker in the CT-INT QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_base.hpp"

#include "dca/linalg/vector.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/kernels_interface.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <class Parameters>
class CtintWalker<linalg::GPU, Parameters> : public CtintWalkerBase<linalg::GPU, Parameters> {
public:
  using this_type = CtintWalker<linalg::GPU, Parameters>;
  using BaseClass = CtintWalkerBase<linalg::GPU, Parameters>;
  using Rng = typename BaseClass::Rng;

public:
  CtintWalker(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
              const DMatrixBuilder<linalg::GPU>& builder_ref, int id = 0);

  AccumulatorConfiguration getConfiguration();
  AccumulatorConfiguration moveConfiguration();
  void setConfiguration(AccumulatorConfiguration&& config);

protected:
  bool tryVertexInsert(bool forced = false);
  bool tryVertexRemoval();

  double insertionProbability(int delta_vertices);

  using MatrixView = typename BaseClass::MatrixView;
  using Matrix = typename BaseClass::Matrix;
  virtual void smallInverse(const MatrixView& in, MatrixView& out, int s);
  virtual void smallInverse(MatrixView& in_out, int s);
  virtual double separateIndexDeterminant(Matrix& m, const std::vector<ushort>& indices, int s);

  // Test handle.
  const auto& getM() {
    for (int s = 0; s < 2; ++s)
      M_host_[s].setAsync(M_[s], streams_[s]);
    for (int s = 0; s < 2; ++s)
      cudaStreamSynchronize(streams_[s]);
    return M_host_;
  }

private:  // Members.
  using BaseClass::parameters_;
  using BaseClass::thread_id_;
  using BaseClass::configuration_;
  using BaseClass::rng_;
  using BaseClass::d_builder_;
  using BaseClass::total_interaction_;
  using BaseClass::sign_;
  using BaseClass::beta_;
  using BaseClass::M_;
  using BaseClass::M_Q_;

protected:
  using BaseClass::det_ratio_;

private:  // Work spaces specific to the GPU implementation.
  std::array<linalg::Matrix<double, linalg::GPU>, 2> S_, Q_, R_;
  std::array<linalg::Matrix<double, linalg::CPU>, 2> S_host_, M_host_;
  std::array<cudaStream_t, 2> streams_;
  std::array<linalg::Vector<ushort, linalg::GPU>, 2> indices_;
  std::array<linalg::Vector<double, linalg::GPU>, 2> dev_det_;
};

template <class Parameters>
CtintWalker<linalg::GPU, Parameters>::CtintWalker(Parameters& parameters_ref, Rng& rng_ref,
                                                  const InteractionVertices& vertices,
                                                  const DMatrixBuilder<linalg::GPU>& builder_ref,
                                                  int id)
    : BaseClass(parameters_ref, rng_ref, vertices, builder_ref, id),
      streams_{linalg::util::getStream(thread_id_, 0), linalg::util::getStream(thread_id_, 1)} {
  for (auto& det : dev_det_)
    det.resize(1);
  while (parameters_.getInitialConfigurationSize() > configuration_.size())
    tryVertexInsert(true);
}

template <class Parameters>
AccumulatorConfiguration CtintWalker<linalg::GPU, Parameters>::getConfiguration() {
  for (int s = 0; s < 2; ++s)
    M_host_[s].setAsync(M_[s], streams_[s]);
  for (int s = 0; s < 2; ++s)
    cudaStreamSynchronize(streams_[s]);
  return AccumulatorConfiguration{sign_, M_host_, configuration_};
}
template <class Parameters>
AccumulatorConfiguration CtintWalker<linalg::GPU, Parameters>::moveConfiguration() {
  throw(std::logic_error("Not implemented."));
}
template <class Parameters>
void CtintWalker<linalg::GPU, Parameters>::setConfiguration(AccumulatorConfiguration&& /*config*/) {
  throw(std::logic_error("Not implemented."));
}

template <class Parameters>
bool CtintWalker<linalg::GPU, Parameters>::tryVertexInsert(bool forced) {
  configuration_.insertRandom(rng_);
  const int delta_vertices = configuration_.lastInsertionSize();

  // Compute the new pieces of the D(= M^-1) matrix.
  d_builder_.buildSQR(S_, Q_, R_, configuration_, thread_id_);

  const double accept_prob = insertionProbability(delta_vertices);
  const bool accept = std::abs(accept_prob) > rng_() or (forced and accept_prob != 0);

  if (accept) {
    if (accept_prob < 0)
      sign_ *= -1;
    BaseClass::applyInsertion(S_, Q_, R_);
  }
  else {
    configuration_.pop(delta_vertices);
  }

  return accept;
}

template <class Parameters>
bool CtintWalker<linalg::GPU, Parameters>::tryVertexRemoval() {
  const double accept_prob = BaseClass::removalProbability();
  const bool accept = std::abs(accept_prob) > rng_();

  if (accept) {
    if (accept_prob < 0)
      sign_ *= -1;

    BaseClass::applyRemoval();
  }
  return accept;
}

template <class Parameters>
double CtintWalker<linalg::GPU, Parameters>::insertionProbability(const int delta_vertices) {
  const int old_size = configuration_.size() - delta_vertices;

  for (int s = 0; s < 2; ++s) {
    const auto& Q = Q_[s];
    if (Q.nrCols() == 0) {
      det_ratio_[s] = 1.;
      continue;
    }
    const auto& R = R_[s];
    auto& S = S_[s];
    auto& M = M_[s];

    if (M.nrRows()) {
      auto& M_Q = M_Q_[s];
      M_Q.resizeNoCopy(Q.size());
      linalg::matrixop::gemm(M_[s], Q, M_Q, thread_id_, s);
      // S <- S_tilde^(-1) = S - R*M*Q
      linalg::matrixop::gemm(-1., R, M_Q, 1., S, thread_id_, s);
    }

    S_host_[s].setAsync(S, streams_[s]);
  }
  for (int s = 0; s < 2; ++s) {
    cudaStreamSynchronize(streams_[s]);
    det_ratio_[s] = details::smallDeterminant(S_host_[s]);
  }

  const int combinatorial_factor =
      delta_vertices == 1 ? old_size + 1
                          : (old_size + 2) * configuration_.occupationNumber(old_size + 1);

  return details::computeAcceptanceProbability(
      delta_vertices, det_ratio_[0] * det_ratio_[1], total_interaction_, beta_,
      configuration_.getStrength(old_size), combinatorial_factor, details::VERTEX_INSERTION);
}

template <class Parameters>
void CtintWalker<linalg::GPU, Parameters>::smallInverse(const MatrixView& in, MatrixView& out,
                                                        const int s) {
  details::smallInverse(in, out, streams_[s]);
}

template <class Parameters>
void CtintWalker<linalg::GPU, Parameters>::smallInverse(MatrixView& in_out, const int s) {
  details::smallInverse(in_out, streams_[s]);
}

template <class Parameters>
double CtintWalker<linalg::GPU, Parameters>::separateIndexDeterminant(
    Matrix& m, const std::vector<ushort>& indices, int s) {
  indices_[s].setAsync(indices, streams_[s]);
  return details::separateIndexDeterminant(MatrixView(m), indices_[s].ptr(), indices_[s].size(),
                                           dev_det_[s].ptr(), streams_[s]);
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_HPP
