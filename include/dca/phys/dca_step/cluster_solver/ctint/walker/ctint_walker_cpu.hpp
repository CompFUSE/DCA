// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class organizes the MC walker in the CT-INT QMC purely on the CPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_HPP

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_base.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/walker_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <class Parameters>
class CtintWalker<linalg::CPU, Parameters> : public CtintWalkerBase<linalg::CPU, Parameters> {
public:
  using this_type = CtintWalker<linalg::CPU, Parameters>;
  using BaseClass = CtintWalkerBase<linalg::CPU, Parameters>;
  using Rng = typename BaseClass::Rng;

public:
  CtintWalker(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
              const DMatrixBuilder<linalg::CPU>& builder_ref, int id = 0);

  AccumulatorConfiguration getConfiguration() const;
  AccumulatorConfiguration moveConfiguration();
  void setConfiguration(AccumulatorConfiguration&& config);

protected:
  bool tryVertexInsert(bool forced = false);
  bool tryVertexRemoval();

  // Test handle.
  const auto& getM() const {
    return M_;
  }

private:
  double insertionProbability(int delta_vertices);

  using MatrixView = typename BaseClass::MatrixView;
  using Matrix = typename BaseClass::Matrix;
  virtual void smallInverse(const MatrixView& in, MatrixView& out, int s);
  virtual void smallInverse(MatrixView& in_out, int s);
  virtual double separateIndexDeterminant(Matrix& m, const std::vector<ushort>& indices, int s);

protected:
  using BaseClass::parameters_;
  using BaseClass::configuration_;
  using BaseClass::rng_;
  using BaseClass::d_builder_;
  using BaseClass::total_interaction_;
  using BaseClass::beta_;
  using BaseClass::sign_;

protected:
  using BaseClass::M_;
  using BaseClass::det_ratio_;

private:
  std::array<linalg::Matrix<double, linalg::CPU>, 2> S_, Q_, R_;
  // work spaces
  using BaseClass::M_Q_;
};

template <class Parameters>
CtintWalker<linalg::CPU, Parameters>::CtintWalker(Parameters& parameters_ref, Rng& rng_ref,
                                                  const InteractionVertices& vertices,
                                                  const DMatrixBuilder<linalg::CPU>& builder_ref,
                                                  int id)
    : BaseClass(parameters_ref, rng_ref, vertices, builder_ref, id) {
  while (parameters_.getInitialConfigurationSize() > configuration_.size())
    tryVertexInsert(true);
}

template <class Parameters>
AccumulatorConfiguration CtintWalker<linalg::CPU, Parameters>::getConfiguration() const {
  return AccumulatorConfiguration{sign_, M_, configuration_};
}
template <class Parameters>
AccumulatorConfiguration CtintWalker<linalg::CPU, Parameters>::moveConfiguration() {
  return AccumulatorConfiguration{sign_, std::move(M_), std::move(configuration_)};
}

template <class Parameters>
void CtintWalker<linalg::CPU, Parameters>::setConfiguration(AccumulatorConfiguration&& config) {
  sign_ = config.sign;
  M_ = std::move(config.M);
  static_cast<MatrixConfiguration&>(configuration_) = std::move(config.matrix_configuration);
}

template <class Parameters>
bool CtintWalker<linalg::CPU, Parameters>::tryVertexInsert(bool forced) {
  configuration_.insertRandom(rng_);
  const int delta_vertices = configuration_.lastInsertionSize();

  // Compute the new pieces of the D(= M^-1) matrix.
  d_builder_.buildSQR(S_, Q_, R_, configuration_);

  const double accept_prob = insertionProbability(delta_vertices);

  const bool accept = std::abs(accept_prob) > rng_() or (forced and accept_prob != 0);

  if (not accept)
    configuration_.pop(delta_vertices);
  else {
    if (accept_prob < 0)
      sign_ *= -1;
    BaseClass::applyInsertion(S_, Q_, R_);
  }
  return accept;
}

template <class Parameters>
bool CtintWalker<linalg::CPU, Parameters>::tryVertexRemoval() {
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
double CtintWalker<linalg::CPU, Parameters>::insertionProbability(const int delta_vertices) {
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
      linalg::matrixop::gemm(M, Q, M_Q);
      // S <- S_tilde^(-1) = S - R*M*Q
      linalg::matrixop::gemm(-1., R, M_Q, 1., S);
    }

    det_ratio_[s] = details::smallDeterminant(S);
  }

  const int combinatorial_factor =
      delta_vertices == 1 ? old_size + 1
                          : (old_size + 2) * configuration_.occupationNumber(old_size + 1);

  return details::computeAcceptanceProbability(
      delta_vertices, det_ratio_[0] * det_ratio_[1], total_interaction_, beta_,
      configuration_.getStrength(old_size), combinatorial_factor, details::VERTEX_INSERTION);
}

template <class Parameters>
void CtintWalker<linalg::CPU, Parameters>::smallInverse(const MatrixView& in, MatrixView& out,
                                                        const int s) {
  details::smallInverse(in, out, det_ratio_[s], BaseClass::ipiv_, BaseClass::v_work_);
}
template <class Parameters>
void CtintWalker<linalg::CPU, Parameters>::smallInverse(MatrixView& in_out, const int s) {
  details::smallInverse(in_out, det_ratio_[s], BaseClass::ipiv_, BaseClass::v_work_);
}

template <class Parameters>
double CtintWalker<linalg::CPU, Parameters>::separateIndexDeterminant(
    Matrix& m, const std::vector<ushort>& indices, int /*s*/) {
  return details::separateIndexDeterminant(m, indices);
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_HPP
