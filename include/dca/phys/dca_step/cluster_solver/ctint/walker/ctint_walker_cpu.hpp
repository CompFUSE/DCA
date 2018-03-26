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

  CtintWalker(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
              const DMatrixBuilder<linalg::CPU>& builder_ref, int id = 0);

protected:
  bool tryVertexInsert(bool forced = false);
  bool tryVertexRemoval();

  using BaseClass::det_ratio_;

private:
  using Matrix = typename BaseClass::Matrix;
  using MatrixView = linalg::MatrixView<double, linalg::CPU>;

  using BaseClass::parameters_;
  using BaseClass::configuration_;
  using BaseClass::rng_;
  using BaseClass::d_builder_;
  using BaseClass::total_interaction_;
  using BaseClass::beta_;
  using BaseClass::sign_;
  using BaseClass::M_;

private:
  MatrixPair<linalg::CPU> S_, Q_, R_;
  // work spaces
  using BaseClass::ipiv_;
  using BaseClass::v_work_;
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
bool CtintWalker<linalg::CPU, Parameters>::tryVertexInsert(bool forced) {
  configuration_.insertRandom(rng_);
  const int delta_vertices = configuration_.lastInsertionSize();

  // Compute the new pieces of the D(= M^-1) matrix.
  d_builder_.buildSQR(S_, Q_, R_, configuration_);

  const double accept_prob = BaseClass::insertionProbability(S_, Q_, R_, delta_vertices);

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

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_HPP
