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

#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration_gpu.hpp"

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

  CtintWalker(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
              const DMatrixBuilder<linalg::GPU>& builder_ref, int id = 0);

private:
  bool tryVertexInsert(bool forced = false);
  bool tryVertexRemoval();

  void pruneInsertionMatrices(int stram_idx, int spin);
  void prepareInsertionMatrices(const int stream);

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
  using BaseClass::det_ratio_;

private:  // Work spaces specific to the GPU implementation.
  using Matrix = dca::linalg::Matrix<double, linalg::CPU>;
  using MatrixView = dca::linalg::MatrixView<double, linalg::CPU>;

  std::array<MatrixPair<linalg::CPU>, 2> Sa_, Qa_, Ra_;
  DeviceWorkspace devspace;
  std::array<std::vector<ushort>, 2> removed_indices_;
  ushort stream_index_ = 0;
  std::array<SolverConfiguration<linalg::GPU>, 2> configuration_copy_;
};

template <class Parameters>
CtintWalker<linalg::GPU, Parameters>::CtintWalker(Parameters& parameters_ref, Rng& rng_ref,
                                                  const InteractionVertices& vertices,
                                                  const DMatrixBuilder<linalg::GPU>& builder_ref,
                                                  int id)
    : BaseClass(parameters_ref, rng_ref, vertices, builder_ref, id),
      configuration_copy_{configuration_, configuration_} {
  prepareInsertionMatrices(stream_index_);
  while (parameters_.getInitialConfigurationSize() > configuration_.size())
    tryVertexInsert(true);
}

template <class Parameters>
void CtintWalker<linalg::GPU, Parameters>::pruneInsertionMatrices(const int mat_idx, const int s) {
  details::removeIndices(Qa_[mat_idx][s], Ra_[mat_idx][s], removed_indices_[s]);
  removed_indices_[s].clear();
}

template <class Parameters>
void CtintWalker<linalg::GPU, Parameters>::prepareInsertionMatrices(const int stream) {
  auto& candidate_config = configuration_copy_[stream];
  auto& running_config = configuration_copy_[not stream];

  candidate_config = configuration_;
  // Reinsert the vertex whose insertion is pending.
  if (running_config.size() > 0)
    candidate_config.push_back(running_config.back());
  candidate_config.insertRandom(rng_);

  // Compute the new pieces of the D(= M^-1) matrix.
  d_builder_.buildSQR(Sa_[stream], Qa_[stream], Ra_[stream], devspace, candidate_config, thread_id_,
                      stream);
}

template <class Parameters>
bool CtintWalker<linalg::GPU, Parameters>::tryVertexInsert(bool forced) {
  auto& candidate_config = configuration_copy_[stream_index_];
  // Prepare the matrix for the next insertion.
  prepareInsertionMatrices(not stream_index_);

  const int delta = candidate_config.lastInsertionSize();
  if (delta == 1)
    configuration_.push_back(candidate_config.back());
  else {
    configuration_.push_back(candidate_config[candidate_config.size() - 2]);
    configuration_.push_back(candidate_config.back());
  }

  const std::array<ushort, 2> last_removal_size = {(ushort)removed_indices_[0].size(),
                                                   (ushort)removed_indices_[1].size()};
  for (int s = 0; s < 2; ++s) {
    linalg::util::syncStream(thread_id_, 2 * stream_index_ + 1);
    pruneInsertionMatrices(stream_index_, s);
  }

  const double accept_prob = BaseClass::insertionProbability(Sa_[stream_index_], Qa_[stream_index_],
                                                             Ra_[stream_index_], delta);
  const bool accept = std::abs(accept_prob) > rng_() or (forced and accept_prob != 0);

  if (accept) {
    if (accept_prob < 0)
      sign_ *= -1;
    BaseClass::applyInsertion(Sa_[stream_index_], Qa_[stream_index_], Ra_[stream_index_]);
  }
  else {
    configuration_.pop(delta);

    const std::array<int, 2> refusal_size = candidate_config.sizeIncrease();
    for (int s = 0; s < 2; ++s) {
      const int sector_size = candidate_config.getSector(s).size();
      for (int idx = sector_size - refusal_size[s] - last_removal_size[s];
           idx < sector_size - last_removal_size[s]; ++idx)
        removed_indices_[s].push_back(idx);
    }
  }

  // Flip the stream_index for the next step.
  stream_index_ = not stream_index_;
  return accept;
}

template <class Parameters>
bool CtintWalker<linalg::GPU, Parameters>::tryVertexRemoval() {
  const double accept_prob = BaseClass::removalProbability();
  const bool accept = std::abs(accept_prob) > rng_();

  if (accept) {
    if (accept_prob < 0)
      sign_ *= -1;

    for (int s = 0; s < 2; ++s)
      removed_indices_[s].insert(removed_indices_[s].end(),
                                 BaseClass::removal_matrix_indices_[s].begin(),
                                 BaseClass::removal_matrix_indices_[s].end());
    BaseClass::applyRemoval();
  }
  return accept;
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_HPP
