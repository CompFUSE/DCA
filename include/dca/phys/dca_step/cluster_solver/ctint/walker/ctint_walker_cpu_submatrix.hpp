// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Jérémie Bouquet (bouquetj@gmail.com).
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.com).

// This class organizes the MC walker in the CT-INT QMC purely on the CPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_base.hpp"

#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <ctime>

#include "dca/linalg/linalg.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/move.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/walker_methods.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {

template <class Parameters, typename Real, bool fix_rng_order>
class CtintWalkerSubmatrixCpu : public CtintWalkerBase<Parameters, Real> {
public:
  using this_type = CtintWalkerSubmatrixCpu;
  using BaseClass = CtintWalkerBase<Parameters, Real>;
  using typename BaseClass::Rng;
  using typename BaseClass::Data;
  using typename BaseClass::Profiler;
  using typename BaseClass::CudaStream;

  CtintWalkerSubmatrixCpu(const Parameters& pars_ref, const Data& /*data*/, Rng& rng_ref, int id = 0);

  virtual ~CtintWalkerSubmatrixCpu() = default;

  virtual void doSweep();

  virtual void computeM(typename BaseClass::MatrixPair& m_accum);

  using BaseClass::order;

protected:
  virtual void setMFromConfig() override;

  void doSteps();

  void generateDelayedMoves(int nbr_of_movesto_delay);

  void mainSubmatrixProcess();

  void updateM();

  void transformM();

  // For testing purposes.
  void doStep(const int nbr_of_movesto_delay);

private:
  virtual void doStep();

  void doSubmatrixUpdate();

  auto computeAcceptanceProbability();

  void updateGammaInv(int s);

  void removeRowAndColOfGammaInv();

  void computeMInit();

  //  void computeG0Init();
  void computeGInit();

  Move generateMoveType();

  void computeInsertionMatrices(const std::vector<int>& indices, const int s);

  void computeRemovalMatrix(int s);

  void computeMixedInsertionAndRemoval(int s);

  void findSectorIndices(const int s);

  void recomputeGammaInv();

  bool recentlyAdded(int move_idx, int s) {
    assert(move_idx < sector_indices_[s].size());
    return sector_indices_[s][move_idx] >= n_init_[s];
  }

protected:
  using MatrixView = linalg::MatrixView<Real, linalg::CPU>;
  using Matrix = linalg::Matrix<Real, linalg::CPU>;

  using BaseClass::parameters_;
  using BaseClass::configuration_;
  using BaseClass::rng_;
  using BaseClass::thread_id_;
  using BaseClass::d_builder_ptr_;
  using BaseClass::total_interaction_;
  using BaseClass::sign_;
  using BaseClass::M_;
  using BaseClass::n_bands_;
  using BaseClass::beta_;

  using BaseClass::thermalized_;

protected:
  int max_submatrix_size_;

  struct DelayedMoveType {
    Move move_type;
    std::array<double, 3> removal_rng{1., 1., 1.};
    Real acceptance_rng;
    std::array<int, 2> indices{-1, -1};
  };

protected:
  using BaseClass::acceptance_prob_;

protected:
  std::vector<DelayedMoveType> delayed_moves_;

  using MatrixPair = std::array<linalg::Matrix<Real, linalg::CPU>, 2>;
  MatrixPair G_;
  MatrixPair G0_;
  MatrixPair Gamma_inv_;      // TODO: don't pin
  MatrixPair Gamma_inv_cpy_;  // TODO: don't pin
  MatrixPair q_;              // TODO: don't pin or store in Gamma inv
  MatrixPair r_;
  MatrixPair s_;
  std::array<std::vector<Real>, 2> gamma_;

  Real det_ratio_;
  std::map<int, std::array<Real, n_bands_>> f_;
  std::map<int, std::array<Real, n_bands_>> prob_const_;
  std::map<std::pair<int, int>, std::array<Real, n_bands_>> gamma_values_;

  using BaseClass::nb_steps_per_sweep_;
  int nbr_of_steps_;
  int nbr_of_moves_to_delay_;
  int max_nbr_of_moves;

  std::array<int, 2> Gamma_size_;
  std::array<std::vector<int>, 2> Gamma_indices_;
  std::array<std::vector<int>, 2> sector_indices_;

  const DelayedMoveType* current_move_;

  std::array<unsigned, 2> nbr_of_indices_;

  std::vector<int> index_;

  // Initial configuration size.
  unsigned config_size_init_;

  // Initial and current sector sizes.
  std::array<unsigned, 2> n_init_;
  unsigned n_;

  // Maximal sector size after submatrix update.

  std::array<unsigned, 2> n_max_;

  std::array<linalg::util::HostVector<int>, 2> move_indices_;
  std::array<linalg::util::HostVector<int>, 2> removal_list_;
  std::array<linalg::util::HostVector<int>, 2> source_list_;
  std::array<std::vector<int>, 2> insertion_list_;
  std::array<std::vector<int>, 2> insertion_Gamma_indices_;

  std::vector<int> conf_removal_list_;

  std::vector<int>::iterator insertion_list_it_;

  std::array<Matrix, 2> Gamma_q_;
  Matrix workspace_;
  Matrix D_;

  using BaseClass::flop_;
};

template <class Parameters, typename Real, bool fix_rng_order>
CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::CtintWalkerSubmatrixCpu(
    const Parameters& parameters_ref, const Data& /*data*/, Rng& rng_ref, int id)
    : BaseClass(parameters_ref, rng_ref, id) {
  for (int b = 0; b < n_bands_; ++b) {
    for (int i = 1; i <= 3; ++i) {
      f_[i][b] = d_builder_ptr_->computeF(i, b);
      f_[-i][b] = d_builder_ptr_->computeF(-i, b);

      gamma_values_[std::make_pair(0, i)][b] = d_builder_ptr_->computeGamma(0, i, b);
      gamma_values_[std::make_pair(0, -i)][b] = d_builder_ptr_->computeGamma(0, -i, b);
      gamma_values_[std::make_pair(i, 0)][b] = d_builder_ptr_->computeGamma(i, 0, b);
      gamma_values_[std::make_pair(-i, 0)][b] = d_builder_ptr_->computeGamma(-i, 0, b);

      prob_const_[i][b] = prob_const_[-i][b] = -1. / (f_[i][b] - 1) / (f_[-i][b] - 1);
    }
    f_[0][b] = 1;
  }

  if (BaseClass::concurrency_.id() == BaseClass::concurrency_.first() && thread_id_ == 0)
    std::cout << "\nCT-INT submatrix walker created." << std::endl;
}

template <class Parameters, typename Real, bool fix_rng_order>

void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::setMFromConfig() {
  BaseClass::setMFromConfig();
  transformM();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::doSweep() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);
  doSteps();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::doSteps() {
  // Here nbr_of_steps is the number of single steps/moves during the current sweep,
  // while nbr_of_submatrix_steps is the number of times the entire submatrix algorithm  is run.

  if (nb_steps_per_sweep_ < 0)  // Not thermalized or fixed.
    nbr_of_steps_ = BaseClass::avgOrder() + 1;
  else
    nbr_of_steps_ = nb_steps_per_sweep_;

  // Get the maximum of Monte Carlo steps/moves that can be performed during one submatrix step.
  max_nbr_of_moves = parameters_.getMaxSubmatrixSize();

  BaseClass::n_steps_ += nbr_of_steps_;
  while (nbr_of_steps_ > 0) {
    nbr_of_moves_to_delay_ = std::min(nbr_of_steps_, max_nbr_of_moves);
    nbr_of_steps_ -= nbr_of_moves_to_delay_;

    doStep();
  }

  BaseClass::updateSweepAverages();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::doStep() {
  generateDelayedMoves(nbr_of_moves_to_delay_);
  doSubmatrixUpdate();
}

// Do one step with arbitrary number of moves. For testing.
template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::doStep(const int nbr_of_movesto_delay) {
  std::cout << "\nStarted doStep() function for testing." << std::endl;

  generateDelayedMoves(nbr_of_movesto_delay);

  std::cout << "\nGenerated " << nbr_of_movesto_delay << " moves for testing.\n" << std::endl;

  doSubmatrixUpdate();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::generateDelayedMoves(
    const int nbr_of_movesto_delay) {
  assert(nbr_of_movesto_delay > 0);

  delayed_moves_.clear();

  for (int s = 0; s < 2; ++s) {
    n_init_[s] = configuration_.getSector(s).size();
  }
  n_ = config_size_init_ = configuration_.size();

  // Generate delayed moves.
  for (int i = 0; i < nbr_of_movesto_delay; ++i) {
    DelayedMoveType delayed_move;

    delayed_move.move_type = generateMoveType();

    switch (delayed_move.move_type) {
      case REMOVAL:
        if constexpr (fix_rng_order) {
          if (configuration_.getDoubleUpdateProb())
            delayed_move.removal_rng = {rng_(), rng_(), rng_()};
          else
            delayed_move.removal_rng[0] = rng_();
        }
        break;

      case INSERTION:
        configuration_.insertRandom(rng_);
        for (int i = 0; i < configuration_.lastInsertionSize(); ++i)
          delayed_move.indices[i] = configuration_.size() - configuration_.lastInsertionSize() + i;
        break;

      default:
        throw(std::logic_error("Unknown move type encountered."));
    }

    delayed_move.acceptance_rng = rng_();
    delayed_moves_.push_back(delayed_move);
  }

  for (int s = 0; s < 2; ++s) {
    n_max_[s] = configuration_.getSector(s).size();
  }

  const int config_size_final = configuration_.size();
  conf_removal_list_.resize(config_size_final - config_size_init_);
  std::iota(conf_removal_list_.begin(), conf_removal_list_.end(), config_size_init_);
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::doSubmatrixUpdate() {
  computeMInit();
  computeGInit();
  mainSubmatrixProcess();
  updateM();
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::mainSubmatrixProcess() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    Gamma_inv_[s].resizeNoCopy(0);
    gamma_[s].clear();

    move_indices_[s].clear();
    insertion_list_[s].clear();
    insertion_Gamma_indices_[s].clear();
  }

  std::vector<int> aux_spin_type, new_aux_spin_type, move_band;

  for (int delay_ind = 0; delay_ind < delayed_moves_.size(); ++delay_ind) {
    current_move_ = &delayed_moves_[delay_ind];
    const auto move_type = current_move_->move_type;

    det_ratio_ = 1.;
    std::array<int, 2> indices_array;

    if (move_type == INSERTION) {
      indices_array = current_move_->indices;
    }
    else {  // move_type == REMOVAL
      if constexpr (fix_rng_order)
        indices_array =
            configuration_.randomRemovalCandidateSlow(delayed_moves_[delay_ind].removal_rng);
      else
        indices_array = configuration_.randomRemovalCandidate(rng_);
    }

    index_.clear();
    for (auto idx : indices_array) {
      if (idx >= 0)
        index_.push_back(idx);
    }

    bool at_least_one_recently_added = false;
    bool all_recently_added = false;
    if (move_type == REMOVAL) {
      // Check if the vertex to remove was inserted during the current submatrix update.
      const auto recently_added = [=](int id) { return id >= config_size_init_; };
      all_recently_added = math::util::all(index_, recently_added);
      at_least_one_recently_added = math::util::any(index_, recently_added);
    }

    for (int s = 0; s < 2; ++s) {
      Gamma_indices_[s].clear();
      new_aux_spin_type.clear();
      aux_spin_type.clear();
      move_band.clear();

      findSectorIndices(s);

      if (move_type == INSERTION) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          new_aux_spin_type.push_back(
              configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));
          aux_spin_type.push_back(0);
        }
      }
      else if (move_type == REMOVAL) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          new_aux_spin_type.push_back(0);
          aux_spin_type.push_back(
              configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));

          if (recentlyAdded(ind, s)) {
            // Find pre-existing Gamma_indices_.
            insertion_list_it_ = std::find(insertion_list_[s].begin(), insertion_list_[s].end(),
                                           sector_indices_[s][ind]);
            Gamma_indices_[s].push_back(insertion_Gamma_indices_[s][std::distance(
                insertion_list_[s].begin(), insertion_list_it_)]);
          }
        }
      }  // endif removal
      for (int index : sector_indices_[s])
        move_band.push_back(configuration_.getSector(s).getLeftB(index));

      for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
        if (move_type == INSERTION || !recentlyAdded(ind, s))
          gamma_[s].push_back(
              gamma_values_[std::make_pair(aux_spin_type[ind], new_aux_spin_type[ind])][move_band[ind]]);
        else {
          // TODO: find instead of adding.
          gamma_[s].push_back(
              gamma_values_[std::make_pair(new_aux_spin_type[ind], aux_spin_type[ind])][move_band[ind]]);
        }
      }

      if (!at_least_one_recently_added)
        computeInsertionMatrices(sector_indices_[s], s);
      else {
        if (all_recently_added)
          computeRemovalMatrix(s);
        else
          computeMixedInsertionAndRemoval(s);
      }

      det_ratio_ *= details::smallDeterminant(s_[s]);
    }  // s loop.

    if (n_ == 0 && move_type == REMOVAL)
      continue;

    // Compute acceptance probability.
    auto [acceptance_prob, mc_weigth_ratio] = computeAcceptanceProbability();
    acceptance_prob_ = acceptance_prob;

    // Note: acceptance and rejection can be forced for testing with the appropriate "acceptance_rng".
    const bool accepted =
        delayed_moves_[delay_ind].acceptance_rng < std::min(std::abs(acceptance_prob_), Real(1.));

    // NB: recomputeGammaInv is just a inefficient alternative to updateGammaInv. Only for testing
    // or debbuging.
    // recomputeGammaInv();

    // Update
    if (accepted) {
      ++BaseClass::n_accepted_;

      BaseClass::mc_log_weight_ += std::log(std::abs(mc_weigth_ratio));

      if (acceptance_prob_ < 0)
        sign_ *= -1;

      // Update GammaInv if necessary.
      if (!at_least_one_recently_added)
        for (int s = 0; s < 2; ++s)
          updateGammaInv(s);
      else {
        removeRowAndColOfGammaInv();
        for (int s = 0; s < 2; ++s) {
          for (int ind = sector_indices_[s].size() - 1; ind >= 0; --ind) {
            if (!recentlyAdded(ind, s))
              continue;
            insertion_list_[s].erase(std::remove(insertion_list_[s].begin(),
                                                 insertion_list_[s].end(), sector_indices_[s][ind]),
                                     insertion_list_[s].end());
          }
          for (int ind = Gamma_indices_[s].size() - 1; ind >= 0; --ind) {
            insertion_Gamma_indices_[s].erase(std::remove(insertion_Gamma_indices_[s].begin(),
                                                          insertion_Gamma_indices_[s].end(),
                                                          Gamma_indices_[s][ind]),
                                              insertion_Gamma_indices_[s].end());
            gamma_[s].erase(gamma_[s].begin() + Gamma_indices_[s][ind]);

            for (int i = 0; i < insertion_Gamma_indices_[s].size(); ++i) {
              if (insertion_Gamma_indices_[s][i] > Gamma_indices_[s][ind])
                --insertion_Gamma_indices_[s][i];
            }
          }
        }
      }

      for (int s = 0; s < 2; ++s) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          if (move_type == INSERTION || !recentlyAdded(ind, s)) {
            move_indices_[s].push_back(sector_indices_[s][ind]);
          }
          else {
            gamma_[s].pop_back();
            move_indices_[s].erase(std::remove(move_indices_[s].begin(), move_indices_[s].end(),
                                               sector_indices_[s][ind]),
                                   move_indices_[s].end());
          }
        }
      }

      if (move_type == INSERTION) {
        n_ += index_.size();
        for (auto idx : index_)
          configuration_.commitInsertion(idx);

        for (int s = 0; s < 2; ++s) {
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            insertion_list_[s].push_back(sector_indices_[s][ind]);
            insertion_Gamma_indices_[s].push_back(Gamma_inv_[s].nrRows() - nbr_of_indices_[s] + ind);
          }
        }

        // TODO: cleanup
        for (int idx : index_)
          conf_removal_list_.erase(
              std::remove(conf_removal_list_.begin(), conf_removal_list_.end(), idx),
              conf_removal_list_.end());
      }
      else if (move_type == REMOVAL) {
        n_ -= index_.size();
        for (auto idx : index_)
          configuration_.markForRemoval(idx);

        // TODO: cleanup.
        for (int idx : index_)
          conf_removal_list_.push_back(idx);
      }
    }
    else {  // The move is rejected:
      for (int s = 0; s < 2; ++s) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
          gamma_[s].pop_back();
      }
      if (at_least_one_recently_added && !all_recently_added)
        for (int s = 0; s < 2; ++s) {
          Gamma_inv_[s].swap(Gamma_inv_cpy_[s]);
        }
    }
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
Move CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::generateMoveType() {
  if (rng_() <= 0.5)
    return INSERTION;
  else
    return REMOVAL;
}

// Extend M by adding non-interacting vertices.
// M is not computed again and should be up-to-date.
template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::computeMInit() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];

    if (delta > 0) {
      Real f_j;
      D_.resize(std::make_pair(delta, n_init_[s]));

      d_builder_ptr_->computeG0(D_, configuration_.getSector(s), n_init_[s], n_max_[s], 0);

      for (int j = 0; j < n_init_[s]; ++j) {
        const auto field_type = configuration_.getSector(s).getAuxFieldType(j);
        const auto b = configuration_.getSector(s).getRightB(j);
        f_j = f_[field_type][b] - 1;

        for (int i = 0; i < delta; ++i) {
          D_(i, j) *= f_j;
        }
      }

      M_[s].resize(n_max_[s]);

      MatrixView M(M_[s], 0, 0, n_init_[s], n_init_[s]);
      MatrixView D_M(M_[s], n_init_[s], 0, delta, n_init_[s]);

      linalg::matrixop::gemm(D_, M, D_M);
      flop_ += 2 * D_.nrRows() * D_.nrCols() * M.nrCols();

      for (int i = 0; i < n_max_[s]; ++i) {
        for (int j = n_init_[s]; j < n_max_[s]; ++j) {
          M_[s](i, j) = 0;
        }
      }

      for (int i = n_init_[s]; i < n_max_[s]; ++i)
        M_[s](i, i) = 1;
    }
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::computeGInit() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];

    Real f;
    Matrix G0(std::make_pair(n_max_[s], delta));

    G_[s].resizeNoCopy(n_max_[s]);

    for (int j = 0; j < n_init_[s]; ++j) {
      const auto field_type = configuration_.getSector(s).getAuxFieldType(j);
      const auto b = configuration_.getSector(s).getRightB(j);
      f = f_[field_type][b];

      for (int i = 0; i < n_max_[s]; ++i) {
        G_[s](i, j) = (M_[s](i, j) * f - Real(i == j)) / (f - 1);
      }
    }

    if (delta > 0) {
      d_builder_ptr_->computeG0(G0, configuration_.getSector(s), n_init_[s], n_max_[s], 1);

      MatrixView G(G_[s], 0, n_init_[s], n_max_[s], delta);

      linalg::matrixop::gemm(M_[s], G0, G);
      flop_ += 2. * M_[s].nrRows() * M_[s].nrCols() * G0.nrCols();
    }
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
auto CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::computeAcceptanceProbability() {
  Real acceptance_probability = det_ratio_;

  Real gamma_factor = 1;
  const auto move_type = current_move_->move_type;
  for (int s = 0; s < 2; ++s) {
    if (!sector_indices_[s].size())
      continue;

    for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
      if (move_type == INSERTION || !recentlyAdded(ind, s))
        gamma_factor *= gamma_[s].at(gamma_[s].size() - nbr_of_indices_[s] + ind);
      else
        gamma_factor /= gamma_[s].at(gamma_[s].size() - nbr_of_indices_[s] + ind);
    }
  }
  acceptance_probability *= gamma_factor;

  const int delta_vertices = index_.size();
  const int interaction_sign = configuration_.getSign(index_[0]);
  const int non_empty_sector = sector_indices_[0].size() ? 0 : 1;

  Real mc_weight_ratio = acceptance_probability;
  Real K = total_interaction_;

  for (int v_id = 0; v_id < delta_vertices; ++v_id) {
    const auto field_type = configuration_.getSector(non_empty_sector)
                                .getAuxFieldType(sector_indices_[non_empty_sector][v_id]);
    const auto b =
        configuration_.getSector(non_empty_sector).getLeftB(sector_indices_[non_empty_sector][v_id]);
    K *= beta_ * prob_const_[field_type][b] * interaction_sign;

    const Real weight_term = prob_const_[field_type][b] * configuration_.getStrength(index_[v_id]);
    if (move_type == INSERTION)
      mc_weight_ratio *= weight_term;
    else
      mc_weight_ratio /= weight_term;
  }

  // Account for combinatorial factor and update acceptance probability.
  if (move_type == INSERTION) {
    if (delta_vertices == 1) {
      acceptance_probability *= K / (n_ + 1);
    }
    else if (delta_vertices == 2) {
      const Real possible_partners = configuration_.possiblePartners(index_[0]);
      const Real combinatorial_factor =
          (n_ + 2) * (configuration_.nPartners(index_[0]) + 1) / possible_partners;
      acceptance_probability *= configuration_.getStrength(index_[1]) * K / combinatorial_factor;
    }
    else
      throw(std::logic_error("Not implemented"));
  }
  else if (move_type == REMOVAL) {
    if (delta_vertices == 1) {
      acceptance_probability *= n_ / K;
    }
    else if (delta_vertices == 2) {
      const Real possible_partners = configuration_.possiblePartners(index_[0]);
      const Real combinatorial_factor = n_ * configuration_.nPartners(index_[0]) / possible_partners;
      acceptance_probability *= combinatorial_factor / (configuration_.getStrength(index_[1]) * K);
    }
    else
      throw(std::logic_error("Not implemented"));
  }

  return std::make_tuple(acceptance_probability, mc_weight_ratio);
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::updateGammaInv(int s) {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);
  const int delta = s_[s].nrRows();
  if (delta == 0)
    return;
  const int old_size = Gamma_inv_[s].nrRows();
  Gamma_inv_[s].resize(old_size + delta);

  if (old_size > 0) {
    MatrixView bulk(Gamma_inv_[s], 0, 0, old_size, old_size);
    MatrixView q_inv(Gamma_inv_[s], 0, old_size, old_size, delta);
    MatrixView r_inv(Gamma_inv_[s], old_size, 0, delta, old_size);
    MatrixView s_inv(Gamma_inv_[s], old_size, old_size, delta, delta);

    details::smallInverse(s_[s], s_inv);

    auto& Gamma_q = Gamma_q_[s];
    linalg::matrixop::gemm(Real(-1.), Gamma_q, s_inv, Real(0.), q_inv);

    auto& r_Gamma = workspace_;
    r_Gamma.resizeNoCopy(r_[s].size());
    linalg::matrixop::gemm(r_[s], bulk, r_Gamma);
    linalg::matrixop::gemm(Real(-1.), s_inv, r_Gamma, Real(0.), r_inv);

    // Gamma_ += Gamma_ * q_ * s_^-1 * r_ * Gamma_
    linalg::matrixop::gemm(Real(-1.), q_inv, r_Gamma, Real(1.), bulk);
  }
  else {
    Gamma_inv_[s].resizeNoCopy(delta);
    details::smallInverse(s_[s], Gamma_inv_[s]);
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::updateM() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    if (Gamma_inv_[s].nrRows() > 0) {
      int p;

      Matrix old_G(std::make_pair(n_max_[s], Gamma_inv_[s].nrRows()));
      Matrix old_M(std::make_pair(Gamma_inv_[s].nrRows(), n_max_[s]));
      Matrix result_matrix(std::make_pair(Gamma_inv_[s].nrRows(), n_max_[s]));

      for (int j = 0; j < Gamma_inv_[s].nrRows(); ++j) {
        p = move_indices_[s][j];

        for (int i = 0; i < n_max_[s]; ++i) {
          old_G(i, j) = G_[s](i, p);
          old_M(j, i) = M_[s](p, i);
        }
      }

      linalg::matrixop::gemm(Gamma_inv_[s], old_M, result_matrix);
      linalg::matrixop::gemm(Real(-1.), old_G, result_matrix, Real(1.), M_[s]);
      flop_ += 2 * Gamma_inv_[s].nrRows() * Gamma_inv_[s].nrCols() * old_M.nrCols();
      flop_ += 2 * old_G.nrRows() * old_G.nrCols() * result_matrix.nrCols();

      p = 0;
      for (auto& i : move_indices_[s]) {
        for (int j = 0; j < n_max_[s]; ++j) {
          M_[s](i, j) /= 1 + gamma_[s][p];
        }
        ++p;
      }
    }
  }

  // Remove "non-interacting" rows and columns.
  configuration_.moveAndShrink(source_list_, removal_list_, conf_removal_list_);

  for (int s = 0; s < 2; ++s) {
    linalg::matrixop::copyRows(M_[s], source_list_[s], M_[s], removal_list_[s]);
    linalg::matrixop::copyCols(M_[s], source_list_[s], M_[s], removal_list_[s]);
    M_[s].resize(configuration_.getSector(s).size());
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::findSectorIndices(const int s) {
  sector_indices_[s].clear();
  for (auto index : index_) {
    configuration_.findIndices(sector_indices_[s], index, s);
  }

  if (index_.size() > 1)
    std::sort(sector_indices_[s].begin(), sector_indices_[s].end());
  nbr_of_indices_[s] = sector_indices_[s].size();
}

// Remove row and column of Gamma_inv with Woodbury's formula.
// Gamma <- Gamma - U.V => GammaInv <- GammaInv + GammaInv.U.(Id-V.GammaInv.U)^(-1).V.GammaInv.
template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::removeRowAndColOfGammaInv() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);
  for (int s = 0; s < 2; ++s) {
    const int delta = Gamma_indices_[s].size();
    if (delta == 0)
      continue;
    const int n_init = Gamma_inv_[s].nrRows();
    const int n = n_init - delta;

    if (n) {
      q_[s].resizeNoCopy(std::make_pair(n, delta));
      r_[s].resizeNoCopy(std::make_pair(delta, n));
      auto& q_s = workspace_;
      q_s.resizeNoCopy(q_[s].size());

      //  TODO: check if gamma indices are ordered.
      for (int j = 0; j < delta; ++j) {
        int idx = 0;
        for (int i = 0; i < n_init; ++i) {
          if (idx == delta || i != Gamma_indices_[s][idx])
            q_[s](i - idx, j) = Gamma_inv_[s](i, Gamma_indices_[s][j]);
          else
            ++idx;
        }
      }

      int idx = 0;
      for (int j = 0; j < n_init; ++j) {
        if (idx == delta || j != Gamma_indices_[s][idx])
          for (int i = 0; i < delta; ++i)
            r_[s](i, j - idx) = Gamma_inv_[s](Gamma_indices_[s][i], j);
        else
          ++idx;
      }

      // TODO: this could be avoided in the case of pure removal.
      s_[s].resizeNoCopy(delta);
      for (int j = 0; j < delta; ++j)
        for (int i = 0; i < delta; ++i)
          s_[s](i, j) = Gamma_inv_[s](Gamma_indices_[s][i], Gamma_indices_[s][j]);

      for (int ind = delta - 1; ind >= 0; --ind)
        linalg::matrixop::removeRowAndCol(Gamma_inv_[s], Gamma_indices_[s][ind]);

      details::smallInverse(s_[s]);

      linalg::matrixop::gemm(q_[s], s_[s], q_s);

      // Gamma_inv_ -= Q*S^-1*R
      linalg::matrixop::gemm(Real(-1.), q_s, r_[s], Real(1.), Gamma_inv_[s]);
    }  // if n
    else {
      Gamma_inv_[s].resizeNoCopy(0);
    }
  }
}

// This method is unused and left to potentially use as a testing reference.
template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::recomputeGammaInv() {
  for (int s = 0; s < 2; ++s) {
    if (Gamma_inv_[s].nrRows() > 0)
      linalg::matrixop::inverse(Gamma_inv_[s]);

    Gamma_inv_[s].resize(Gamma_inv_[s].nrRows() + nbr_of_indices_[s]);

    for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
      for (int i = 0; i < Gamma_inv_[s].nrRows(); ++i) {
        Gamma_inv_[s](Gamma_inv_[s].nrRows() + ind, i) =
            G_[s](sector_indices_[s][ind], move_indices_[s][i]);
        Gamma_inv_[s](i, Gamma_inv_[s].nrRows() + ind) =
            G_[s](move_indices_[s][i], sector_indices_[s][ind]);
      }

      for (int ind_2 = 0; ind_2 < nbr_of_indices_[s]; ++ind_2) {
        Gamma_inv_[s](Gamma_inv_[s].nrRows() + ind, Gamma_inv_[s].nrRows() + ind_2) =
            G_[s](sector_indices_[s][ind], sector_indices_[s][ind_2]);
      }
      Gamma_inv_[s](Gamma_inv_[s].nrRows() + ind, Gamma_inv_[s].nrRows() + ind) -=
          (1 + gamma_[s].end()[-nbr_of_indices_[s] + ind]) /
          gamma_[s].end()[-nbr_of_indices_[s] + ind];
    }

    Gamma_inv_[s].nrRows() = Gamma_inv_[s].nrRows();

    if (Gamma_inv_[s].nrRows() > 0)
      details::smallInverse(Gamma_inv_[s]);
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::transformM() {
  for (int s = 0; s < 2; ++s) {
    for (int j = 0; j < M_[s].size().second; ++j) {
      for (int i = 0; i < M_[s].size().first; ++i) {
        const auto field_type = configuration_.getSector(s).getAuxFieldType(i);
        const auto b = configuration_.getSector(s).getLeftB(i);
        const Real f_i = -(f_[field_type][b] - 1);
        M_[s](i, j) /= f_i;
      }
    }
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::computeM(
    typename BaseClass::MatrixPair& m_accum) {
  for (int s = 0; s < 2; ++s) {
    m_accum[s].resizeNoCopy(M_[s].size());

    for (int j = 0; j < M_[s].size().second; ++j) {
      for (int i = 0; i < M_[s].size().first; ++i) {
        const auto field_type = configuration_.getSector(s).getAuxFieldType(i);
        const auto b = configuration_.getSector(s).getLeftB(i);
        const Real factor = -(f_[field_type][b] - 1.);
        m_accum[s](i, j) = M_[s](i, j) * factor;
      }
    }
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::computeInsertionMatrices(
    const std::vector<int>& insertion_indices, const int s) {
  const int delta = insertion_indices.size();

  s_[s].resizeNoCopy(delta);
  if (delta == 0)
    return;

  for (int i = 0; i < delta; ++i)
    for (int j = 0; j < delta; ++j) {
      s_[s](i, j) = G_[s](insertion_indices[i], insertion_indices[j]);
      if (i == j) {
        const Real gamma_val = gamma_[s].at(gamma_[s].size() + i - nbr_of_indices_[s]);
        s_[s](i, j) -= (1 + gamma_val) / gamma_val;
      }
    }

  assert(Gamma_inv_[s].nrRows() == move_indices_[s].size());
  if (Gamma_inv_[s].nrRows() > 0) {
    q_[s].resizeNoCopy(std::make_pair(Gamma_inv_[s].nrRows(), delta));
    r_[s].resizeNoCopy(std::make_pair(delta, Gamma_inv_[s].nrRows()));

    for (int i = 0; i < Gamma_inv_[s].nrRows(); ++i)
      for (int j = 0; j < delta; ++j) {
        q_[s](i, j) = G_[s](move_indices_[s][i], insertion_indices[j]);
        r_[s](j, i) = G_[s](insertion_indices[j], move_indices_[s][i]);
      }

    auto& Gamma_q = Gamma_q_[s];
    Gamma_q.resizeNoCopy(q_[s].size());
    linalg::matrixop::gemm(Gamma_inv_[s], q_[s], Gamma_q);
    linalg::matrixop::gemm(Real(-1.), r_[s], Gamma_q, Real(1.), s_[s]);
  }
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::computeRemovalMatrix(const int s) {
  const int delta = Gamma_indices_[s].size();
  s_[s].resizeNoCopy(delta);
  for (int j = 0; j < delta; ++j)
    for (int i = 0; i < delta; ++i)
      s_[s](i, j) = Gamma_inv_[s](Gamma_indices_[s][i], Gamma_indices_[s][j]);
}

template <class Parameters, typename Real, bool fix_rng_order>
void CtintWalkerSubmatrixCpu<Parameters, Real, fix_rng_order>::computeMixedInsertionAndRemoval(int s) {
  Gamma_inv_cpy_[s] = Gamma_inv_[s];

  if (sector_indices_[s].size() == 0) {
    s_[s].resizeNoCopy(0);
    return;
  }

  // TODO: avoid reallocation.
  std::vector<int> insertion_indices;
  for (int i = 0; i < nbr_of_indices_[s]; ++i)
    if (!recentlyAdded(i, s))
      insertion_indices.push_back(sector_indices_[s][i]);
  assert(Gamma_indices_[s].size() + insertion_indices.size() == sector_indices_.size());

  computeInsertionMatrices(insertion_indices, s);
  det_ratio_ *= details::smallDeterminant(s_[s]);
  updateGammaInv(s);

  computeRemovalMatrix(s);
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP
