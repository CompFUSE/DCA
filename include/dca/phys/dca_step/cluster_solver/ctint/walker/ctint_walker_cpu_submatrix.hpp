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
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"
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

template <class Parameters>
class CtintWalkerSubmatrix<linalg::CPU, Parameters> : public CtintWalkerBase<Parameters> {
public:
  using this_type = CtintWalkerSubmatrix<linalg::CPU, Parameters>;
  using BaseClass = CtintWalkerBase<Parameters>;
  using Rng = typename BaseClass::Rng;
  using Profiler = typename Parameters::profiler_type;

  CtintWalkerSubmatrix(const Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
                       const DMatrixBuilder<linalg::CPU>& builder_ref, int id = 0);

  virtual ~CtintWalkerSubmatrix() = default;

  virtual void doSweep();

  using BaseClass::order;

protected:
  void doSteps();
  void generateDelayedMoves(int nbr_of_movesto_delay);
  void mainSubmatrixProcess();
  void updateM();

  // For testing purposes.
  void doStep(const int nbr_of_movesto_delay);

private:
  virtual void doStep();
  void doSubmatrixUpdate();
  double computeAcceptanceProbability();
  void updateGammaInv();
  void removeRowAndColOfGammaInv(const int s);
  void pushToEnd();
  void computeMInit();
  //  void computeG0Init();
  void computeGInit();
  Move generateMoveType();

  void findSectorIndices(const int s);
  void recomputeGammaInv();

protected:
  using MatrixView = linalg::MatrixView<double, linalg::CPU>;
  using Matrix = linalg::Matrix<double, linalg::CPU>;

  using BaseClass::parameters_;
  using BaseClass::configuration_;
  using BaseClass::rng_;
  using BaseClass::thread_id_;
  using BaseClass::d_builder_;
  using BaseClass::total_interaction_;
  using BaseClass::sign_;
  using BaseClass::M_;
  using BaseClass::n_bands_;

  double max_tau_ = BaseClass::beta_;

  using BaseClass::thermalized_;

protected:
  int max_submatrix_size_;

  struct DelayedMoveType {
    Move move_type_;

    double removal_rng_;
    double acceptance_rng_;

    int index_;
    //    int partner_index_;
    std::array<std::vector<int>, 2> sector_indices_;
    //    std::array<std::vector<int>, 2> partner_sector_indices_;
  };

protected:
  using BaseClass::acceptance_prob_;

protected:
  std::vector<DelayedMoveType> delayed_moves_;

  using MatrixPair = std::array<linalg::Matrix<double, linalg::CPU>, 2>;
  MatrixPair G_;
  MatrixPair G0_;
  MatrixPair Gamma_inv_;
  MatrixPair q_;
  MatrixPair r_;
  MatrixPair s_;
  std::array<std::vector<double>, 2> gamma_;

  std::map<int, double> f_;
  std::map<int, double> prob_const_;
  std::map<std::pair<int, int>, double> gamma_values_;

  using BaseClass::nb_steps_per_sweep_;
  int nbr_of_steps_;
  int nbr_of_submatrix_steps_;
  int nbr_of_moves_to_delay_;
  int max_nbr_of_moves;

  std::array<int, 2> Gamma_size_;
  std::array<std::vector<int>, 2> Gamma_indices_;
  std::array<std::vector<int>, 2> sector_indices_;
  Move move_type_;

  std::array<int, 2> nbr_of_indices_;
  std::array<bool, 2> found_sector_index_;

  int index_;
  double tau_;

  // Initial configuration size.
  int config_size_init_;

  // Initial and current sector sizes.
  std::array<int, 2> n_init_;
  std::array<int, 2> n_;

  // Maximal sector size after submatrix update.

  std::array<int, 2> n_max_;

  std::array<std::vector<int>, 2> move_indices_;
  std::vector<int> existing_indices_;
  std::array<std::vector<int>, 2> removal_list_;
  std::array<std::vector<int>, 2> insertion_list_;
  std::array<std::vector<int>, 2> insertion_Gamma_indices_;

  std::vector<int> conf_removal_list_;

  std::vector<int>::iterator insertion_list_it_;
  bool recently_added_;
  bool accepted_;

  std::array<Matrix, 2> Gamma_q_;
  Matrix workspace_;
  Matrix D_;
};

template <class Parameters>
CtintWalkerSubmatrix<linalg::CPU, Parameters>::CtintWalkerSubmatrix(
    const Parameters& parameters_ref, Rng& rng_ref, const InteractionVertices& vertices,
    const DMatrixBuilder<linalg::CPU>& builder_ref, int id)
    : BaseClass(parameters_ref, rng_ref, vertices, builder_ref, id) {
  for (int i = 1; i <= 3; ++i) {
    f_[i] = d_builder_.computeF(i);
    f_[-i] = d_builder_.computeF(-i);

    gamma_values_[std::make_pair(0, i)] = d_builder_.computeGamma(0, i);
    gamma_values_[std::make_pair(0, -i)] = d_builder_.computeGamma(0, -i);
    gamma_values_[std::make_pair(i, 0)] = d_builder_.computeGamma(i, 0);
    gamma_values_[std::make_pair(-i, 0)] = d_builder_.computeGamma(-i, 0);

    prob_const_[i] = prob_const_[-i] = -max_tau_ / (f_[i] - 1) / (f_[-i] - 1);
  }
  f_[0] = 1;

  if (BaseClass::concurrency_.id() == BaseClass::concurrency_.first() && thread_id_ == 0)
    std::cout << "\nCT-INT submatrix walker created." << std::endl;
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::doSweep() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);
  double f_i;

  for (int s = 0; s < 2; ++s) {
    for (int i = 0; i < M_[s].size().first; ++i) {
      f_i = f_[configuration_.getSector(s).getAuxFieldType(i)];

      for (int j = 0; j < M_[s].size().second; ++j) {
        M_[s](i, j) /= -(f_i - 1);
      }
    }
  }

  doSteps();

  for (int s = 0; s < 2; ++s) {
    for (int i = 0; i < M_[s].size().first; ++i) {
      f_i = f_[configuration_.getSector(s).getAuxFieldType(i)];

      for (int j = 0; j < M_[s].size().second; ++j) {
        M_[s](i, j) *= -(f_i - 1);
      }
    }
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::doSteps() {
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

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::doStep() {
  generateDelayedMoves(nbr_of_moves_to_delay_);
  doSubmatrixUpdate();
}

// Do one step with arbitrary number of moves. For testing.
template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::doStep(const int nbr_of_movesto_delay) {
  std::cout << "\nStarted doStep() function for testing." << std::endl;

  double f_i;

  for (int s = 0; s < 2; ++s) {
    for (int i = 0; i < M_[s].size().first; ++i) {
      f_i = f_[configuration_.getSector(s).getAuxFieldType(i)];

      for (int j = 0; j < M_[s].size().second; ++j) {
        M_[s](i, j) /= -(f_i - 1);
      }
    }
  }

  generateDelayedMoves(nbr_of_movesto_delay);

  std::cout << "\nGenerated " << nbr_of_movesto_delay << " moves for testing.\n" << std::endl;

  doSubmatrixUpdate();

  for (int s = 0; s < 2; ++s) {
    for (int i = 0; i < M_[s].size().first; ++i) {
      f_i = f_[configuration_.getSector(s).getAuxFieldType(i)];

      for (int j = 0; j < M_[s].size().second; ++j) {
        M_[s](i, j) *= -(f_i - 1);
      }
    }
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::generateDelayedMoves(const int nbr_of_movesto_delay) {
  assert(nbr_of_movesto_delay > 0);

  delayed_moves_.clear();

  int nbr_of_moves = 0;

  for (int s = 0; s < 2; ++s) {
    n_init_[s] = configuration_.getSector(s).size();
    n_[s] = n_init_[s];
  }

  config_size_init_ = configuration_.size();
  assert(config_size_init_ == (n_init_[0] + n_init_[1]) / 2);

  // Generate delayed moves.

  while (nbr_of_moves < nbr_of_movesto_delay) {
    DelayedMoveType delayed_move;

    delayed_move.move_type_ = generateMoveType();

    switch (delayed_move.move_type_) {
      case REMOVAL:
        delayed_move.removal_rng_ = rng_();
        break;

      case INSERTION:
        configuration_.insertRandom(rng_);

        delayed_move.index_ = configuration_.size() - 1;

        for (int s = 0; s < 2; ++s) {
          delayed_move.sector_indices_[s].push_back(configuration_.getSector(s).size() - 1);
        }

        // Depending on the parameters, insertRandom() may sometimes insert two vertices at once.
        if (configuration_.lastInsertionSize() > 1) {
          throw(std::logic_error("\nDouble insertion not yet supported."));
        }
        break;

      default:
        throw(std::logic_error("Unkown move type encountered."));
    }
    delayed_move.acceptance_rng_ = rng_();
    delayed_moves_.push_back(delayed_move);

    ++nbr_of_moves;
  }

  existing_indices_.clear();

  existing_indices_.resize(config_size_init_);
  std::iota(existing_indices_.begin(), existing_indices_.end(), 0);

  for (int s = 0; s < 2; ++s) {
    n_max_[s] = configuration_.getSector(s).size();

    removal_list_[s].resize(n_max_[s] - n_init_[s]);
    std::iota(removal_list_[s].begin(), removal_list_[s].end(), n_init_[s]);
  }

  const int config_size_final = configuration_.size();
  conf_removal_list_.resize(config_size_final - config_size_init_);
  std::iota(conf_removal_list_.begin(), conf_removal_list_.end(), config_size_init_);
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::doSubmatrixUpdate() {
  computeMInit();
  computeGInit();
  mainSubmatrixProcess();
  updateM();
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::mainSubmatrixProcess() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    Gamma_size_[s] = 0;
    Gamma_inv_[s].resizeNoCopy(0);
    gamma_[s].clear();

    move_indices_[s].clear();
    insertion_list_[s].clear();
    insertion_Gamma_indices_[s].clear();
  }

  std::vector<int> aux_spin_type, new_aux_spin_type;
  bool do_nothing;

  for (int delay_ind = 0; delay_ind < delayed_moves_.size(); ++delay_ind) {
    move_type_ = delayed_moves_[delay_ind].move_type_;
    do_nothing = false;
    recently_added_ = false;

    if (move_type_ == INSERTION)
      index_ = delayed_moves_[delay_ind].index_;
    else {  // move_type == REMOVAL
      // Do nothing if there is no vertex to remove.
      if (existing_indices_.size() == 0)
        do_nothing = true;
      else
        index_ =
            existing_indices_[int(delayed_moves_[delay_ind].removal_rng_ * existing_indices_.size())];
      // Check if the vertex to remove was inserted during the current submatrix update.
      if (index_ >= config_size_init_)
        recently_added_ = true;
    }

    for (int s = 0; s < 2; ++s) {
      assert(Gamma_size_[s] == Gamma_inv_[s].nrRows());

      Gamma_indices_[s].clear();
      new_aux_spin_type.clear();
      aux_spin_type.clear();

      if (do_nothing)
        continue;

      findSectorIndices(s);
      // Continue to next sector/move if there is nothing to change in current sector.
      if (!found_sector_index_[s])
        continue;

      if (move_type_ == INSERTION) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          new_aux_spin_type.push_back(
              configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));
          aux_spin_type.push_back(0);
        }
      }
      else if (move_type_ == REMOVAL) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          new_aux_spin_type.push_back(0);
          aux_spin_type.push_back(
              configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));
        }

        if (recently_added_) {
          // Compute Gamma_indices_.
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            insertion_list_it_ = std::find(insertion_list_[s].begin(), insertion_list_[s].end(),
                                           sector_indices_[s][ind]);
            Gamma_indices_[s].push_back(insertion_Gamma_indices_[s][std::distance(
                insertion_list_[s].begin(), insertion_list_it_)]);
          }
        }
      }  // else if(removal).

      //(What follows is executed for both insertions or removals).

      for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
        if (!recently_added_)
          gamma_[s].push_back(
              gamma_values_[std::make_pair(aux_spin_type[ind], new_aux_spin_type[ind])]);
        else {
          // TODO: do not add.
          gamma_[s].push_back(
              gamma_values_[std::make_pair(new_aux_spin_type[ind], aux_spin_type[ind])]);
        }
      }

      if (!recently_added_) {
        const int delta = nbr_of_indices_[s];
        s_[s].resizeNoCopy(delta);
        for (int i = 0; i < delta; ++i)
          for (int j = 0; j < delta; ++j) {
            s_[s](i, j) = G_[s](sector_indices_[s][i], sector_indices_[s][j]);
            if (i == j) {
              const double gamma_val = gamma_[s].end()[i - nbr_of_indices_[s]];
              s_[s](i, j) -= (1 + gamma_val) / gamma_val;
            }
          }

        assert(Gamma_size_[s] == move_indices_[s].size());
        if (Gamma_size_[s] > 0) {
          q_[s].resizeNoCopy(std::make_pair(Gamma_size_[s], delta));
          r_[s].resizeNoCopy(std::make_pair(delta, Gamma_size_[s]));

          for (int i = 0; i < Gamma_size_[s]; ++i)
            for (int j = 0; j < delta; ++j) {
              q_[s](i, j) = G_[s](move_indices_[s][i], sector_indices_[s][j]);
              r_[s](j, i) = G_[s](sector_indices_[s][j], move_indices_[s][i]);
            }

          auto& Gamma_q = Gamma_q_[s];
          Gamma_q.resizeNoCopy(q_[s].size());
          linalg::matrixop::gemm(Gamma_inv_[s], q_[s], Gamma_q);
          linalg::matrixop::gemm(-1., r_[s], Gamma_q, 1., s_[s]);
        }
      }
      else {
        const int delta = Gamma_indices_[s].size();
        s_[s].resizeNoCopy(delta);
        for (int j = 0; j < delta; ++j)
          for (int i = 0; i < delta; ++i)
            s_[s](i, j) = Gamma_inv_[s](Gamma_indices_[s][i], Gamma_indices_[s][j]);
      }
    }  // s loop.

    // Directly go to next move if do_nothing is true.

    if (do_nothing)
      continue;

    // Compute acceptance probability.
    acceptance_prob_ = computeAcceptanceProbability();

    // Acceptance can be forced (for testing).
    accepted_ = delayed_moves_[delay_ind].acceptance_rng_ < std::min(std::abs(acceptance_prob_), 1.);

    // NB: recomputeGammaInv is just a inefficient alternative to updateGammaInv. Only for testing
    // or debbuging.
    //
    // recomputeGammaInv();

    // Update other objects.

    if (accepted_) {
      ++BaseClass::n_accepted_;
      if (acceptance_prob_ < 0)
        sign_ *= -1;

      // Update GammaInv if necessary.
      if (!recently_added_) {
        updateGammaInv();
      }
      else {
        for (int s = 0; s < 2; ++s) {
          removeRowAndColOfGammaInv(s);

          // TODO: write a more general method for indices removal.
          // TODO: if Gamma_inv_ is not physically shrunk, do not resize insertion_Gamma_indices_.
          for (int ind = nbr_of_indices_[s] - 1; ind >= 0; --ind) {
            insertion_list_[s].erase(std::remove(insertion_list_[s].begin(),
                                                 insertion_list_[s].end(), sector_indices_[s][ind]),
                                     insertion_list_[s].end());

            insertion_Gamma_indices_[s].erase(std::remove(insertion_Gamma_indices_[s].begin(),
                                                          insertion_Gamma_indices_[s].end(),
                                                          Gamma_indices_[s][ind]),
                                              insertion_Gamma_indices_[s].end());
            gamma_[s].erase(gamma_[s].begin() + Gamma_indices_[s][ind]);
            //            gamma_[s].erase(gamma_[s].begin() + Gamma_indices_[s][ind] - ind);

            for (int i = 0; i < insertion_Gamma_indices_[s].size(); ++i) {
              if (insertion_Gamma_indices_[s][i] > Gamma_indices_[s][ind])
                --insertion_Gamma_indices_[s][i];
            }
          }
        }
      }

      if (!recently_added_)
        for (int s = 0; s < 2; ++s)
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
            move_indices_[s].push_back(sector_indices_[s][ind]);
      else {
        for (int s = 0; s < 2; ++s)
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            gamma_[s].pop_back();

            move_indices_[s].erase(std::remove(move_indices_[s].begin(), move_indices_[s].end(),
                                               sector_indices_[s][ind]),
                                   move_indices_[s].end());
          }
      }

      if (move_type_ == INSERTION) {
        for (int s = 0; s < 2; ++s) {
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            ++n_[s];
            insertion_list_[s].push_back(sector_indices_[s][ind]);
            insertion_Gamma_indices_[s].push_back(Gamma_size_[s] - nbr_of_indices_[s] + ind);

            removal_list_[s].erase(std::remove(removal_list_[s].begin(), removal_list_[s].end(),
                                               sector_indices_[s][ind]),
                                   removal_list_[s].end());
          }
        }

        conf_removal_list_.erase(
            std::remove(conf_removal_list_.begin(), conf_removal_list_.end(), index_),
            conf_removal_list_.end());

        existing_indices_.push_back(index_);
      }
      else if (move_type_ == REMOVAL) {
        for (int s = 0; s < 2; ++s) {
          // TODO: append without loop.
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            --n_[s];
            removal_list_[s].push_back(sector_indices_[s][ind]);
          }
        }

        conf_removal_list_.push_back(index_);

        existing_indices_.erase(
            std::remove(existing_indices_.begin(), existing_indices_.end(), index_),
            existing_indices_.end());
      }
    }

    // If the move is rejected:
    else {
      for (int s = 0; s < 2; ++s) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
          // TODO: do not pop if not added.
          gamma_[s].pop_back();
      }
    }
  }
}

template <class Parameters>
Move CtintWalkerSubmatrix<linalg::CPU, Parameters>::generateMoveType() {
  if (rng_() <= 0.5)
    return INSERTION;
  else
    return REMOVAL;
}

// Extend M by adding non-interacting vertices.
// M is not computed again and should be up-to-date.
template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::computeMInit() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];

    if (delta > 0) {
      double f_j;
      D_.resize(std::make_pair(delta, n_init_[s]));

      d_builder_.computeG0(D_, configuration_.getSector(s), n_init_[s], n_max_[s], 0);

      for (int j = 0; j < n_init_[s]; ++j) {
        f_j = f_[configuration_.getSector(s).getAuxFieldType(j)] - 1;

        for (int i = 0; i < delta; ++i) {
          D_(i, j) *= f_j;
        }
      }

      M_[s].resize(n_max_[s]);

      MatrixView M(M_[s], 0, 0, n_init_[s], n_init_[s]);
      MatrixView D_M(M_[s], n_init_[s], 0, delta, n_init_[s]);

      linalg::matrixop::gemm(D_, M, D_M);

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

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::computeGInit() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];

    double f;
    Matrix G0(std::make_pair(n_max_[s], delta));

    G_[s].resizeNoCopy(n_max_[s]);

    for (int j = 0; j < n_init_[s]; ++j) {
      f = f_[configuration_.getSector(s).getAuxFieldType(j)];

      for (int i = 0; i < n_max_[s]; ++i) {
        G_[s](i, j) = (M_[s](i, j) * f - double(i == j)) / (f - 1);
      }
    }

    if (delta > 0) {
      d_builder_.computeG0(G0, configuration_.getSector(s), n_init_[s], n_max_[s], 1);

      MatrixView G(G_[s], 0, n_init_[s], n_max_[s], delta);

      linalg::matrixop::gemm(M_[s], G0, G);
    }
  }
}

// template <class Parameters>
// void CtintWalkerSubmatrix<linalg::CPU, Parameters>::computeG0Init() {
//  for (int s = 0; s < 2; ++s) {
//    d_builder_.computeG0Init(G0_[s], configuration_.getSector(s), n_init_[s], n_max_[s]);
//  }
//}

template <class Parameters>
double CtintWalkerSubmatrix<linalg::CPU, Parameters>::computeAcceptanceProbability() {
  double acceptance_probability = 1;
  double gamma_factor = 1;

  for (int s = 0; s < 2; ++s) {
    if (!sector_indices_[s].size())
      continue;

    acceptance_probability *= details::smallDeterminant(s_[s]);

    for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
      gamma_factor *= gamma_[s].end()[-nbr_of_indices_[s] + ind];
  }

  if (!recently_added_)
    acceptance_probability *= gamma_factor;
  else
    acceptance_probability /= gamma_factor;

  double K;
  if (sector_indices_[0].size() != 0)
    K = total_interaction_ *
        prob_const_.at(configuration_.getSector(0).getAuxFieldType(sector_indices_[0][0]));
  else
    K = total_interaction_ *
        prob_const_.at(configuration_.getSector(1).getAuxFieldType(sector_indices_[1][0]));

  int n = (n_[0] + n_[1]) / 2;  // TODO: use only one n.

  if (move_type_ == INSERTION)
    acceptance_probability *= K / (n + 1);
  else if (move_type_ == REMOVAL)
    acceptance_probability *= n / K;

  return acceptance_probability;
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::updateGammaInv() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);
  for (int s = 0; s < 2; ++s) {
    const int delta = nbr_of_indices_[s];
    if (delta == 0)
      continue;
    Gamma_inv_[s].resize(Gamma_size_[s] + delta);

    if (Gamma_size_[s] > 0) {
      MatrixView bulk(Gamma_inv_[s], 0, 0, Gamma_size_[s], Gamma_size_[s]);
      MatrixView q_inv(Gamma_inv_[s], 0, Gamma_size_[s], Gamma_size_[s], delta);
      MatrixView r_inv(Gamma_inv_[s], Gamma_size_[s], 0, delta, Gamma_size_[s]);
      MatrixView s_inv(Gamma_inv_[s], Gamma_size_[s], Gamma_size_[s], delta, delta);

      details::smallInverse(s_[s], s_inv);

      auto& Gamma_q = Gamma_q_[s];
      Gamma_q.resizeNoCopy(q_[s].size());
      linalg::matrixop::gemm(bulk, q_[s], Gamma_q);
      linalg::matrixop::gemm(-1., Gamma_q, s_inv, 0., q_inv);

      // TODO: reuse previous result.
      auto& r_Gamma = workspace_;
      r_Gamma.resizeNoCopy(r_[s].size());
      linalg::matrixop::gemm(r_[s], bulk, r_Gamma);
      linalg::matrixop::gemm(-1., s_inv, r_Gamma, 0., r_inv);

      // Gamma_ += Gamma_ * q_ * s_^-1 * r_ * Gamma_
      linalg::matrixop::gemm(-1., q_inv, r_Gamma, 1., bulk);
    }
    else {
      Gamma_inv_[s].resizeNoCopy(delta);
      details::smallInverse(s_[s], Gamma_inv_[s]);
    }
  }

  Gamma_size_[0] = Gamma_inv_[0].nrRows();
  Gamma_size_[1] = Gamma_inv_[1].nrRows();
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::updateM() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    if (Gamma_size_[s] > 0) {
      int p;

      Matrix old_G(std::make_pair(n_max_[s], Gamma_size_[s]));
      Matrix old_M(std::make_pair(Gamma_size_[s], n_max_[s]));
      Matrix result_matrix(std::make_pair(Gamma_size_[s], n_max_[s]));

      for (int j = 0; j < Gamma_size_[s]; ++j) {
        p = move_indices_[s][j];

        for (int i = 0; i < n_max_[s]; ++i) {
          old_G(i, j) = G_[s](i, p);
          old_M(j, i) = M_[s](p, i);
        }
      }

      linalg::matrixop::gemm(Gamma_inv_[s], old_M, result_matrix);
      linalg::matrixop::gemm(-1., old_G, result_matrix, 1., M_[s]);

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

  pushToEnd();

  for (int s = 0; s < 2; ++s) {
    M_[s].resize(n_max_[s] - removal_list_[s].size());
  }

  while (configuration_.size() > (M_[0].size().first + M_[1].size().first) / 2)
    configuration_.pop();
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::findSectorIndices(const int s) {
  // TODO: use configuration_.findIndices(index_, s);
  // TODO: use sector indices for insertion.
  tau_ = configuration_[index_].tau;

  found_sector_index_[s] = false;
  sector_indices_[s].clear();

  // NB: In multiband case, two sector indices may correspond to same tau (?).
  // Then the corresponding sector must be extended/shrunk by 2 rows and columns.

  for (int i = 0; i < configuration_.getSector(s).size(); ++i) {
    if (configuration_.getSector(s).getTau(i) == tau_) {
      found_sector_index_[s] = true;
      sector_indices_[s].push_back(i);
      if (sector_indices_[s].size() == 2)
        break;
    }
  }

  nbr_of_indices_[s] = sector_indices_[s].size();
}

// Remove row and column of Gamma_inv with Woodbury's formula.

// Gamma <- Gamma - U.V => GammaInv <- GammaInv + GammaInv.U.(Id-V.GammaInv.U)^(-1).V.GammaInv.

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::removeRowAndColOfGammaInv(const int s) {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);
  // TODO: optimize.
  const int delta = Gamma_indices_[s].size();
  if (delta == 0)
    return;
  const int n_init = Gamma_size_[s];
  const int n = n_init - delta;

  if (n) {
    // TODO: MAYBE do not resize Gamma but set sector to id.
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

    for (int ind = nbr_of_indices_[s] - 1; ind >= 0; --ind)
      linalg::matrixop::removeRowAndCol(Gamma_inv_[s], Gamma_indices_[s][ind]);

    details::smallInverse(s_[s]);

    linalg::matrixop::gemm(q_[s], s_[s], q_s);

    // Gamma_inv_ -= Q*S^-1*R
    linalg::matrixop::gemm(-1., q_s, r_[s], 1., Gamma_inv_[s]);
  }  // if n
  else {
    Gamma_inv_[s].resizeNoCopy(0);
  }

  Gamma_size_[s] = Gamma_inv_[s].nrRows();
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::pushToEnd() {
  int source;
  int destination;

  for (int s = 0; s < 2; ++s) {
    // Sort in reverse order.
    std::sort(removal_list_[s].rbegin(), removal_list_[s].rend());

    destination = n_max_[s] - 1;

    for (int i = 0; i < removal_list_[s].size(); ++i) {
      source = removal_list_[s][i];

      removal_list_[s][i] = destination;

      linalg::matrixop::swapRowAndCol(M_[s], source, destination);
      configuration_.swapSectorLabels(source, destination, s);

      --destination;
    }
  }

  std::sort(conf_removal_list_.rbegin(), conf_removal_list_.rend());

  destination = configuration_.size() - 1;

  for (int i = 0; i < conf_removal_list_.size(); ++i) {
    source = conf_removal_list_[i];

    conf_removal_list_[i] = destination;

    configuration_.swapVertices(source, destination);

    --destination;
  }
}

// This method is unused and left to potentially use as a testing reference.
template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::recomputeGammaInv() {
  for (int s = 0; s < 2; ++s) {
    if (Gamma_size_[s] > 0)
      linalg::matrixop::inverse(Gamma_inv_[s]);

    Gamma_inv_[s].resize(Gamma_size_[s] + nbr_of_indices_[s]);

    for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
      for (int i = 0; i < Gamma_size_[s]; ++i) {
        Gamma_inv_[s](Gamma_size_[s] + ind, i) = G_[s](sector_indices_[s][ind], move_indices_[s][i]);
        Gamma_inv_[s](i, Gamma_size_[s] + ind) = G_[s](move_indices_[s][i], sector_indices_[s][ind]);
      }

      for (int ind_2 = 0; ind_2 < nbr_of_indices_[s]; ++ind_2) {
        Gamma_inv_[s](Gamma_size_[s] + ind, Gamma_size_[s] + ind_2) =
            G_[s](sector_indices_[s][ind], sector_indices_[s][ind_2]);
      }
      Gamma_inv_[s](Gamma_size_[s] + ind, Gamma_size_[s] + ind) -=
          (1 + gamma_[s].end()[-nbr_of_indices_[s] + ind]) /
          gamma_[s].end()[-nbr_of_indices_[s] + ind];
    }

    Gamma_size_[s] = Gamma_inv_[s].nrRows();

    if (Gamma_size_[s] > 0)
      details::smallInverse(Gamma_inv_[s]);
  }
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP
