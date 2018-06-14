// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Jérémie Bouquet (bouquetj@gmail.com).
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.com).

// This class organizes the MC walker in the CT-INT QMC purely on the CPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_base.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <ctime>

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/walker_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/move.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
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

  CtintWalkerSubmatrix(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
                       const DMatrixBuilder<linalg::CPU>& builder_ref, int id = 0);

public:
  void doSweep();

  using BaseClass::order;

protected:
  void doSteps();
  void generateDelayedMoves(int nbr_of_moves_to_delay);
  void mainSubmatrixProcess();
  void updateM();

  // For testing purposes.
  void doStep(const int nbr_of_moves_to_delay);

private:
  virtual void doStep();
  void doSubmatrixUpdate();
  double computeAcceptanceProbability() const;
  void updateGammaInv();
  void removeRowAndColOfGammaInv(const int s);
  void pushToEnd();
  void computeMInit();
  void computeG0Init();
  void computeGInit();
  void setMFromConfig();
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
  using BaseClass::det_ratio_;
  using BaseClass::n_bands_;

  double max_tau_ = BaseClass::beta_;

  using BaseClass::thermalized_;

protected:
  int max_submatrix_size_;

  struct DelayedMoveType {
    Move move_type_;

    double removal_rng_;

    int index_;
    int partner_index_;
    std::array<std::vector<int>, 2> sector_indices_;
    std::array<std::vector<int>, 2> partner_sector_indices_;

    bool double_move_ = false;
  };

protected:
  std::vector<DelayedMoveType> delayed_moves_;

  template <linalg::DeviceType matrix_device>
  using MatrixPair = std::array<linalg::Matrix<double, matrix_device>, 2>;
  MatrixPair<linalg::CPU> G_;
  MatrixPair<linalg::CPU> G0_;
  MatrixPair<linalg::CPU> Gamma_inv_;
  MatrixPair<linalg::CPU> s_;
  MatrixPair<linalg::CPU> w_;
  MatrixPair<linalg::CPU> s_2_;
  MatrixPair<linalg::CPU> w_2_;
  std::array<double, 2> d_;
  std::array<double, 2> d_1_2_;
  std::array<double, 2> d_2_1_;
  std::array<double, 2> d_2_2_;
  std::array<double, 2> beta_;
  std::array<double, 2> beta_2_;
  std::array<std::vector<double>, 2> gamma_;

  std::map<int, double> f_;
  std::map<int, double> prob_const_;
  std::map<std::pair<int, int>, double> gamma_values_;

  using BaseClass::nb_steps_per_sweep_;
  int nbr_of_steps_;
  int nbr_of_submatrix_steps_;
  int nbr_of_moves_to_delay_;
  int max_nbr_of_moves_;

  std::array<int, 2> Gamma_size_;
  std::array<std::vector<int>, 2> Gamma_indices_;
  std::array<std::vector<int>, 2> sector_indices_;
  Move move_type_;

  std::array<int, 2> nbr_of_indices_;
  std::array<bool, 2> found_sector_index_;

  int index_;
  double tau_;

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

  bool double_;

  Matrix result_matrix_1_, result_matrix_2_, result_matrix_3_, result_matrix_4_, result_matrix_5_;
  Matrix result_matrix_6_, result_matrix_7_, result_matrix_8_, result_matrix_9_, result_matrix_10_;
  Matrix D_;

  // For testing.

  bool force_acceptance_ = false;
};

template <class Parameters>
CtintWalkerSubmatrix<linalg::CPU, Parameters>::CtintWalkerSubmatrix(
    Parameters& parameters_ref, Rng& rng_ref, const InteractionVertices& vertices,
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

  if (not thermalized_)
    nbr_of_steps_ = BaseClass::avgOrder() + 1;
  else
    nbr_of_steps_ = nb_steps_per_sweep_ * parameters_.get_sweeps_per_measurement();

  // Get the maximum of Monte Carlo steps/moves that can be performed during one submatrix step.
  max_nbr_of_moves_ = parameters_.getMaxSubmatrixSize();

  BaseClass::n_steps_ += nbr_of_steps_;
  while (nbr_of_steps_ > 0) {
    nbr_of_moves_to_delay_ = std::min(nbr_of_steps_, max_nbr_of_moves_);
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
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::doStep(const int nbr_of_moves_to_delay) {
  std::cout << "\nStarted doStep() function for testing." << std::endl;

  force_acceptance_ = true;

  double f_i;

  for (int s = 0; s < 2; ++s) {
    for (int i = 0; i < M_[s].size().first; ++i) {
      f_i = f_[configuration_.getSector(s).getAuxFieldType(i)];

      for (int j = 0; j < M_[s].size().second; ++j) {
        M_[s](i, j) /= -(f_i - 1);
      }
    }
  }

  generateDelayedMoves(nbr_of_moves_to_delay);

  std::cout << "\nGenerated " << nbr_of_moves_to_delay << " moves for testing.\n" << std::endl;

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
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::generateDelayedMoves(int nbr_of_moves_to_delay) {
  assert(nbr_of_moves_to_delay > 0);

  delayed_moves_.clear();

  int nbr_of_moves_ = 0;

  for (int s = 0; s < 2; ++s) {
    n_init_[s] = configuration_.getSector(s).size();
    n_[s] = n_init_[s];
  }

  // Generate delayed moves.

  while (nbr_of_moves_ < nbr_of_moves_to_delay) {
    DelayedMoveType delayed_move;

    delayed_move.move_type_ = generateMoveType();

    if (delayed_move.move_type_ == REMOVAL) {
      delayed_move.removal_rng_ = rng_();
    }

    else if (delayed_move.move_type_ == INSERTION) {
      configuration_.insertRandom(rng_);
      const int delta = configuration_.lastInsertionSize();

      delayed_move.index_ = configuration_.size() - 1;

      for (int s = 0; s < 2; ++s) {
        delayed_move.sector_indices_[s].push_back(configuration_.getSector(s).size() - 1);
      }

      // Depending on the parameters, insertRandom() may sometimes insert two vertices at once.

      if (delta > 1) {
        std::cout << "\nWARNING Double insertion.\n" << std::endl;
        // ....
      }
    }

    else {
      std::cout << "\nWARNING Unkown move type encountered and ignored." << std::endl;
      continue;
    }

    delayed_moves_.push_back(delayed_move);

    ++nbr_of_moves_;
  }

  existing_indices_.clear();

  for (int i = 0; i < (n_init_[0] + n_init_[1]) / 2; ++i) {
    existing_indices_.push_back(i);
  }

  for (int s = 0; s < 2; ++s) {
    n_max_[s] = configuration_.getSector(s).size();

    removal_list_[s].clear();

    for (int i = n_init_[s]; i < n_max_[s]; ++i) {
      removal_list_[s].push_back(i);
    }
  }

  conf_removal_list_.clear();

  for (int i = (n_init_[0] + n_init_[1]) / 2; i < (n_max_[0] + n_max_[1]) / 2; ++i) {
    conf_removal_list_.push_back(i);
  }
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
    for (int s = 0; s < 2; ++s) {
      assert(Gamma_size_[s] == Gamma_inv_[s].size().first);

      move_type_ = delayed_moves_[delay_ind].move_type_;
      do_nothing = false;
      recently_added_ = false;
      double_ = false;

      if (move_type_ == INSERTION) {
        index_ = delayed_moves_[delay_ind].index_;

        // Index in configuration may differ from index in each sector,
        // therefore look for corresponding indices.

        findSectorIndices(s);

        // Continue to next sector/move if there is nothing to change in current sector.

        if (!found_sector_index_[s])
          continue;

        new_aux_spin_type.clear();
        aux_spin_type.clear();
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          new_aux_spin_type.push_back(
              configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));
          aux_spin_type.push_back(0);
        }
      }  // if(insertion)

      else if (move_type_ == REMOVAL) {
        // Do nothing if there is no vertex to remove.

        if (existing_indices_.size() == 0) {
          do_nothing = true;
          break;
        }
        else {
          index_ =
              existing_indices_[int(delayed_moves_[delay_ind].removal_rng_ * existing_indices_.size())];

          // Index in configuration may differ from index in each sector,
          // therefore look for corresponding indices.

          findSectorIndices(s);

          // Continue to next sector/move if there is nothing to change in current sector.

          if (!found_sector_index_[s])
            continue;

          new_aux_spin_type.clear();
          aux_spin_type.clear();
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            new_aux_spin_type.push_back(0);
            aux_spin_type.push_back(
                configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));
          }

          // Check if the vertex to remove was inserted during the current submatrix update.

          if (sector_indices_[s][0] >= n_init_[s]) {
            recently_added_ = true;

            Gamma_indices_[s].clear();
            for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
              insertion_list_it_ = std::find(insertion_list_[s].begin(), insertion_list_[s].end(),
                                             sector_indices_[s][ind]);
              Gamma_indices_[s].push_back(insertion_Gamma_indices_[s][std::distance(
                  insertion_list_[s].begin(), insertion_list_it_)]);
            }
          }
        }
      }  // else if(removal).

      //(What follows is executed for both insertions or removals).

      for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
        if (!recently_added_)
          gamma_[s].push_back(
              gamma_values_[std::make_pair(aux_spin_type[ind], new_aux_spin_type[ind])]);
        else {
          gamma_[s].push_back(
              gamma_values_[std::make_pair(new_aux_spin_type[ind], aux_spin_type[ind])]);

          move_indices_[s].erase(
              std::remove(move_indices_[s].begin(), move_indices_[s].end(), sector_indices_[s][ind]),
              move_indices_[s].end());

          --Gamma_size_[s];
        }
      }

      if (Gamma_size_[s] > 0) {
        if (nbr_of_indices_[s] == 1 || nbr_of_indices_[s] == 2) {
          s_[s].resize(std::make_pair(Gamma_size_[s], 1));
          w_[s].resize(std::make_pair(1, Gamma_size_[s]));
  // TODO: this looks fishy.
          for (int i = 0; i < Gamma_size_[s]; ++i) {
            s_[s](i, 0) = G_[s](move_indices_[s][i], sector_indices_[s][0]);
            w_[s](0, i) = G_[s](sector_indices_[s][0], move_indices_[s][i]);
          }

          d_[s] = G_[s](sector_indices_[s][0], sector_indices_[s][0]) -
                  (1 + gamma_[s].end()[-nbr_of_indices_[s]]) / gamma_[s].end()[-nbr_of_indices_[s]];

          if (recently_added_)
            removeRowAndColOfGammaInv(s);

          result_matrix_1_.resizeNoCopy(s_[s].size());
          linalg::matrixop::gemm(Gamma_inv_[s], s_[s], result_matrix_1_);
          result_matrix_2_.resizeNoCopy(std::make_pair(1, 1));
          linalg::matrixop::gemm(w_[s], result_matrix_1_, result_matrix_2_);

          beta_[s] = d_[s] - result_matrix_2_(0, 0);
        }
        else {
          beta_[s] = 1;
        }

        // In case of a two-indices update a given sector, obtain "beta" from direct formula
        // to avoid having to compute an intermediate state of GammaInv.

        // TODO: add code for update of more than 2 indices per sector, if this is supposed to occur
        // (?).

        if (nbr_of_indices_[s] == 2) {
          s_2_[s].resize(std::make_pair(Gamma_size_[s], 1));
          w_2_[s].resize(std::make_pair(1, Gamma_size_[s]));

          for (int i = 0; i < Gamma_size_[s]; ++i) {
            s_2_[s](i, 0) = G_[s](move_indices_[s][i], sector_indices_[s][1]);
            w_2_[s](0, i) = G_[s](sector_indices_[s][1], move_indices_[s][i]);
          }

          d_2_2_[s] = G_[s](sector_indices_[s][1], sector_indices_[s][1]) -
                      (1 + gamma_[s].back()) / gamma_[s].back();
          d_1_2_[s] = G_[s](sector_indices_[s][0], sector_indices_[s][1]);
          d_2_1_[s] = G_[s](sector_indices_[s][1], sector_indices_[s][1]);

          linalg::matrixop::gemm(Gamma_inv_[s], s_2_[s], result_matrix_1_);
          result_matrix_3_.resizeNoCopy(std::make_pair(1, 1));
          linalg::matrixop::gemm(w_[s], result_matrix_1_, result_matrix_3_);

          result_matrix_5_.resizeNoCopy(std::make_pair(1, 1));
          linalg::matrixop::gemm(w_2_[s], result_matrix_1_, result_matrix_5_);

          linalg::matrixop::gemm(Gamma_inv_[s], s_[s], result_matrix_1_);
          result_matrix_4_.resizeNoCopy(std::make_pair(1, 1));
          linalg::matrixop::gemm(w_2_[s], result_matrix_1_, result_matrix_4_);

          beta_2_[s] = (d_2_2_[s] - result_matrix_5_(0, 0));
          beta_2_[s] -=
              (d_1_2_[s] - result_matrix_3_(0, 0)) * (d_2_1_[s] - result_matrix_4_(0, 0)) / beta_[s];
        }
        else
          beta_2_[s] = 1;
      }
      else {
        if (nbr_of_indices_[s] == 1 || nbr_of_indices_[s] == 2) {
          Gamma_inv_[s].resize(0);

          d_[s] = beta_[s] =
              G_[s](sector_indices_[s][0], sector_indices_[s][0]) -
              (1 + gamma_[s].end()[-nbr_of_indices_[s]]) / gamma_[s].end()[-nbr_of_indices_[s]];
        }
        else {
          beta_[s] = 1;
        }

        // Case of a two-indices update.

        if (nbr_of_indices_[s] == 2) {
          s_2_[s].resize(std::make_pair(1, 1));
          w_2_[s].resize(std::make_pair(1, 1));

          s_2_[s](0, 0) = G_[s](sector_indices_[s][0], sector_indices_[s][1]);
          w_2_[s](0, 0) = G_[s](sector_indices_[s][1], sector_indices_[s][0]);

          d_2_2_[s] = G_[s](sector_indices_[s][1], sector_indices_[s][1]) -
                      (1 + gamma_[s].back()) / gamma_[s].back();

          beta_2_[s] = d_2_2_[s] - w_2_[s](0, 0) * s_2_[s](0, 0) / beta_[s];
        }
        else
          beta_2_[s] = 1;
      }

      if (recently_added_) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          insertion_list_[s].erase(std::remove(insertion_list_[s].begin(), insertion_list_[s].end(),
                                               sector_indices_[s][ind]),
                                   insertion_list_[s].end());

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
    }  // s loop.

    // Directly go to next move if do_nothing is true.

    if (do_nothing)
      continue;

    // Compute acceptance probability.

    const double acceptance_probability = computeAcceptanceProbability();

    accepted_ = rng_() < std::abs(acceptance_probability);

    // Acceptance can be forced (for testing).

    if (force_acceptance_)
      accepted_ = true;

    // Update GammaInv if necessary.

    updateGammaInv();

    // NB: recomputeGammaInv is just a inefficient alternative to updateGammaInv. Only for testing
    // or debbuging.
    //
    // recomputeGammaInv();

    // Update other objects.

    if (accepted_) {
      if (acceptance_probability < 0)
        sign_ *= -1;

      if (!recently_added_)
        for (int s = 0; s < 2; ++s)
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
            move_indices_[s].push_back(sector_indices_[s][ind]);
      else {
        for (int s = 0; s < 2; ++s)
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
            gamma_[s].pop_back();
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
      if (!recently_added_) {
        for (int s = 0; s < 2; ++s) {
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
            gamma_[s].pop_back();
        }
      }
      else {
        for (int s = 0; s < 2; ++s) {
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            insertion_list_[s].push_back(sector_indices_[s][ind]);
            insertion_Gamma_indices_[s].push_back(Gamma_size_[s] - nbr_of_indices_[s] + ind);

            move_indices_[s].push_back(sector_indices_[s][ind]);
          }
        }
        existing_indices_.push_back(index_);
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

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::computeG0Init() {
  for (int s = 0; s < 2; ++s) {
    d_builder_.computeG0Init(G0_[s], configuration_.getSector(s), n_init_[s], n_max_[s]);
  }
}

template <class Parameters>
double CtintWalkerSubmatrix<linalg::CPU, Parameters>::computeAcceptanceProbability() const {
  double acceptance_probability = 1;

  for (int s = 0; s < 2; ++s) {
    acceptance_probability *= beta_[s] * beta_2_[s];
    for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
      acceptance_probability *= gamma_[s].end()[-nbr_of_indices_[s] + ind];
  }

  double K;

  if (recently_added_)
    acceptance_probability = 1 / acceptance_probability;

  if (sector_indices_[0].size() != 0)
    K = total_interaction_ *
        prob_const_.at(configuration_.getSector(0).getAuxFieldType(sector_indices_[0][0]));
  else
    K = total_interaction_ *
        prob_const_.at(configuration_.getSector(1).getAuxFieldType(sector_indices_[1][0]));

  int n = (n_[0] + n_[1]) / 2;

  if (move_type_ == INSERTION)
    acceptance_probability *= K / (n + 1);
  else if (move_type_ == REMOVAL)
    acceptance_probability *= n / K;

  return acceptance_probability;
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::updateGammaInv() {
  if ((!recently_added_ && accepted_) || (recently_added_ && !accepted_)) {
    for (int s = 0; s < 2; ++s) {
      if (Gamma_size_[s] > 0) {
        if (nbr_of_indices_[s] == 1 || nbr_of_indices_[s] == 2) {
          Gamma_inv_[s].resize(Gamma_size_[s] + 1);

          MatrixView bulk(Gamma_inv_[s], 0, 0, Gamma_size_[s], Gamma_size_[s]);
          MatrixView s_inv(Gamma_inv_[s], 0, Gamma_size_[s], Gamma_size_[s], 1);
          MatrixView w_inv(Gamma_inv_[s], Gamma_size_[s], 0, 1, Gamma_size_[s]);

          linalg::matrixop::gemm(-1 / beta_[s], bulk, s_[s], 0., s_inv);
          linalg::matrixop::gemm(-1 / beta_[s], w_[s], bulk, 0., w_inv);
          linalg::matrixop::gemm(beta_[s], s_inv, w_inv, 1., bulk);

          Gamma_inv_[s](Gamma_size_[s], Gamma_size_[s]) = 1 / beta_[s];

          if (nbr_of_indices_[s] == 2) {
            Gamma_inv_[s].resize(Gamma_size_[s] + 2);

            MatrixView bulk(Gamma_inv_[s], 0, 0, Gamma_size_[s] + 1, Gamma_size_[s] + 1);
            MatrixView s_inv(Gamma_inv_[s], 0, Gamma_size_[s] + 1, Gamma_size_[s] + 1, 1);
            MatrixView w_inv(Gamma_inv_[s], Gamma_size_[s] + 1, 0, 1, Gamma_size_[s] + 1);

            s_2_[s].resize(std::make_pair(Gamma_size_[s] + 1, 1));
            s_2_[s](Gamma_size_[s], 0) = d_1_2_[s];
            w_2_[s].resize(std::make_pair(1, Gamma_size_[s] + 1));
            w_2_[s](0, Gamma_size_[s]) = d_2_1_[s];

            linalg::matrixop::gemm(-1 / beta_2_[s], bulk, s_2_[s], 0., s_inv);
            linalg::matrixop::gemm(-1 / beta_2_[s], w_2_[s], bulk, 0., w_inv);
            linalg::matrixop::gemm(beta_2_[s], s_inv, w_inv, 1., bulk);

            Gamma_inv_[s](Gamma_size_[s] + 1, Gamma_size_[s] + 1) = 1 / beta_2_[s];
          }
        }
        else if (nbr_of_indices_[s] > 0) {
          // TODO: Add code of update of more than 2 indices if needed...
          std::cout << "\n\t\t\tERROR nbr_of_indices_[" << s << "] > 2 (" << nbr_of_indices_[s]
                    << ")." << std::endl;
        }
      }
      else {
        if (nbr_of_indices_[s] == 1 || nbr_of_indices_[s] == 2) {
          Gamma_inv_[s].resizeNoCopy(nbr_of_indices_[s]);
          Gamma_inv_[s](0, 0) = 1 / beta_[s];

          if (nbr_of_indices_[s] == 2) {
            Gamma_inv_[s](0, 0) += s_2_[s](0, 0) * w_2_[s](0, 0) / pow(beta_[s], 2) / beta_2_[s];
            Gamma_inv_[s](0, 1) = -s_2_[s](0, 0) / beta_[s] / beta_2_[s];
            Gamma_inv_[s](1, 0) = -w_2_[s](0, 0) / beta_[s] / beta_2_[s];
            Gamma_inv_[s](1, 1) = 1 / beta_2_[s];
          }
        }
        else if (nbr_of_indices_[s] > 0) {
          // TODO: Add code of update of more than 2 indices if needed...
          std::cout << "\n\t\t\tERROR nbr_of_indices_[" << s << "] > 2 (" << nbr_of_indices_[s]
                    << ")." << std::endl;
        }
      }

      Gamma_size_[s] += nbr_of_indices_[s];
    }
  }
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

  //  if (nbr_of_indices_[s] > 2) {
  //    std::cout << "\n\n\t\t\t**********ERROR**********\n\n" << std::endl;
  //    std::cout << "Nbr of indices = " << nbr_of_indices_[s] << std::endl;
  //    std::cout << "s              = " << s << std::endl;
  //    for (int i = 0; i < nbr_of_indices_[s]; ++i)
  //      std::cout << sector_indices_[s][i] << std::endl;
  //  }
}

// Remove row and column of Gamma_inv with Woodbury's formula.

// Gamma <- Gamma - U.V => GammaInv <- GammaInv + GammaInv.U.(Id-V.GammaInv.U)^(-1).V.GammaInv.

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::removeRowAndColOfGammaInv(const int s) {
  Matrix U(std::make_pair(Gamma_size_[s] + nbr_of_indices_[s], 2 * nbr_of_indices_[s]));
  Matrix V(std::make_pair(2 * nbr_of_indices_[s], Gamma_size_[s] + nbr_of_indices_[s]));

  result_matrix_1_.resize(std::make_pair(Gamma_size_[s] + nbr_of_indices_[s], 2 * nbr_of_indices_[s]));
  result_matrix_2_.resize(std::make_pair(2 * nbr_of_indices_[s], Gamma_size_[s] + nbr_of_indices_[s]));
  result_matrix_3_.resize(std::make_pair(2 * nbr_of_indices_[s], 2 * nbr_of_indices_[s]));

  if (nbr_of_indices_[s] == 1) {
    for (int i = 0; i < Gamma_indices_[s][0]; ++i) {
      U(i, 0) = 0;
      U(i, 1) = s_[s](i, 0);
      V(0, i) = w_[s](0, i);
      V(1, i) = 0;
    }

    for (int i = Gamma_indices_[s][0] + 1; i < Gamma_size_[s] + 1; ++i) {
      U(i, 0) = 0;
      U(i, 1) = s_[s](i - 1, 0);
      V(0, i) = w_[s](0, i - 1);
      V(1, i) = 0;
    }

    U(Gamma_indices_[s][0], 0) = 1;
    U(Gamma_indices_[s][0], 1) = -1;
    V(0, Gamma_indices_[s][0]) = d_[s];
    V(1, Gamma_indices_[s][0]) = 1;

    linalg::matrixop::gemm(Gamma_inv_[s], U, result_matrix_1_);
    linalg::matrixop::gemm(V, Gamma_inv_[s], result_matrix_2_);
    linalg::matrixop::gemm(V, result_matrix_1_, result_matrix_3_);

    result_matrix_3_(0, 0) = 1 - result_matrix_3_(0, 0);
    result_matrix_3_(1, 1) = 1 - result_matrix_3_(1, 1);
    result_matrix_3_(0, 1) *= -1;
    result_matrix_3_(1, 0) *= -1;

    linalg::matrixop::inverse(result_matrix_3_);

    linalg::matrixop::gemm(result_matrix_1_, result_matrix_3_, U);
    linalg::matrixop::gemm(1., U, result_matrix_2_, 1., Gamma_inv_[s]);

    linalg::matrixop::removeRowAndCol(Gamma_inv_[s], Gamma_indices_[s][0]);
  }

  // Code for any nbr_of_indices_[s] > 0.

  else {
    int i_max;
    for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
      for (int i = 0; i < Gamma_indices_[s][ind]; ++i) {
        U(i, ind) = 0;
        U(i, nbr_of_indices_[s] + ind) = G_[s](move_indices_[s][i], sector_indices_[s][ind]);
        V(ind, i) = G_[s](sector_indices_[s][ind], move_indices_[s][i]);
        V(nbr_of_indices_[s] + ind, i) = 0;
      }

      if (ind < nbr_of_indices_[s] - 1)
        i_max = Gamma_indices_[s][ind + 1];
      else
        i_max = Gamma_size_[s] + nbr_of_indices_[s];

      for (int i = Gamma_indices_[s][ind] + 1; i < i_max; ++i) {
        U(i, ind) = 0;
        U(i, nbr_of_indices_[s] + ind) =
            G_[s](move_indices_[s][i - 1 - ind], sector_indices_[s][ind]);
        V(ind, i) = G_[s](sector_indices_[s][ind], move_indices_[s][i - 1 - ind]);
        V(nbr_of_indices_[s] + ind, i) = 0;
      }

      U(Gamma_indices_[s][ind], ind) = 1;
      U(Gamma_indices_[s][ind], nbr_of_indices_[s] + ind) = -1;
      V(ind, Gamma_indices_[s][ind]) = G_[s](sector_indices_[s][ind], sector_indices_[s][ind]) -
                                       (1 + gamma_[s].end()[-nbr_of_indices_[s] + ind]) /
                                           gamma_[s].end()[-nbr_of_indices_[s] + ind];
      V(nbr_of_indices_[s] + ind, Gamma_indices_[s][ind]) = 1;
    }

    linalg::matrixop::gemm(Gamma_inv_[s], U, result_matrix_1_);
    linalg::matrixop::gemm(V, Gamma_inv_[s], result_matrix_2_);
    linalg::matrixop::gemm(V, result_matrix_1_, result_matrix_3_);

    for (int i = 0; i < result_matrix_3_.size().first; ++i) {
      for (int j = 0; j < result_matrix_3_.size().second; ++j) {
        result_matrix_3_(i, j) *= -1;
      }
      result_matrix_3_(i, i) += 1;
    }

    linalg::matrixop::inverse(result_matrix_3_);

    linalg::matrixop::gemm(result_matrix_1_, result_matrix_3_, U);
    linalg::matrixop::gemm(1., U, result_matrix_2_, 1., Gamma_inv_[s]);

    for (int ind = nbr_of_indices_[s] - 1; ind >= 0; --ind) {
      linalg::matrixop::removeRowAndCol(Gamma_inv_[s], Gamma_indices_[s][ind]);
    }
  }
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

template <class Parameters>
void CtintWalkerSubmatrix<linalg::CPU, Parameters>::recomputeGammaInv() {
  if ((!recently_added_ && accepted_) || (recently_added_ && !accepted_)) {
    for (int s = 0; s < 2; ++s) {
      if (Gamma_size_[s] > 0)
        linalg::matrixop::inverse(Gamma_inv_[s]);

      Gamma_inv_[s].resize(Gamma_size_[s] + nbr_of_indices_[s]);

      for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
        for (int i = 0; i < Gamma_size_[s]; ++i) {
          Gamma_inv_[s](Gamma_size_[s] + ind, i) =
              G_[s](sector_indices_[s][ind], move_indices_[s][i]);
          Gamma_inv_[s](i, Gamma_size_[s] + ind) =
              G_[s](move_indices_[s][i], sector_indices_[s][ind]);
        }

        for (int ind_2 = 0; ind_2 < nbr_of_indices_[s]; ++ind_2) {
          Gamma_inv_[s](Gamma_size_[s] + ind, Gamma_size_[s] + ind_2) =
              G_[s](sector_indices_[s][ind], sector_indices_[s][ind_2]);
        }
        Gamma_inv_[s](Gamma_size_[s] + ind, Gamma_size_[s] + ind) -=
            (1 + gamma_[s].end()[-nbr_of_indices_[s] + ind]) /
            gamma_[s].end()[-nbr_of_indices_[s] + ind];
      }

      Gamma_size_[s] += nbr_of_indices_[s];

      if (Gamma_size_[s] > 0)
        linalg::matrixop::inverse(Gamma_inv_[s]);
    }
  }
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP
