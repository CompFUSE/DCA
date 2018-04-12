// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Jérémie Bouquet (bouquetj@gmail.com).

// This class organizes the MC walker in the CT-INT QMC purely on the CPU.

#ifdef DCA_HAVE_CUDA
#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_base.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <ctime>

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/walker_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/move.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/util/integer_division.hpp"


namespace dca {
namespace phys {
namespace solver {
namespace ctint {

template <class Parameters>
class CtintWalkerSubmatrix<linalg::GPU, Parameters>
    : public CtintWalkerBase<linalg::CPU, Parameters> {
public:
  using this_type = CtintWalkerSubmatrix<linalg::GPU, Parameters>;
  using BaseClass = CtintWalkerBase<linalg::CPU, Parameters>;
  using Rng = typename BaseClass::Rng;

  CtintWalkerSubmatrix(Parameters& pars_ref, Rng& rng_ref, const InteractionVertices& vertices,
                       const DMatrixBuilder<linalg::CPU>& builder_ref, int id = 0);

public:
  void doSweep();

  using BaseClass::order;
  using BaseClass::avgOrder;

protected:
  bool tryVertexInsert(bool forced = false);
  bool tryVertexRemoval();

  void doStep();
  void doStep(const int nbr_of_moves_to_delay);
  void generateDelayedMoves(int nbr_of_moves_to_delay);
  void doSubmatrixUpdate();
  void computeAcceptanceProbability();
  void updateGammaInv();
  void updateM();
  void removeRowAndColOfGammaInv(const int s);
  void pushToEnd();
  void computeMInit();
  void computeG0Init();
  void computeGInit();
  void setMFromConfig();
  Move generateMoveType();

  using BaseClass::det_ratio_;

private:
  template <linalg::DeviceType dev>
  using MatrixView = linalg::MatrixView<double, dev>;
  using Matrix = typename BaseClass::Matrix;

  using BaseClass::parameters_;
  using BaseClass::configuration_;
  using BaseClass::rng_;
  using BaseClass::d_builder_;
  using BaseClass::total_interaction_;
  using BaseClass::sign_;
  using BaseClass::M_;

  double max_tau_ = BaseClass::beta_;

  using BaseClass::thermalized_;
  using BaseClass::total_steps_;
  using BaseClass::order_sum_;
  using BaseClass::total_sweeps_;
  using BaseClass::partial_order_sum_;
  using BaseClass::partial_num_steps_;

private:
  int max_submatrix_size_;

  struct DelayedMoveType {
    Move move_type_;

    double removal_rng_;

    int index_;
    int partner_index_;

    bool double_move_ = false;
  };

private:
  //   MatrixPair<linalg::CPU> S_, Q_, R_;
  //   // work spaces
  //   using BaseClass::ipiv_;
  //   using BaseClass::v_work_;
  //   using BaseClass::M_Q_;

  std::vector<DelayedMoveType> delayed_moves_;

  MatrixPair<linalg::CPU> G_;
  MatrixPair<linalg::CPU> G0_;
  MatrixPair<linalg::CPU> D_;
  MatrixPair<linalg::GPU> M_dev_;
  MatrixPair<linalg::GPU> D_dev_;
  MatrixPair<linalg::CPU> Gamma_inv_;
  MatrixPair<linalg::CPU> s_;
  MatrixPair<linalg::CPU> w_;
  std::array<double, 2> d_;
  std::array<double, 2> beta_;
  std::array<std::vector<double>, 2> gamma_;

  std::map<int, double> f_;
  std::map<int, double> prob_const_;
  std::map<std::pair<int, int>, double> gamma_values_;

  int nbr_of_steps_;
  int nbr_of_moves_to_delay_;

  int nbr_of_moves_;
  int Gamma_size_;
  int Gamma_index_;
  int index_;
  Move move_type_;

  using BaseClass::concurrency_;
  using BaseClass::thread_id_;
  std::array<cudaStream_t, 2> stream_;

  // Initial and current sector sizes.

  int n_init_;
  int n_;

  // Maximal sector size after submatrix update.

  int n_max_;

  std::vector<int> move_indices_;
  std::vector<int> existing_indices_;
  std::vector<int> removal_list_;
  std::vector<int> insertion_list_;
  std::vector<int> insertion_Gamma_indices_;

  std::vector<int>::iterator insertion_list_it_;
  bool recently_added_;
  bool accepted_;

  double acceptance_probability_;
  bool do_nothing_;
  bool double_;

  // For testing.

  bool force_acceptance_ = false;
};

template <class Parameters>
CtintWalkerSubmatrix<linalg::GPU, Parameters>::CtintWalkerSubmatrix(
    Parameters& parameters_ref, Rng& rng_ref, const InteractionVertices& vertices,
    const DMatrixBuilder<linalg::CPU>& builder_ref, int id)
    : BaseClass(parameters_ref, rng_ref, vertices, builder_ref, id) {
  for (int s = 0; s < 2; ++s)
    stream_[s] = linalg::util::getStream(thread_id_, s);

  while (parameters_.getInitialConfigurationSize() > configuration_.size())
    configuration_.insertRandom(rng_);
  setMFromConfig();

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

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\nCT-INT GPU submatrix walker created." << std::endl;
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::doSweep() {
  nbr_of_moves_to_delay_ = parameters_.getMaxSubmatrixSize();

  if (not thermalized_)
    nbr_of_steps_ = avgOrder() / nbr_of_moves_to_delay_ + 1;
  else
    nbr_of_steps_ = BaseClass::nb_steps_per_sweep_ * parameters_.get_sweeps_per_measurement() /
                    nbr_of_moves_to_delay_;

  double f_i;

  for (int s = 0; s < 2; ++s) {
    for (int i = 0; i < M_[s].size().first; ++i) {
      f_i = f_[configuration_.getSector(s).getAuxFieldType(i)];

      for (int j = 0; j < M_[s].size().second; ++j) {
        M_[s](i, j) /= -(f_i - 1);
      }
    }
  }

  for (int i = 0; i < nbr_of_steps_; i++) {
    doStep();
    ++total_steps_;
    order_sum_ += order();
  }
  ++total_sweeps_;

  for (int s = 0; s < 2; ++s) {
    for (int i = 0; i < M_[s].size().first; ++i) {
      f_i = f_[configuration_.getSector(s).getAuxFieldType(i)];

      for (int j = 0; j < M_[s].size().second; ++j) {
        M_[s](i, j) *= -(f_i - 1);
      }
    }
  }

  // Keep tha average after half thermalization for deciding the order.
  if ((not thermalized_) and (total_sweeps_ == parameters_.get_warm_up_sweeps() / 2)) {
    partial_order_sum_ = order_sum_;
    partial_num_steps_ = total_steps_;
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::doStep() {
  generateDelayedMoves(nbr_of_moves_to_delay_);
  doSubmatrixUpdate();
}

// Do one step with arbitrary number of moves. For testing.

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::doStep(const int nbr_of_moves_to_delay) {
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
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::generateDelayedMoves(int nbr_of_moves_to_delay) {
  assert(nbr_of_moves_to_delay > 0);

  const int max_nbr_of_moves = parameters_.getMaxSubmatrixSize();

  nbr_of_moves_to_delay = std::min(nbr_of_moves_to_delay, max_nbr_of_moves);

  delayed_moves_.clear();

  nbr_of_moves_ = 0;

  n_init_ = configuration_.size();
  n_ = n_init_;

  // assert(std::make_pair(n_init_, n_init_) == M_[0].size() && M_[0].size() == M_[1].size());

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

      // Depending on the parameters, insertRandom() may sometimes insert two vertices at once.

      if (delta > 1) {
        std::cout << "\nWARNING Double insertion.\n" << std::endl;
        delayed_move.double_move_ = true;
        delayed_move.partner_index_ = configuration_.size() - 2;
      }
    }

    else {
      std::cout << "\nWARNING Unkown move type encountered and ignored." << std::endl;
      continue;
    }

    delayed_moves_.push_back(delayed_move);

    ++nbr_of_moves_;
  }

  n_max_ = configuration_.size();

  existing_indices_.clear();

  for (int i = 0; i < n_init_; ++i) {
    existing_indices_.push_back(i);
  }

  removal_list_.clear();

  for (int i = n_init_; i < n_max_; ++i) {
    removal_list_.push_back(i);
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::doSubmatrixUpdate() {
  computeMInit();
  computeGInit();

  Gamma_size_ = 0;
  for (int s = 0; s < 2; ++s) {
    Gamma_inv_[s].resizeNoCopy(0);
    gamma_[s].clear();
  }

  move_indices_.clear();
  insertion_list_.clear();
  insertion_Gamma_indices_.clear();

  int aux_spin_type, new_aux_spin_type;

  Matrix result_matrix_1, result_matrix_2;

  for (int delay_ind = 0; delay_ind < delayed_moves_.size(); ++delay_ind) {
    for (int s = 0; s < 2; ++s) {
      auto& gamma = gamma_[s];
      auto& sector = configuration_.getSector(s);
      move_type_ = delayed_moves_[delay_ind].move_type_;
      do_nothing_ = false;
      recently_added_ = false;
      double_ = false;

      if (move_type_ == INSERTION) {
        index_ = delayed_moves_[delay_ind].index_;
        new_aux_spin_type = sector.getAuxFieldType(index_);

        aux_spin_type = 0;
      }
      else if (move_type_ == REMOVAL) {
        // Do nothing if there is no vertex to remove.

        if (existing_indices_.size() == 0) {
          do_nothing_ = true;
          break;
        }
        else {
          index_ =
              existing_indices_[int(delayed_moves_[delay_ind].removal_rng_ * existing_indices_.size())];
          new_aux_spin_type = 0;
          aux_spin_type = sector.getAuxFieldType(index_);

          // Check if the vertex to remove was inserted during the current submatrix update.

          if (index_ >= n_init_) {
            insertion_list_it_ = std::find(insertion_list_.begin(), insertion_list_.end(), index_);

            recently_added_ = true;
            Gamma_index_ =
                insertion_Gamma_indices_[std::distance(insertion_list_.begin(), insertion_list_it_)];
          }
        }
      }

      // if (!recently_added_)
      // 	gamma.push_back(d_builder_.computeGamma(aux_spin_type, new_aux_spin_type));
      // else {
      // 	gamma.push_back(d_builder_.computeGamma(new_aux_spin_type, aux_spin_type));

      if (!recently_added_)
        gamma.push_back(gamma_values_[std::make_pair(aux_spin_type, new_aux_spin_type)]);
      else {
        gamma.push_back(gamma_values_[std::make_pair(new_aux_spin_type, aux_spin_type)]);

        if (s == 0) {
          move_indices_.erase(std::remove(move_indices_.begin(), move_indices_.end(), index_),
                              move_indices_.end());

          --Gamma_size_;
        }
      }

      if (Gamma_size_ > 0) {
        s_[s].resize(std::make_pair(Gamma_size_, 1));
        w_[s].resize(std::make_pair(1, Gamma_size_));

        for (int i = 0; i < Gamma_size_; ++i) {
          s_[s](i, 0) = G_[s](move_indices_[i], index_);
          w_[s](0, i) = G_[s](index_, move_indices_[i]);
        }

        d_[s] = G_[s](index_, index_) - (1 + gamma.back()) / gamma.back();

        if (recently_added_)
          removeRowAndColOfGammaInv(s);

        result_matrix_1.resizeNoCopy(s_[s].size());
        linalg::matrixop::gemm(Gamma_inv_[s], s_[s], result_matrix_1);
        result_matrix_2.resizeNoCopy(std::make_pair(1, 1));
        linalg::matrixop::gemm(w_[s], result_matrix_1, result_matrix_2);

        beta_[s] = d_[s] - result_matrix_2(0, 0);
      }
      else {
        Gamma_inv_[s].resize(0);

        beta_[s] = G_[s](index_, index_) - (1 + gamma.back()) / gamma.back();
      }
    }

    // Directly go to next move if do_nothing_ is true.

    if (do_nothing_)
      continue;

    if (recently_added_) {
      insertion_list_.erase(std::remove(insertion_list_.begin(), insertion_list_.end(), index_),
                            insertion_list_.end());

      insertion_Gamma_indices_.erase(std::remove(insertion_Gamma_indices_.begin(),
                                                 insertion_Gamma_indices_.end(), Gamma_index_),
                                     insertion_Gamma_indices_.end());

      for (int s = 0; s < 2; ++s) {
        gamma_[s].erase(gamma_[s].begin() + Gamma_index_);
      }

      for (int i = 0; i < insertion_Gamma_indices_.size(); ++i) {
        if (insertion_Gamma_indices_[i] > Gamma_index_)
          --insertion_Gamma_indices_[i];
      }
    }

    // Compute acceptance probability.

    computeAcceptanceProbability();

    accepted_ = rng_() < std::abs(acceptance_probability_);

    if (force_acceptance_)
      accepted_ = true;

    // Update GammaInv if necessary.

    updateGammaInv();

    // Update other objects.

    if (accepted_) {
      if (acceptance_probability_ < 0)
        sign_ *= -1;

      if (!recently_added_)
        move_indices_.push_back(index_);
      else {
        for (int s = 0; s < 2; ++s)
          gamma_[s].pop_back();
      }

      if (move_type_ == INSERTION) {
        ++n_;
        insertion_list_.push_back(index_);
        insertion_Gamma_indices_.push_back(Gamma_size_ - 1);

        existing_indices_.push_back(index_);

        removal_list_.erase(std::remove(removal_list_.begin(), removal_list_.end(), index_),
                            removal_list_.end());
      }
      else if (move_type_ == REMOVAL) {
        --n_;
        removal_list_.push_back(index_);
        existing_indices_.erase(
            std::remove(existing_indices_.begin(), existing_indices_.end(), index_),
            existing_indices_.end());
      }
    }

    // If the move is rejected:

    else {
      if (!recently_added_) {
        for (int s = 0; s < 2; ++s) {
          gamma_[s].pop_back();
        }
      }
      else {
        insertion_list_.push_back(index_);
        insertion_Gamma_indices_.push_back(Gamma_size_ - 1);

        move_indices_.push_back(index_);
        existing_indices_.push_back(index_);
      }
    }
  }

  updateM();
}

template <class Parameters>
Move CtintWalkerSubmatrix<linalg::GPU, Parameters>::generateMoveType() {
  if (rng_() <= 0.5)
    return INSERTION;
  else
    return REMOVAL;
}

// Extend M by adding non-interacting vertices.
// M is not computed again and should be up-to-date.
template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::computeMInit() {
  const int delta = n_max_ - n_init_;

  if (delta > 0) {
    double f_j;

    for (int s = 0; s < 2; ++s) {
      auto& D = D_[s];
      D.resizeNoCopy(std::make_pair(delta, n_init_));
      d_builder_.computeG0(D, configuration_.getSector(s), n_init_, n_max_, 0);

      for (int j = 0; j < n_init_; ++j) {
        f_j = f_[configuration_.getSector(s).getAuxFieldType(j)] - 1;

        for (int i = 0; i < delta; ++i) {
          D(i, j) *= f_j;
        }
      }
      D_dev_[s].setAsync(D, stream_[s]);

      M_[s].resize(n_max_);
      M_dev_[s].setAsync(M_[s], stream_[s]);

      MatrixView<linalg::GPU> M(M_dev_[s], 0, 0, n_init_, n_init_);
      MatrixView<linalg::GPU> D_M(M_dev_[s], n_init_, 0, delta, n_init_);

      linalg::matrixop::gemm(D_dev_[s], M, D_M, thread_id_, s);

      M_[s].setAsync(M_dev_[s], stream_[s]);
    }

    for (int s = 0; s < 2; ++s) {
      cudaStreamSynchronize(stream_[s]);

      for (int i = 0; i < n_max_; ++i) {
        for (int j = n_init_; j < n_max_; ++j) {
          M_[s](i, j) = 0;
        }
      }

      for (int i = n_init_; i < n_max_; ++i)
        M_[s](i, i) = 1;
    }
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::computeGInit() {
  const int delta = n_max_ - n_init_;

  double f;
  Matrix G0(std::make_pair(n_max_, delta));

  for (int s = 0; s < 2; ++s) {
    G_[s].resizeNoCopy(n_max_);

    for (int j = 0; j < n_init_; ++j) {
      f = f_[configuration_.getSector(s).getAuxFieldType(j)];

      for (int i = 0; i < n_max_; ++i) {
        G_[s](i, j) = (M_[s](i, j) * f - double(i == j)) / (f - 1);
      }
    }

    if (delta > 0) {
      d_builder_.computeG0(G0, configuration_.getSector(s), n_init_, n_max_, 1);

      MatrixView<linalg::CPU> G(G_[s], 0, n_init_, n_max_, delta);

      linalg::matrixop::gemm(M_[s], G0, G);
    }
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::computeG0Init() {
  for (int s = 0; s < 2; ++s) {
    d_builder_.computeG0Init(G0_[s], configuration_.getSector(s), n_init_, n_max_);
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::computeAcceptanceProbability() {
  double K = total_interaction_ * prob_const_[configuration_.getSector(0).getAuxFieldType(index_)];

  acceptance_probability_ = 1;

  for (int s = 0; s < 2; ++s) {
    acceptance_probability_ *= beta_[s] * gamma_[s].back();
  }

  if (recently_added_)
    acceptance_probability_ = 1 / acceptance_probability_;

  if (move_type_ == INSERTION)
    acceptance_probability_ *= K / (n_ + 1);
  else if (move_type_ == REMOVAL)
    acceptance_probability_ *= n_ / K;
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::updateGammaInv() {
  if ((!recently_added_ && accepted_) || (recently_added_ && !accepted_)) {
    if (Gamma_size_ > 0) {
      for (int s = 0; s < 2; ++s) {
        Gamma_inv_[s].resize(Gamma_size_ + 1);

        MatrixView<linalg::CPU> bulk(Gamma_inv_[s], 0, 0, Gamma_size_, Gamma_size_);
        MatrixView<linalg::CPU> s_inv(Gamma_inv_[s], 0, Gamma_size_, Gamma_size_, 1);
        MatrixView<linalg::CPU> w_inv(Gamma_inv_[s], Gamma_size_, 0, 1, Gamma_size_);

        linalg::matrixop::gemm(-1 / beta_[s], bulk, s_[s], 0., s_inv);
        linalg::matrixop::gemm(-1 / beta_[s], w_[s], bulk, 0., w_inv);
        linalg::matrixop::gemm(beta_[s], s_inv, w_inv, 1., bulk);

        Gamma_inv_[s](Gamma_size_, Gamma_size_) = 1 / beta_[s];
      }
    }
    else {
      for (int s = 0; s < 2; ++s) {
        Gamma_inv_[s].resizeNoCopy(1);
        Gamma_inv_[s](0, 0) = 1 / beta_[s];
      }
    }

    ++Gamma_size_;
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::updateM() {
  // assert(Gamma_size_ == gamma_[0].size() && gamma_[0].size() == gamma_[1].size());

  if (Gamma_size_ > 0) {
    int p;

    Matrix old_G(std::make_pair(n_max_, Gamma_size_));
    Matrix old_M(std::make_pair(Gamma_size_, n_max_));
    Matrix result_matrix(std::make_pair(Gamma_size_, n_max_));

    for (int s = 0; s < 2; ++s) {
      for (int j = 0; j < Gamma_size_; ++j) {
        p = move_indices_[j];

        for (int i = 0; i < n_max_; ++i) {
          old_G(i, j) = G_[s](i, p);
          old_M(j, i) = M_[s](p, i);
        }
      }

      linalg::matrixop::gemm(Gamma_inv_[s], old_M, result_matrix);
      linalg::matrixop::gemm(-1., old_G, result_matrix, 1., M_[s]);

      p = 0;
      for (auto& i : move_indices_) {
        for (int j = 0; j < n_max_; ++j) {
          M_[s](i, j) /= 1 + gamma_[s][p];
        }
        ++p;
      }
    }
  }

  // Remove "non-interacting" rows and columns.

  pushToEnd();

  for (int i = 0; i < removal_list_.size(); ++i)
    configuration_.pop();

  for (int s = 0; s < 2; ++s) {
    M_[s].resize(n_max_ - removal_list_.size());

    // G0_[s].resize(n_max_ - removal_list_.size());
  }
}

// Remove row and column of Gamma_inv with Woodbury's formula.

// Gamma <- Gamma - U.V => GammaInv <- GammaInv + GammaInv.U.(Id-V.GammaInv.U)^(-1).V.GammaInv.

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::removeRowAndColOfGammaInv(const int s) {
  Matrix U(std::make_pair(Gamma_size_ + 1, 2));
  Matrix V(std::make_pair(2, Gamma_size_ + 1));

  Matrix matrix_1(std::make_pair(Gamma_size_ + 1, 2));
  Matrix matrix_2(std::make_pair(2, Gamma_size_ + 1));
  Matrix matrix_3(std::make_pair(2, 2));

  for (int i = 0; i < Gamma_index_; ++i) {
    U(i, 0) = 0;
    U(i, 1) = s_[s](i, 0);
    V(0, i) = w_[s](0, i);
    V(1, i) = 0;
  }

  for (int i = Gamma_index_ + 1; i < Gamma_size_ + 1; ++i) {
    U(i, 0) = 0;
    U(i, 1) = s_[s](i - 1, 0);
    V(0, i) = w_[s](0, i - 1);
    V(1, i) = 0;
  }

  U(Gamma_index_, 0) = 1;
  U(Gamma_index_, 1) = -1;
  V(0, Gamma_index_) = d_[s];
  V(1, Gamma_index_) = 1;

  linalg::matrixop::gemm(Gamma_inv_[s], U, matrix_1);
  linalg::matrixop::gemm(V, Gamma_inv_[s], matrix_2);
  linalg::matrixop::gemm(V, matrix_1, matrix_3);

  matrix_3(0, 0) = 1 - matrix_3(0, 0);
  matrix_3(1, 1) = 1 - matrix_3(1, 1);
  matrix_3(0, 1) *= -1;
  matrix_3(1, 0) *= -1;

  linalg::matrixop::inverse(matrix_3);

  linalg::matrixop::gemm(matrix_1, matrix_3, U);
  linalg::matrixop::gemm(1., U, matrix_2, 1., Gamma_inv_[s]);

  linalg::matrixop::removeRowAndCol(Gamma_inv_[s], Gamma_index_);
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::pushToEnd() {
  // Sort in reverse order.

  std::sort(removal_list_.rbegin(), removal_list_.rend());

  for (int s = 0; s < 2; ++s) {
    int source;
    int destination = n_max_ - 1;

    for (int i = 0; i < removal_list_.size(); ++i) {
      source = removal_list_[i];

      if (s == 1)
        removal_list_[i] = destination;

      linalg::matrixop::swapRowAndCol(M_[s], source, destination);
      // linalg::matrixop::swapRowAndCol(G0_[s], source, destination);

      configuration_.swapSectorLabels(source, destination, s);

      configuration_.swapVertices(source, destination);

      --destination;
    }
  }
}

template <class Parameters>
void CtintWalkerSubmatrix<linalg::GPU, Parameters>::setMFromConfig() {
  sign_ = 1;
  for (int s = 0; s < 2; ++s) {
    const auto& sector = configuration_.getSector(s);
    auto& M = M_[s];
    const int n = sector.size();
    M.resize(n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        M(i, j) = d_builder_.computeD(i, j, sector);
      }
    }
    linalg::matrixop::inverse(M);

    if (linalg::matrixop::determinant(M) < 0)
      sign_ *= -1;
  }
}

// Useless.

template <class Parameters>
bool CtintWalkerSubmatrix<linalg::GPU, Parameters>::tryVertexInsert(bool forced) {
  return forced;
}

template <class Parameters>
bool CtintWalkerSubmatrix<linalg::GPU, Parameters>::tryVertexRemoval() {
  return 0;
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_GPU_SUBMATRIX_HPP
#endif  // DCA_HAVE_CUDA
