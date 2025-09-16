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
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <ctime>

#include "dca/linalg/linalg.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/move.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/walker_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_submatrix_base.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {

template <class Parameters, DistType DIST>
class CtintWalkerSubmatrixCpu : public CtintWalkerSubmatrixBase<Parameters, DIST> {
public:
  using this_type = CtintWalkerSubmatrixCpu;
  using SubmatrixBase = CtintWalkerSubmatrixBase<Parameters, DIST>;
  using BaseClass = CtintWalkerBase<Parameters, DIST>;
  using typename BaseClass::Rng;
  using typename BaseClass::Data;
  using typename BaseClass::Profiler;
  using typename BaseClass::GpuStream;
  using typename BaseClass::Real;
  using typename BaseClass::Scalar;

  using Resource = DMatrixBuilder<linalg::CPU, Scalar>;

  CtintWalkerSubmatrixCpu(const Parameters& pars_ref, const Data& /*data*/, Rng& rng_ref,
                          DMatrixBuilder<linalg::CPU, Scalar>& d_matrix_builder, int id = 0);

  virtual ~CtintWalkerSubmatrixCpu() = default;

  void computeM(typename BaseClass::MatrixPair& m_accum);

  using BaseClass::order;

  void setMFromConfig() override;

  void markThermalized() override;

protected:
  void doSteps();
  void generateDelayedMoves(int nbr_of_movesto_delay);

  void updateM() override;

  DMatrixBuilder<linalg::CPU, Scalar>& d_matrix_builder_;

  BaseClass::MatrixPair getM();

  /** The following methods are really only here to get decent unit testing
      they shouldn't really be called outside of the base implementations
  */
  void computeMInit() override;
  void computeGInit() override;
  BaseClass::MatrixPair getRawM();
  BaseClass::MatrixPair getRawG();

private:
  /** This does a bunch of things, this is the majority of a step
   *  For each delayed_move it:
   *  * computes the acceptance prob
   *  * compares std::abs(acceptance prob) and random value
   *    * Does an update of the log_weight for a configuration
   *        mc_weight_ratio = acceptance_probability;
   *      For a new vertex
   *        mc_weight ratio *= weight_term
   *      For a removed vertes
   *        mc_weight_ratio /= weight_term
   *      with
   *        weight_term = prob_const_[field_type][b] * the vertex interaction strength;
   *        BaseClass::mc_log_weight_ += std::log(std::abs(mc_weight_ratio));
   */
  void doSubmatrixUpdate();

  /** returns [acceptance_probability , mc_weight_ratio ]
   */
  auto computeAcceptanceProbability();

  void updateGammaInv(int s);

  void removeRowAndColOfGammaInv();

  //  void computeG0Init();

  Move generateMoveType();

  void computeInsertionMatrices(const std::vector<int>& indices, const int s);

  void computeRemovalMatrix(int s);

  void computeMixedInsertionAndRemoval(int s);

  void findSectorIndices(const int s);

  void recomputeGammaInv();

protected:
  using MatrixView = linalg::MatrixView<Scalar, linalg::CPU>;
  using Matrix = linalg::Matrix<Scalar, linalg::CPU>;

  using BaseClass::parameters_;
  using BaseClass::configuration_;
  using BaseClass::rng_;
  using BaseClass::thread_id_;
  using BaseClass::total_interaction_;
  using BaseClass::phase_;
  using BaseClass::M_;
  using BaseClass::n_bands_;
  using BaseClass::beta_;

  using BaseClass::thermalized_;
  using BaseClass::sweeps_per_meas_;
  using BaseClass::partial_order_avg_;
  using BaseClass::thermalization_steps_;
  using BaseClass::order_avg_;
  using BaseClass::sign_avg_;
  using BaseClass::n_steps_;
  using BaseClass::mc_log_weight_;
  using BaseClass::n_accepted_;

protected:
  using SubmatrixBase::max_submatrix_size_;

  using typename SubmatrixBase::DelayedMoveType;
  using SubmatrixBase::f_;
  using SubmatrixBase::gamma_values_;
  using SubmatrixBase::prob_const_;
  using SubmatrixBase::nb_steps_per_sweep_;
  using SubmatrixBase::n_max_;
  using SubmatrixBase::n_init_;
  using SubmatrixBase::D_;
  using SubmatrixBase::G_;
  using SubmatrixBase::M_;
  using SubmatrixBase::s_;
  using SubmatrixBase::Gamma_indices_;
  using SubmatrixBase::Gamma_inv_;
  using SubmatrixBase::Gamma_inv_cpy_;
  using SubmatrixBase::Gamma_q_;
  using SubmatrixBase::workspace_;
  using SubmatrixBase::r_;
  using SubmatrixBase::move_indices_;
  using SubmatrixBase::gamma_;
  using SubmatrixBase::source_list_;
  using SubmatrixBase::removal_list_;
  using SubmatrixBase::conf_removal_list_;
  using SubmatrixBase::sector_indices_;
  using SubmatrixBase::index_;
  using SubmatrixBase::nbr_of_indices_;
  using SubmatrixBase::q_;
  using SubmatrixBase::det_ratio_;

protected:
  using BaseClass::acceptance_prob_;

protected:
  using BaseClass::flop_;
};

template <class Parameters, DistType DIST>
CtintWalkerSubmatrixCpu<Parameters, DIST>::CtintWalkerSubmatrixCpu(
    const Parameters& parameters_ref, const Data& data, Rng& rng_ref,
    DMatrixBuilder<linalg::CPU, Scalar>& d_matrix_builder, int id)
    : SubmatrixBase(parameters_ref, data, rng_ref, id), d_matrix_builder_(d_matrix_builder) {
  for (int b = 0; b < n_bands_; ++b) {
    for (int i = 1; i <= 3; ++i) {
      f_[i][b] = d_matrix_builder_.computeF(i, b);
      f_[-i][b] = d_matrix_builder_.computeF(-i, b);

      gamma_values_[std::make_pair(0, i)][b] = d_matrix_builder_.computeGamma(0, i, b);
      gamma_values_[std::make_pair(0, -i)][b] = d_matrix_builder_.computeGamma(0, -i, b);
      gamma_values_[std::make_pair(i, 0)][b] = d_matrix_builder_.computeGamma(i, 0, b);
      gamma_values_[std::make_pair(-i, 0)][b] = d_matrix_builder_.computeGamma(-i, 0, b);

      prob_const_[i][b] = prob_const_[-i][b] = -1. / (f_[i][b] - 1) / (f_[-i][b] - 1);
    }
    f_[0][b] = 1;
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::setMFromConfig() {
  BaseClass::setMFromConfigImpl(d_matrix_builder_);
  SubmatrixBase::transformM();
#ifdef DEBUG_SUBMATRIX
  std::cout << "cpu M post setMFromConfigImpl: \n";
  M_[0].set_name("cpu_M_0");
  M_[1].set_name("cpu_M_1");
  M_[0].print();
  M_[1].print();
#endif
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::markThermalized() {
  thermalized_ = true;

  nb_steps_per_sweep_ = std::max(1., std::ceil(sweeps_per_meas_ * partial_order_avg_.mean()));
  thermalization_steps_ = n_steps_;

  order_avg_.reset();
  sign_avg_.reset();
  n_accepted_ = 0;

  // Recompute the Monte Carlo weight.
  setMFromConfig();
#ifndef NDEBUG
  // writeAlphas();
#endif
}

// Extend M by adding non-interacting vertices.
// M is not computed again and should be up-to-date.
template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::computeMInit() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    const int delta = n_max_[s] - n_init_[s];

    if (delta > 0) {
      Scalar f_j;
      D_.resize(std::make_pair(delta, n_init_[s]));

      if (delta == 0 || n_init_[s] == 0)
        throw std::runtime_error(
            "expansion factor dropped to 0 or below use a higher beta or larger interaction!");

      d_matrix_builder_.computeG0(D_, configuration_.getSector(s), n_init_[s], n_max_[s], 0);

      std::array<linalg::Vector<Real, linalg::CPU>, 2> f_values;
      f_values[s].resize(n_init_[s]);
      for (int j = 0; j < n_init_[s]; ++j) {
        const auto field_type = configuration_.getSector(s).getAuxFieldType(j);
        const auto b = configuration_.getSector(s).getRightB(j);
        f_j = f_[field_type][b] - 1;
        f_values[s][j] = f_[field_type][b];
        for (int i = 0; i < delta; ++i) {
          D_(i, j) *= f_j;
        }
      }

#ifdef DEBUG_SUBMATRIX
      f_values[0].set_name("cpu_f_values_0");
      f_values[1].set_name("cpu_f_values_1");
      f_values[0].print();
      f_values[1].print();
      std::cout << "cpu D post f factor mult\n";
      D_.print();
      using namespace dca::addt_str_oper;
      std::cout << "M_[" << s << "] size: " << M_[s].size() << '\n';
#endif
      M_[s].resize(n_max_[s]);

      MatrixView M(M_[s], 0, 0, n_init_[s], n_init_[s]);
      MatrixView D_M(M_[s], n_init_[s], 0, delta, n_init_[s]);

#ifdef DEBUG_SUBMATRIX
      std::cout << "cpu M pre gemm\n";
      M_[s].print();
#endif

      linalg::matrixop::gemm(D_, M, D_M);
      flop_ += 2 * D_.nrRows() * D_.nrCols() * M.nrCols();

#ifdef DEBUG_SUBMATRIX
      D_M.print();
#endif

      for (int i = 0; i < n_max_[s]; ++i) {
        for (int j = n_init_[s]; j < n_max_[s]; ++j) {
          M_[s](i, j) = 0;
        }
      }

      for (int i = n_init_[s]; i < n_max_[s]; ++i)
        M_[s](i, i) = 1;
#ifdef DEBUG_SUBMATRIX
      M_[s].print();
#endif
    }
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::computeGInit() {
  Profiler profiler(__FUNCTION__, "CT-INT walker", __LINE__, thread_id_);

  for (int s = 0; s < 2; ++s) {
    // Doing subtraction with unsigned types...
    const int delta = n_max_[s] - n_init_[s];

    // f?  thanks
    Scalar f;
    Matrix G0(std::make_pair(n_max_[s], delta));

    G_[s].resizeNoCopy(n_max_[s]);

    for (int j = 0; j < n_init_[s]; ++j) {
      const auto field_type = configuration_.getSector(s).getAuxFieldType(j);
      const auto b = configuration_.getSector(s).getRightB(j);
      f = f_[field_type][b];

      for (int i = 0; i < n_max_[s]; ++i) {
        G_[s](i, j) = (M_[s](i, j) * f - Real(i == j)) / (f - 1.0);
      }
    }

    if (delta > 0) {
      d_matrix_builder_.computeG0(G0, configuration_.getSector(s), n_init_[s], n_max_[s], 1);

      MatrixView G(G_[s], 0, n_init_[s], n_max_[s], delta);

      linalg::matrixop::gemm(M_[s], G0, G);
      flop_ += 2. * M_[s].nrRows() * M_[s].nrCols() * G0.nrCols();
    }
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::updateGammaInv(int s) {
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
    linalg::matrixop::gemm(Scalar(-1.), Gamma_q, s_inv, Scalar(0.), q_inv);

    auto& r_Gamma = workspace_;
    r_Gamma.resizeNoCopy(r_[s].size());
    linalg::matrixop::gemm(r_[s], bulk, r_Gamma);
    linalg::matrixop::gemm(Scalar(-1.), s_inv, r_Gamma, Scalar(0.), r_inv);

    // Gamma_ += Gamma_ * q_ * s_^-1 * r_ * Gamma_
    linalg::matrixop::gemm(Scalar(-1.), q_inv, r_Gamma, Scalar(1.), bulk);
  }
  else {
    Gamma_inv_[s].resizeNoCopy(delta);
    details::smallInverse(s_[s], Gamma_inv_[s]);
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::updateM() {
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
      linalg::matrixop::gemm(Scalar(-1.), old_G, result_matrix, Scalar(1.), M_[s]);
      flop_ += 2 * Gamma_inv_[s].nrRows() * Gamma_inv_[s].nrCols() * old_M.nrCols();
      flop_ += 2 * old_G.nrRows() * old_G.nrCols() * result_matrix.nrCols();

      p = 0;
      for (auto& i : move_indices_[s]) {
        for (int j = 0; j < n_max_[s]; ++j)
          M_[s](i, j) /= 1.0 + gamma_[s][p];
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

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::findSectorIndices(const int s) {
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
template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::removeRowAndColOfGammaInv() {
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
      linalg::matrixop::gemm(Scalar(-1.), q_s, r_[s], Scalar(1.), Gamma_inv_[s]);
    }  // if n
    else {
      Gamma_inv_[s].resizeNoCopy(0);
    }
  }
}

// This method is unused and left to potentially use as a testing reference.
template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::recomputeGammaInv() {
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

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::computeM(typename BaseClass::MatrixPair& m_accum) {
  for (int s = 0; s < 2; ++s) {
    m_accum[s].resizeNoCopy(M_[s].size());

    for (int j = 0; j < M_[s].size().second; ++j) {
      for (int i = 0; i < M_[s].size().first; ++i) {
        const auto field_type = configuration_.getSector(s).getAuxFieldType(i);
        const auto b = configuration_.getSector(s).getLeftB(i);
        const Scalar factor = -(f_[field_type][b] - 1.);
        m_accum[s](i, j) = M_[s](i, j) * factor;
      }
    }
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::computeInsertionMatrices(
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
    linalg::matrixop::gemm(Scalar(-1.), r_[s], Gamma_q, Scalar(1.), s_[s]);
  }
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::computeRemovalMatrix(const int s) {
  const int delta = Gamma_indices_[s].size();
  s_[s].resizeNoCopy(delta);
  for (int j = 0; j < delta; ++j)
    for (int i = 0; i < delta; ++i)
      s_[s](i, j) = Gamma_inv_[s](Gamma_indices_[s][i], Gamma_indices_[s][j]);
}

template <class Parameters, DistType DIST>
void CtintWalkerSubmatrixCpu<Parameters, DIST>::computeMixedInsertionAndRemoval(int s) {
  Gamma_inv_cpy_[s] = Gamma_inv_[s];

  if (sector_indices_[s].size() == 0) {
    s_[s].resizeNoCopy(0);
    return;
  }

  // TODO: avoid reallocation.
  std::vector<int> insertion_indices;
  for (int i = 0; i < nbr_of_indices_[s]; ++i)
    if (!SubmatrixBase::recentlyAdded(i, s))
      insertion_indices.push_back(sector_indices_[s][i]);
  assert(Gamma_indices_[s].size() + insertion_indices.size() == sector_indices_.size());

  computeInsertionMatrices(insertion_indices, s);
  det_ratio_ *= details::smallDeterminant(s_[s]);
  updateGammaInv(s);

  computeRemovalMatrix(s);
}

template <class Parameters, DistType DIST>
CtintWalkerSubmatrixCpu<Parameters, DIST>::BaseClass::MatrixPair CtintWalkerSubmatrixCpu<
    Parameters, DIST>::getRawM() {
  typename BaseClass::MatrixPair M;
  M = M_;
  M[0].set_name("subMatrixCPU::M[0]");
  M[1].set_name("subMatrixCPU::M[1]");
  return M;
}

template <class Parameters, DistType DIST>
CtintWalkerSubmatrixCpu<Parameters, DIST>::BaseClass::MatrixPair CtintWalkerSubmatrixCpu<
    Parameters, DIST>::getRawG() {
  typename BaseClass::MatrixPair G;
  G = G_;
  G[0].set_name("subMatrixCPU::G[0]");
  G[1].set_name("subMatrixCPU::G[1]");
  return G;
}

template <class Parameters, DistType DIST>
CtintWalkerSubmatrixCpu<Parameters, DIST>::BaseClass::MatrixPair CtintWalkerSubmatrixCpu<Parameters,
                                                                                         DIST>::getM() {
  typename BaseClass::MatrixPair M;
  computeM(M);
  return M;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CPU_SUBMATRIX_HPP
