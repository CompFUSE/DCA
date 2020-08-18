// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// Cluster Monte Carlo integrator based on a CT-INT algorithm.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_CTINT_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_CTINT_CLUSTER_SOLVER_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

#include "dca/config/mc_options.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/statistics/util.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/accumulator/ctint_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/domains/common_domains.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_choice.hpp"
//#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/time_correlator.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <linalg::DeviceType device_t, class Parameters, bool use_submatrix = false>
class CtintClusterSolver {
public:
  using Real = typename config::McOptions::MCScalar;

  using Data = DcaData<Parameters>;
  static constexpr linalg::DeviceType device = device_t;

  CtintClusterSolver(Parameters& parameters_ref, Data& Data_ref);

  ~CtintClusterSolver();

  // Initialize g0_interpolation and reset internal state. Must be called before integrate.
  void initialize(int dca_iteration = 0);

  // Monte Carlo walk and accumulation.
  void integrate();

  // gather the walker's measurements and average them across the processes.
  // Then it computes the final integration results.
  // Postcondition: DcaData contains Sigma, G_r_t, G_r_w, G_k_w, G_k_t
  double finalize();

  // Calls finalize(). In addition:
  // Postcondition: dca_info_struct contains metadata on the integration.
  // Returns: L2 difference between Sigma and Sigma_cluster.
  double finalize(DcaLoopData<Parameters>& dca_info_struct);

  double computeDensity() const;

  template <class Writer>
  void write(Writer& /*writer*/) {}

  // For testing purposes.
  // Returns the function G(k,w) without averaging across MPI ranks.
  auto local_G_k_w() const;

protected:  // thread jacket interface.
  using ParametersType = Parameters;
  using DataType = Data;
  using Rng = typename Parameters::random_number_generator;
  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename Parameters::concurrency_type;
  using Lattice = typename Parameters::lattice_type;

  using Walker = ctint::CtintWalkerChoice<device_t, Parameters, use_submatrix, Real>;
  using Accumulator = ctint::CtintAccumulator<Parameters, device_t, Real>;

private:
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using KDmn = typename CDA::KClusterDmn;
  using RDmn = typename CDA::RClusterDmn;
  using Nu = func::dmn_variadic<BDmn, SDmn>;
  using TDmn = func::dmn_0<domains::time_domain>;
  using WDmn = func::dmn_0<domains::frequency_domain>;
  using LabelDomain = func::dmn_variadic<BDmn, BDmn, RDmn>;

  using SpGreensFunction =
      func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, KDmn, WDmn>>;

protected:  // Protected for testing purposes.
  void warmUp();

  void measure();

  void computeG_k_w(const SpGreensFunction& G0, const SpGreensFunction& M_k_w,
                    SpGreensFunction& G_k_w) const;

  void computeSigma(const SpGreensFunction& G, const SpGreensFunction& G0,
                    SpGreensFunction& Sigma) const;

  // Returns: average sign.
  double gatherMAndG4(SpGreensFunction& M, bool compute_error) const;

  double L2Difference() const;

  void computeErrorBars() const {}

protected:
  Parameters& parameters_;
  Concurrency& concurrency_;
  Data& data_;

  // Stores the integration result
  Accumulator accumulator_;

  double total_time_ = 0;
  double warm_up_time_ = 0;
  int dca_iteration_ = 0;

private:
  bool perform_tp_accumulation_;
  const LabelDomain label_dmn_;
  std::unique_ptr<Walker> walker_;
  // Walker input.
  ctint::G0Interpolation<device_t, Real> g0_;
  Rng rng_;
};

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
CtintClusterSolver<device_t, Parameters, use_submatrix>::CtintClusterSolver(Parameters& parameters_ref,
                                                                            Data& data_ref)
    : parameters_(parameters_ref),
      concurrency_(parameters_.get_concurrency()),
      data_(data_ref),

      accumulator_(parameters_, data_),

      rng_(concurrency_.id(), concurrency_.number_of_processors(), parameters_.get_seed()) {
  Walker::setDMatrixBuilder(g0_);
  //  TimeCorrelator<Parameters, Real, device_t>::setG0(g0_);
  Walker::setInteractionVertices(data_, parameters_);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-INT Integrator is born \n\n";
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
CtintClusterSolver<device_t, Parameters, use_submatrix>::~CtintClusterSolver() {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-INT Integrator has died \n\n";
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
void CtintClusterSolver<device_t, Parameters, use_submatrix>::initialize(int dca_iteration) {
  dca_iteration_ = dca_iteration;

  g0_.initialize(ctint::details::shrinkG0(data_.G0_r_t_cluster_excluded));

  Walker::setDMatrixAlpha(parameters_.getAlphas(), parameters_.adjustAlphaDd());

  perform_tp_accumulation_ =
      parameters_.isAccumulatingG4() && dca_iteration == parameters_.get_dca_iterations() - 1;
  accumulator_.initialize(dca_iteration_);
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
void CtintClusterSolver<device_t, Parameters, use_submatrix>::integrate() {
  walker_ = std::make_unique<Walker>(parameters_, data_, rng_, 0);
  walker_->initialize(dca_iteration_);

  dca::profiling::WallTime start_time;
  auto getTime = [&]() {
    dca::profiling::Duration split_time(dca::profiling::WallTime(), start_time);
    return split_time.sec + 1.e-6 * split_time.usec;
  };

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t Warm up has started.\n" << std::endl;
  warmUp();
  warm_up_time_ = getTime();

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t Measuring has started.\n\n";
  measure();

  total_time_ = getTime();

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\tMeasuring has ended. Done " << parameters_.get_measurements()[dca_iteration_]
              << " measurements.\n";
    walker_->printSummary();
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
double CtintClusterSolver<device_t, Parameters, use_submatrix>::finalize() {
  bool compute_error = false;
  if (dca_iteration_ == parameters_.get_dca_iterations() - 1) {
    if (parameters_.get_error_computation_type() == ErrorComputationType::JACK_KNIFE) {
      if (concurrency_.id() == concurrency_.first())
        std::cout << "Computing jack knife error.\n\n";
      compute_error = true;
    }
    else if (parameters_.get_error_computation_type() == ErrorComputationType::STANDARD_DEVIATION)
      std::cout << "CT-INT does not support ErrorComputationType::STANDARD_DEVIATION.\n"
                << "Error computation will be skipped.\n";
  }

  SpGreensFunction M;

  // average M across ranks.
  double avg_sign = gatherMAndG4(M, compute_error);

  // compute G_r_t and save it into data_.
  computeG_k_w(data_.G0_k_w_cluster_excluded, M, data_.G_k_w);
  symmetrize::execute<Lattice>(data_.G_k_w);

  // transform  G_k_w and save into data_.
  math::transform::FunctionTransform<KDmn, RDmn>::execute(data_.G_k_w, data_.G_r_w);
  symmetrize::execute<Lattice>(data_.G_r_w);

  // compute and  save Sigma into data_
  // TODO: check if a better estimate exists
  computeSigma(data_.G_k_w, data_.G0_k_w_cluster_excluded, data_.Sigma);

  // compute error
  if (compute_error) {
    data_.get_Sigma_error() = concurrency_.jackknifeError(data_.Sigma);
    data_.get_G_k_w_error() = concurrency_.jackknifeError(data_.G_k_w);
    if (perform_tp_accumulation_) {
      for (int channel = 0; channel < data_.get_G4().size(); ++channel)
        data_.get_G4_error()[channel] = concurrency_.jackknifeError(data_.get_G4()[channel]);
    }
  }

  // Fourier transform the Green's function.
  // TODO check and write function
  auto G_k_w_copy = data_.G_k_w;
  G_k_w_copy -= data_.G0_k_w_cluster_excluded;
  math::transform::FunctionTransform<WDmn, TDmn>::execute(G_k_w_copy, data_.G_k_t);
  data_.G_k_t += data_.G0_k_t_cluster_excluded;
  math::transform::FunctionTransform<KDmn, RDmn>::execute(data_.G_k_t, data_.G_r_t);

  auto local_time = total_time_;
  concurrency_.sum(total_time_);
  auto gflop = accumulator_.getFLOPs() * 1e-9;
  concurrency_.sum(gflop);

  if (concurrency_.id() == 0) {
    std::cout << "\n\t\t Collected measurements \t" << dca::util::print_time() << "\n"
              << "\n\t\t\t QMC-local-time : " << local_time << " [sec]"
              << "\n\t\t\t QMC-total-time : " << total_time_ << " [sec]"
              << "\n\t\t\t Gflop   : " << gflop << " [Gf]"
              << "\n\t\t\t Gflop/s   : " << gflop / local_time << " [Gf/s]"
              << "\n\t\t\t sign     : " << avg_sign << "\n\t\t\t Density = " << computeDensity()
              << "\n"
              << std::endl;
  }

  return avg_sign;
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
double CtintClusterSolver<device_t, Parameters, use_submatrix>::finalize(
    DcaLoopData<Parameters>& loop_data) {
  double avg_sign = finalize();
  // Compute and save into loop_data Sigma_zero_mom and std deviation
  for (int nu = 0; nu < Nu::dmn_size(); nu++) {
    for (int k = 0; k < KDmn::dmn_size(); k++) {
      std::vector<double> x;
      for (int l = 0; l < WDmn::dmn_size() / 4; l++)
        x.push_back(real(data_.Sigma(nu, nu, k, l)));

      loop_data.Sigma_zero_moment(nu, k, dca_iteration_) = math::statistics::util::mean(x);
      loop_data.standard_deviation(nu, k, dca_iteration_) =
          math::statistics::util::standard_deviation(x);
    }
  }

  loop_data.average_expansion_order(dca_iteration_) = accumulator_.avgOrder();
  loop_data.sign(dca_iteration_) = avg_sign;                            // This is already averaged.
  loop_data.MC_integration_per_mpi_task(dca_iteration_) = total_time_;  // This is already averaged.
  loop_data.thermalization_per_mpi_task(dca_iteration_) = warm_up_time_;

  if (dca_iteration_ == parameters_.get_dca_iterations() - 1) {
    concurrency_.delayedSum(loop_data.Sigma_zero_moment);
    concurrency_.delayedSum(loop_data.standard_deviation);
    concurrency_.delayedSum(loop_data.average_expansion_order);
    concurrency_.delayedSum(loop_data.thermalization_per_mpi_task);

    concurrency_.resolveSums();

    loop_data.Sigma_zero_moment /= concurrency_.number_of_processors();
    loop_data.Sigma_zero_moment /= concurrency_.number_of_processors();
    loop_data.standard_deviation /= concurrency_.number_of_processors();
    loop_data.average_expansion_order /= concurrency_.number_of_processors();
    loop_data.thermalization_per_mpi_task /= concurrency_.number_of_processors();
  }

  // Free walker memory for the dca loop.
  walker_.release();
  // Compute and return L2 difference between Sigma and Sigma cluster.
  return loop_data.L2_Sigma_difference(dca_iteration_) = L2Difference();
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
void CtintClusterSolver<device_t, Parameters, use_submatrix>::warmUp() {
  const int n_sweep = parameters_.get_warm_up_sweeps();
  for (int i = 0; i < n_sweep; i++) {
    walker_->doSweep();

    walker_->updateShell(i, n_sweep);
  }

  walker_->markThermalized();
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
void CtintClusterSolver<device_t, Parameters, use_submatrix>::measure() {
  const int n_meas = parallel::util::getWorkload(parameters_.get_measurements()[dca_iteration_], 1,
                                                 0, concurrency_);

  for (int i = 0; i < n_meas; i++) {
    {
      Profiler profiler("updating", "QMCI", __LINE__);
      walker_->doSweep();
    }
    {
      Profiler profiler("measurements", "QMCI", __LINE__);
      accumulator_.accumulate(*walker_);
    }
    walker_->updateShell(i, n_meas);
  }

  accumulator_.finalize();
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
void CtintClusterSolver<device_t, Parameters, use_submatrix>::computeG_k_w(
    const SpGreensFunction& G0, const SpGreensFunction& M_k_w, SpGreensFunction& G_k_w) const {
  const int matrix_dim = Nu::dmn_size();
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_M(matrix_dim, matrix_dim);

  const char op = 'N';
  const double one_over_beta = 1. / parameters_.get_beta();

  // G(w) = g0_(w) - g0_(w) *M(w)* g0_(w) / beta
  G_k_w = G0;
  for (int k_ind = 0; k_ind < KDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < WDmn::dmn_size(); w_ind++) {
      // G0_M <- G0 * M
      dca::linalg::blas::gemm(&op, &op, matrix_dim, matrix_dim, matrix_dim, 1.,
                              &G0(0, 0, k_ind, w_ind), matrix_dim, &M_k_w(0, 0, k_ind, w_ind),
                              matrix_dim, 0., G0_M.ptr(), G0_M.leadingDimension());

      // G -= G0 M G0 / beta
      dca::linalg::blas::gemm(&op, &op, matrix_dim, matrix_dim, matrix_dim, -one_over_beta,
                              G0_M.ptr(), G0_M.leadingDimension(), &G0(0, 0, k_ind, w_ind),
                              matrix_dim, 1., &G_k_w(0, 0, k_ind, w_ind), matrix_dim);
    }
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
void CtintClusterSolver<device_t, Parameters, use_submatrix>::computeSigma(
    const SpGreensFunction& G, const SpGreensFunction& G0, SpGreensFunction& Sigma) const {
  const int matrix_dim = Nu::dmn_size();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inv(matrix_dim);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_inv(matrix_dim);
  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  // Sigma = 1/G0 - 1/G
  for (int k_ind = 0; k_ind < KDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < WDmn::dmn_size(); w_ind++) {
      dca::linalg::matrixop::copyArrayToMatrix(matrix_dim, matrix_dim, &G(0, 0, k_ind, w_ind),
                                               matrix_dim, G_inv);
      dca::linalg::matrixop::smallInverse(G_inv, ipiv, work);

      dca::linalg::matrixop::copyArrayToMatrix(matrix_dim, matrix_dim, &G0(0, 0, k_ind, w_ind),
                                               matrix_dim, G0_inv);
      dca::linalg::matrixop::smallInverse(G0_inv, ipiv, work);

      for (int nu2 = 0; nu2 < matrix_dim; ++nu2)
        for (int nu1 = 0; nu1 < matrix_dim; ++nu1)
          Sigma(nu1, nu2, k_ind, w_ind) = G0_inv(nu1, nu2) - G_inv(nu1, nu2);
    }
  }

  symmetrize::execute<Lattice>(data_.Sigma, data_.H_symmetry);
  // TODO : if it is needed implement.
  //   if (parameters_.adjust_self_energy_for_double_counting())
  //    adjust_self_energy_for_double_counting();
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
double CtintClusterSolver<device_t, Parameters, use_submatrix>::L2Difference() const {
  const double alpha = parameters_.get_self_energy_mixing_factor();

  for (int l = 0; l < data_.Sigma.size(); l++)
    data_.Sigma(l) = alpha * data_.Sigma(l) + (1. - alpha) * data_.Sigma_cluster(l);

  const int offset = std::min(1, WDmn::dmn_size() / 2);

  double diff_L2 = 0.;
  for (int w_ind = WDmn::dmn_size() / 2; w_ind < WDmn::dmn_size() / 2 + offset; w_ind++) {
    for (int k_ind = 0; k_ind < KDmn::dmn_size(); k_ind++) {
      for (int l1 = 0; l1 < Nu::dmn_size(); l1++) {
        diff_L2 += std::pow(
            std::abs(data_.Sigma(l1, l1, k_ind, w_ind) - data_.Sigma_cluster(l1, l1, k_ind, w_ind)),
            2);
      }
    }
  }

  double L2_error = std::sqrt(diff_L2) / double(KDmn::dmn_size());
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t |Sigma_QMC - Sigma_cg|_2 ~ " << L2_error << "\n\n";

  return L2_error;
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
double CtintClusterSolver<device_t, Parameters, use_submatrix>::computeDensity() const {
  double result(0.);
  const int t0_minus = TDmn::dmn_size() / 2 - 1;
  for (int i = 0; i < Nu::dmn_size(); i++)
    result += data_.G_r_t(i, i, RDmn::parameter_type::origin_index(), t0_minus);

  return result;
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
double CtintClusterSolver<device_t, Parameters, use_submatrix>::gatherMAndG4(SpGreensFunction& M,
                                                                             bool compute_error) const {
  const auto& M_r = accumulator_.get_sign_times_M_r_w();
  math::transform::FunctionTransform<RDmn, KDmn>::execute(M_r, M);

  double sign = accumulator_.get_accumulated_sign();

  symmetrize::execute<Lattice>(M, data_.H_symmetry);

  // TODO: delay sum.
  auto collect = [&](auto& f) {
    if (compute_error)
      concurrency_.leaveOneOutSum(f);
    else
      concurrency_.sum(f);
  };

  collect(M);
  collect(sign);

  std::size_t n_meas = accumulator_.get_number_of_measurements();
  concurrency_.sum(n_meas);

  M /= std::complex<double>(sign, 0.);

  if (perform_tp_accumulation_) {
    for (int channel = 0; channel < data_.get_G4().size(); ++channel) {
      auto& G4 = data_.get_G4()[channel];
      G4 = accumulator_.get_sign_times_G4()[channel];
      collect(G4);
      G4 /= sign * parameters_.get_beta() * parameters_.get_beta();
    }
  }

  return sign / double(n_meas);
}

template <dca::linalg::DeviceType device_t, class Parameters, bool use_submatrix>
auto CtintClusterSolver<device_t, Parameters, use_submatrix>::local_G_k_w() const {
  const auto& M_r = accumulator_.get_sign_times_M_r_w();
  SpGreensFunction M;
  math::transform::FunctionTransform<RDmn, KDmn>::execute(M_r, M);

  const double sign = accumulator_.get_accumulated_sign();

  M /= sign;
  SpGreensFunction G_k_w("G_k_w");
  computeG_k_w(data_.G0_k_w_cluster_excluded, M, G_k_w);

  return G_k_w;
}

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_CTINT_CLUSTER_SOLVER_HPP
