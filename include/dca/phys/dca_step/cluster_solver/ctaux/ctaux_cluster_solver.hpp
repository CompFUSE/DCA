// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Cluster Monte Carlo integrator based on a continuous-time auxilary field (CT-AUX) expansion.
//
// TODO: Cleanup the computation of Sigma, error bars, etc. and have the same work flow independent
//       of whether the thread jacket (stdthread qmci) is used.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_CLUSTER_SOLVER_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/config/mc_options.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/statistics/util.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_walker.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
class CtauxClusterSolver {
public:
  using DataType = Data;
  using ParametersType = Parameters;
  using Lattice = typename Parameters::lattice_type;

  using Rng = typename Parameters::random_number_generator;

  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename Parameters::concurrency_type;

  using MCScalar = typename dca::config::McOptions::MCScalar;
  using Walker = ctaux::CtauxWalker<device_t, Parameters, Data, MCScalar>;
  using Accumulator = ctaux::CtauxAccumulator<device_t, Parameters, Data, MCScalar>;

  static constexpr linalg::DeviceType device = device_t;

private:
  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  using NuNuKClusterWDmn = func::dmn_variadic<nu, nu, KClusterDmn, w>;
  using NuNuRClusterWDmn = func::dmn_variadic<nu, nu, RClusterDmn, w>;

public:
  CtauxClusterSolver(Parameters& parameters_ref, Data& MOMS_ref);

  template <typename Writer>
  void write(Writer& writer);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

  // Computes and returns the local value of the Green's function G(k, \omega), i.e. without
  // averaging it across processes.
  // For testing purposes.
  // Precondition: The accumulator_ data has not been averaged, i.e. finalize has not been called.
  auto local_G_k_w() const;

protected:
  void warmUp(Walker& walker);

  void measure(Walker& walker);

  void computeErrorBars();

private:
  void symmetrize_measurements();

  // Sums/averages the quantities measured by the individual MPI ranks.
  void collect_measurements();

  void compute_G_k_w_from_M_r_w();

  double compute_S_k_w_from_G_k_w();

  void compute_G_k_w_new(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& M_k_w_new,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new) const;

  void compute_S_k_w_new(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& S_k_w_new);

  void set_non_interacting_bands_to_zero();

  void adjust_self_energy_for_double_counting();

  double mix_self_energy(double alpha);

protected:
  Parameters& parameters_;
  Data& data_;
  Concurrency& concurrency_;

  double total_time_;

  Accumulator accumulator_;

  int dca_iteration_;

private:
  Rng rng_;

  double thermalization_time_;
  double mc_integration_time_;

  func::function<std::complex<double>, NuNuKClusterWDmn> Sigma_old_;
  func::function<std::complex<double>, NuNuKClusterWDmn> Sigma_new_;

  double accumulated_sign_;
  func::function<std::complex<double>, NuNuRClusterWDmn> M_r_w_;
  func::function<std::complex<double>, NuNuRClusterWDmn> M_r_w_squared_;

  bool averaged_;
  bool compute_jack_knife_;
};

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
CtauxClusterSolver<device_t, Parameters, Data>::CtauxClusterSolver(Parameters& parameters_ref,
                                                                   Data& data_ref)
    : parameters_(parameters_ref),
      data_(data_ref),
      concurrency_(parameters_.get_concurrency()),

      total_time_(0),

      accumulator_(parameters_, data_, 0),

      dca_iteration_(-1),

      rng_(concurrency_.id(), concurrency_.number_of_processors(), parameters_.get_seed()),

      thermalization_time_(0),
      mc_integration_time_(0),

      Sigma_old_("Self-Energy-n-1-iteration"),
      Sigma_new_("Self-Energy-n-0-iteration"),

      M_r_w_("M_r_w"),
      M_r_w_squared_("M_r_w_squared"),

      averaged_(false) {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-AUX Integrator is born \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <typename Writer>
void CtauxClusterSolver<device_t, Parameters, Data>::write(Writer& writer) {
  writer.open_group("CT-AUX-SOLVER-functions");

  writer.execute(Sigma_old_);
  writer.execute(Sigma_new_);

  accumulator_.write(writer);

  writer.close_group();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::initialize(int dca_iteration) {
  dca_iteration_ = dca_iteration;

  Sigma_old_ = data_.Sigma;

  accumulator_.initialize(dca_iteration_);

  averaged_ = false;
  compute_jack_knife_ =
      (dca_iteration == parameters_.get_dca_iterations() - 1) &&
      (parameters_.get_error_computation_type() == ErrorComputationType::JACK_KNIFE);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-AUX Integrator has initialized (DCA-iteration : " << dca_iteration
              << ")\n\n";
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::integrate() {
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "QMC integration has started: " << dca::util::print_time() << std::endl;
  }

  Walker walker(parameters_, data_, rng_, 0);

  walker.initialize();

  {
    dca::profiling::WallTime start_time;

    warmUp(walker);

    dca::profiling::WallTime mid_time;

    measure(walker);

    dca::profiling::WallTime end_time;

    dca::profiling::Duration ther_time(mid_time, start_time);
    dca::profiling::Duration meas_time(end_time, mid_time);

    dca::profiling::Duration tot_time(end_time, start_time);

    thermalization_time_ = ther_time.sec + 1.e-6 * ther_time.usec;
    mc_integration_time_ = meas_time.sec + 1.e-6 * meas_time.usec;
    total_time_ = tot_time.sec + 1.e-6 * tot_time.usec;
  }

  accumulator_.get_error_distribution() += walker.get_error_distribution();

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "On-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters_.get_measurements() << std::endl;

    walker.printSummary();
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <typename dca_info_struct_t>
double CtauxClusterSolver<device_t, Parameters, Data>::finalize(dca_info_struct_t& dca_info_struct) {
  collect_measurements();
  symmetrize_measurements();

  // Compute new Sigma.
  compute_G_k_w_from_M_r_w();

  // FT<k_DCA,r_DCA>::execute(data_.G_k_w, data_.G_r_w);
  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(data_.G_k_w, data_.G_r_w);

  dca_info_struct.L2_Sigma_difference(dca_iteration_) = compute_S_k_w_from_G_k_w();

  for (int i = 0; i < b::dmn_size() * s::dmn_size(); i++) {
    for (int j = 0; j < KClusterDmn::dmn_size(); j++) {
      std::vector<double> x;
      for (int l = 0; l < w::dmn_size() / 4; l++)
        x.push_back(real(data_.Sigma(i, i, j, l)));

      dca_info_struct.Sigma_zero_moment(i, j, dca_iteration_) =
          math::statistics::util::mean(x);  // real(data_.Sigma(i,i,j,0));
      dca_info_struct.standard_deviation(i, j, dca_iteration_) =
          math::statistics::util::standard_deviation(x);
    }
  }

  if (compute_jack_knife_ && parameters_.isAccumulatingG4()) {
    for (std::size_t channel = 0; channel < data_.get_G4_error().size(); ++channel)
      data_.get_G4_error()[channel] = concurrency_.jackknifeError(data_.get_G4()[channel], true);
  }

  double total = 1.e-6, integral = 0;

  for (int l = 0; l < accumulator_.get_visited_expansion_order_k().size(); l++) {
    total += accumulator_.get_visited_expansion_order_k()(l);
    integral += accumulator_.get_visited_expansion_order_k()(l) * l;
  }

  dca_info_struct.average_expansion_order(dca_iteration_) = integral / total;

  dca_info_struct.sign(dca_iteration_) =
      double(accumulator_.get_accumulated_sign()) / accumulator_.get_number_of_measurements();

  dca_info_struct.thermalization_per_mpi_task(dca_iteration_) =
      thermalization_time_ / double(concurrency_.number_of_processors());
  dca_info_struct.MC_integration_per_mpi_task(dca_iteration_) =
      mc_integration_time_ / double(concurrency_.number_of_processors());

  dca_info_struct.times_per_mpi_task(dca_iteration_) =
      total_time_ / double(concurrency_.number_of_processors());
  dca_info_struct.Gflop_per_mpi_task(dca_iteration_) =
      accumulator_.get_Gflop() / double(concurrency_.number_of_processors());

  dca_info_struct.Gflops_per_mpi_task(dca_iteration_) =
      dca_info_struct.Gflop_per_mpi_task(dca_iteration_) /
      dca_info_struct.times_per_mpi_task(dca_iteration_);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-AUX Integrator has finalized \n" << std::endl;

  return dca_info_struct.L2_Sigma_difference(dca_iteration_);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::warmUp(Walker& walker) {
  Profiler profiler("thermalization", "QMCI", __LINE__);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t warm-up has started\n" << std::endl;

  for (int i = 0; i < parameters_.get_warm_up_sweeps(); i++) {
    walker.doSweep();
    walker.updateShell(i, parameters_.get_warm_up_sweeps());
  }

  walker.is_thermalized() = true;

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t warm-up has ended\n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::measure(Walker& walker) {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t measuring has started \n" << std::endl;

  const int n_meas = parallel::util::getWorkload(parameters_.get_measurements(), concurrency_);

  for (int i = 0; i < n_meas; i++) {
    {
      Profiler profiler("updating", "QMCI", __LINE__);
      walker.doSweep();
    }

    {
      Profiler profiler("measurements", "QMCI", __LINE__);
      accumulator_.updateFrom(walker);
      accumulator_.measure();
    }

    walker.updateShell(i, n_meas);
  }

  accumulator_.finalize();

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t measuring has ended \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::computeErrorBars() {
  if (!accumulator_.compute_std_deviation())
    return;
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t compute-error-bars on Self-energy\t" << dca::util::print_time() << "\n\n";

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w_new(
      "G_k_w_new");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> M_r_w_new(
      "M_r_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> M_k_w_new(
      "M_k_w_new");

  accumulator_.finalize();

  M_r_w_new = accumulator_.get_sign_times_M_r_w();
  M_r_w_new /= accumulator_.get_accumulated_sign();

  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);
  compute_S_k_w_new(G_k_w_new, Sigma_new_);

  concurrency_.average_and_compute_stddev(Sigma_new_, data_.get_Sigma_stdv());
  concurrency_.average_and_compute_stddev(G_k_w_new, data_.get_G_k_w_stdv());

  // sum G4
  if (dca_iteration_ == parameters_.get_dca_iterations() - 1 && parameters_.isAccumulatingG4()) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t compute-error-bars on G4\t" << dca::util::print_time() << "\n\n";

    // This creates a copy!
    auto G4 = accumulator_.get_sign_times_G4();

    for (std::size_t channel = 0; channel < G4.size(); ++channel) {
      G4[channel] /=
          parameters_.get_beta() * parameters_.get_beta() * accumulator_.get_accumulated_sign();
      concurrency_.average_and_compute_stddev(G4[channel], data_.get_G4_stdv()[channel]);
    }
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::collect_measurements() {
  auto collect_delayed = [&](auto& f) {
    if (compute_jack_knife_)
      concurrency_.leaveOneOutSum(f, true);
    else
      concurrency_.delayedSum(f);
  };

  const double local_time = total_time_;
  const bool accumulate_g4 =
      parameters_.isAccumulatingG4() && dca_iteration_ == parameters_.get_dca_iterations() - 1;

  {
    Profiler profiler("QMC-collectives", "CT-AUX solver", __LINE__);
    concurrency_.delayedSum(total_time_);
    concurrency_.delayedSum(accumulator_.get_Gflop());
    accumulated_sign_ = accumulator_.get_accumulated_sign();
    collect_delayed(accumulated_sign_);

    M_r_w_ = accumulator_.get_sign_times_M_r_w();
    collect_delayed(M_r_w_);

    if (accumulator_.compute_std_deviation()) {
      M_r_w_squared_ = accumulator_.get_sign_times_M_r_w_sqr();
      concurrency_.delayedSum(M_r_w_squared_);
    }

    if (parameters_.additional_time_measurements()) {
      Profiler profiler("Additional time measurements.", "QMC-collectives", __LINE__);
      concurrency_.delayedSum(accumulator_.get_G_r_t());
      concurrency_.delayedSum(accumulator_.get_G_r_t_stddev());
      concurrency_.delayedSum(accumulator_.get_charge_cluster_moment());
      concurrency_.delayedSum(accumulator_.get_magnetic_cluster_moment());
      concurrency_.delayedSum(accumulator_.get_dwave_pp_correlator());
    }

    // sum G4
    if (accumulate_g4) {
      for (int channel = 0; channel < data_.get_G4().size(); ++channel) {
        auto& G4 = data_.get_G4()[channel];
        // function operator = will reset this G4 size to other G4 size if they are not equal
        G4 = accumulator_.get_sign_times_G4()[channel];
        if (compute_jack_knife_)
          concurrency_.leaveOneOutSum(G4);
        else {
#ifdef DCA_WITH_NVLINK
            concurrency_.gatherv(G4, concurrency_.first());
#else
            concurrency_.localSum(G4, concurrency_.first());
#endif
        }
      }
    }

    concurrency_.delayedSum(accumulator_.get_visited_expansion_order_k());
    concurrency_.delayedSum(accumulator_.get_error_distribution());

    concurrency_.resolveSums();
  }

  M_r_w_ /= accumulated_sign_;
  M_r_w_squared_ /= accumulated_sign_;
  if (accumulate_g4) {
    for (auto& G4 : data_.get_G4())
      G4 /= accumulated_sign_ * parameters_.get_beta() * parameters_.get_beta();
  }

  if (parameters_.additional_time_measurements()) {
    accumulator_.get_G_r_t() /= accumulated_sign_;
    data_.G_r_t = accumulator_.get_G_r_t();
    accumulator_.get_G_r_t_stddev() /= accumulated_sign_ * std::sqrt(parameters_.get_measurements());
    accumulator_.get_charge_cluster_moment() /= accumulated_sign_;
    accumulator_.get_magnetic_cluster_moment() /= accumulated_sign_;
    accumulator_.get_dwave_pp_correlator() /= accumulated_sign_;
  }

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t Collect measurements \t" << dca::util::print_time() << "\n"
              << "\n\t\t\t QMC-local-time : " << local_time << " [sec]"
              << "\n\t\t\t QMC-total-time : " << total_time_ << " [sec]"
              << "\n\t\t\t Gflop   : " << accumulator_.get_Gflop() << " [Gf]"
              << "\n\t\t\t Gflop/s   : " << accumulator_.get_Gflop() / local_time << " [Gf/s]"
              << "\n\t\t\t sign     : " << accumulated_sign_ / parameters_.get_measurements()
              << " \n";

  averaged_ = true;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::symmetrize_measurements() {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t symmetrize measurements has started \t" << dca::util::print_time() << "\n";

  symmetrize::execute<Lattice>(M_r_w_, data_.H_symmetry);
  symmetrize::execute<Lattice>(M_r_w_squared_, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::compute_G_k_w_from_M_r_w() {
  func::function<std::complex<double>, NuNuKClusterWDmn> M_k_w;
  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_, M_k_w);

  int matrix_size = b::dmn_size() * s::dmn_size() * b::dmn_size() * s::dmn_size();
  int matrix_dim = b::dmn_size() * s::dmn_size();

  std::complex<double>* G_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_cluster_excluded_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* M_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_times_M_matrix = new std::complex<double>[matrix_size];

  // G = G0 - G0*M*G0/beta

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      memset(static_cast<void*>(G_matrix), 0, sizeof(std::complex<double>) * matrix_size);
      memset(static_cast<void*>(G0_times_M_matrix), 0, sizeof(std::complex<double>) * matrix_size);

      memcpy(G0_cluster_excluded_matrix, &data_.G0_k_w_cluster_excluded(0, 0, 0, 0, k_ind, w_ind),
             sizeof(std::complex<double>) * matrix_size);
      memcpy(M_matrix, &M_k_w(0, 0, 0, 0, k_ind, w_ind), sizeof(std::complex<double>) * matrix_size);

      // G0 * M --> G0_times_M_matrix
      dca::linalg::blas::gemm("N", "N", matrix_dim, matrix_dim, matrix_dim, 1.,
                              G0_cluster_excluded_matrix, matrix_dim, M_matrix, matrix_dim, 0.,
                              G0_times_M_matrix, matrix_dim);

      // - G0_times_M_matrix * G0 / beta --> G_matrix
      dca::linalg::blas::gemm("N", "N", matrix_dim, matrix_dim, matrix_dim,
                              -1. / parameters_.get_beta(), G0_times_M_matrix, matrix_dim,
                              G0_cluster_excluded_matrix, matrix_dim, 0., G_matrix, matrix_dim);

      // G_matrix + G0_cluster_excluded_matrix --> G_matrix
      for (int l = 0; l < matrix_size; l++)
        G_matrix[l] = G_matrix[l] + G0_cluster_excluded_matrix[l];

      memcpy(&data_.G_k_w(0, 0, 0, 0, k_ind, w_ind), G_matrix,
             sizeof(std::complex<double>) * matrix_size);
    }
  }

  symmetrize::execute<Lattice>(data_.G_k_w, data_.H_symmetry);

  delete[] G_matrix;
  delete[] G0_cluster_excluded_matrix;
  delete[] M_matrix;
  delete[] G0_times_M_matrix;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
double CtauxClusterSolver<device_t, Parameters, Data>::compute_S_k_w_from_G_k_w() {
  static double alpha = parameters_.get_self_energy_mixing_factor();
  //     double L2_difference_norm = 0;
  //     double L2_Sigma_norm      = 0;

  int matrix_dim = b::dmn_size() * s::dmn_size();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inverted_matrix(matrix_dim);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_cluster_excluded_inverted_matrix(
      matrix_dim);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> sigma_matrix(matrix_dim);

  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  // Sigma = 1/G0 - 1/G

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      dca::linalg::matrixop::copyArrayToMatrix(matrix_dim, matrix_dim,
                                               &data_.G_k_w(0, 0, 0, 0, k_ind, w_ind), matrix_dim,
                                               G_inverted_matrix);
      dca::linalg::matrixop::inverse(G_inverted_matrix, ipiv, work);

      dca::linalg::matrixop::copyArrayToMatrix(
          matrix_dim, matrix_dim, &data_.G0_k_w_cluster_excluded(0, 0, 0, 0, k_ind, w_ind),
          matrix_dim, G0_cluster_excluded_inverted_matrix);
      dca::linalg::matrixop::inverse(G0_cluster_excluded_inverted_matrix, ipiv, work);

      for (int j = 0; j < sigma_matrix.nrCols(); ++j)
        for (int i = 0; i < sigma_matrix.nrRows(); ++i)
          sigma_matrix(i, j) = G0_cluster_excluded_inverted_matrix(i, j) - G_inverted_matrix(i, j);

      dca::linalg::matrixop::copyMatrixToArray(sigma_matrix, &data_.Sigma(0, 0, 0, 0, k_ind, w_ind),
                                               matrix_dim);
    }
  }

  // Compute error on G and Self Energy.
  if (compute_jack_knife_) {
    data_.get_G_k_w_error() = concurrency_.jackknifeError(data_.G_k_w, true);
    data_.get_G_r_w_error() = concurrency_.jackknifeError(data_.G_r_w, true);
    data_.get_Sigma_error() = concurrency_.jackknifeError(data_.Sigma, true);
  }

  // set_non_interacting_bands_to_zero();

  symmetrize::execute<Lattice>(data_.Sigma, data_.H_symmetry);

  if (parameters_.adjust_self_energy_for_double_counting())
    adjust_self_energy_for_double_counting();

  double L2_norm = mix_self_energy(alpha);

  return L2_norm;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::compute_G_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& M_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new) const {
  //     if(concurrency_.id()==0)
  //       std::cout << "\n\t\t compute-G_k_w_new\t" << dca::util::print_time() << "\n\n";

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix("G_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix("G0_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> M_matrix("M_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_M_matrix("M_matrix", nu::dmn_size());

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = data_.G0_k_w_cluster_excluded(i, j, k_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          M_matrix(i, j) = M_k_w_new(i, j, k_ind, w_ind);

      dca::linalg::matrixop::gemm(G0_matrix, M_matrix, G0_M_matrix);
      dca::linalg::matrixop::gemm(G0_M_matrix, G0_matrix, G_matrix);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_k_w_new(i, j, k_ind, w_ind) = G0_matrix(i, j) - G_matrix(i, j) / parameters_.get_beta();
    }
  }

  symmetrize::execute<Lattice>(G_k_w_new, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::compute_S_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& S_k_w_new) {
  //     if(concurrency_.id()==0)
  //       std::cout << "\n\t\t start compute-S_k_w\t" << dca::util::print_time() << "\n\n";

  int N = nu::dmn_size();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix(N);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix(N);

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = data_.G0_k_w_cluster_excluded(i, j, k_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_matrix(i, j) = G_k_w_new(i, j, k_ind, w_ind);

      dca::linalg::matrixop::inverse(G_matrix);
      dca::linalg::matrixop::inverse(G0_matrix);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          S_k_w_new(i, j, k_ind, w_ind) = G0_matrix(i, j) - G_matrix(i, j);
    }
  }

  if (parameters_.adjust_self_energy_for_double_counting())
    adjust_self_energy_for_double_counting();

  //     if(concurrency_.id()==0)
  //       std::cout << "\n\t\t end compute-S_k_w\t" << dca::util::print_time() << "\n\n";

  symmetrize::execute<Lattice>(S_k_w_new, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::set_non_interacting_bands_to_zero() {
  //  for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++){
  //      for(int l2=0; l2<b::dmn_size(); l2++){
  //        for(int l1=0; l1<b::dmn_size(); l1++){
  //
  //          if( !(parameters_.is_interacting_band()[l1] and
  //                parameters_.is_interacting_band()[l2]))
  //            {
  //              data_.Sigma(l1,0,l2,0,k_ind, w_ind) = 0.;
  //              data_.Sigma(l1,1,l2,1,k_ind, w_ind) = 0.;
  //            }
  //        }
  //      }
  //    }
  //  }
  //
  //  for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //      for(int l2=0; l2<2*b::dmn_size(); l2++)
  //        for(int l1=0; l1<2*b::dmn_size(); l1++)
  //          if( !(l1==l2 and parameters_.is_interacting_band()[l1]) )
  //            data_.Sigma(l1,l2,k_ind, w_ind) = 0.;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::adjust_self_energy_for_double_counting() {
  set_non_interacting_bands_to_zero();

  //  func::function<double, nu> d_0;
  //  for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //      for(int w_ind=0; w_ind<32; w_ind++)
  //        d_0(l1) += real(data_.Sigma(l1,l1,k_ind,w_ind));
  //
  //  d_0 /= double(32.*KClusterDmn::dmn_size());
  //
  //  for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //        data_.Sigma(l1,l1,k_ind,w_ind) -= d_0(l1);
  //
  //  if(parameters_.get_double_counting_method()=="constant")
  //    {
  //      std::vector<int>& interacting_bands = parameters_.get_interacting_orbitals();
  //
  //      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //        for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //          for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
  //            for(int b_ind=0; b_ind<interacting_bands.size(); b_ind++)
  //              data_.Sigma(interacting_bands[b_ind], s_ind,
  //                         interacting_bands[b_ind], s_ind,
  //                         k_ind                   , w_ind) -=
  //  parameters_.get_double_counting_correction();
  //    }
  //
  //  if(parameters_.get_double_counting_method()=="adaptive")
  //    {
  //      std::vector<int>& interacting_bands = parameters_.get_interacting_orbitals();
  //
  //      for(int b_ind=0; b_ind<interacting_bands.size(); b_ind++)
  //        for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++){
  //          for(int s_ind=0; s_ind<s::dmn_size(); s_ind++){
  //
  //            double value = real(data_.Sigma(interacting_bands[b_ind], s_ind,
  //                                           interacting_bands[b_ind], s_ind,
  //                                           k_ind                   , 0));
  //
  //            for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
  //
  //              data_.Sigma(interacting_bands[b_ind], s_ind,
  //                         interacting_bands[b_ind], s_ind,
  //                         k_ind                   , w_ind) -= value;
  //            }
  //          }
  //        }
  //    }

  symmetrize::execute<Lattice>(data_.Sigma, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
double CtauxClusterSolver<device_t, Parameters, Data>::mix_self_energy(double alpha) {
  symmetrize::execute<Lattice>(data_.Sigma, data_.H_symmetry);
  symmetrize::execute<Lattice>(data_.Sigma_cluster, data_.H_symmetry);

  for (int l = 0; l < data_.Sigma.size(); l++)
    data_.Sigma(l) = alpha * data_.Sigma(l) + (1. - alpha) * data_.Sigma_cluster(l);

  int offset = std::min(1, w::dmn_size() / 2);

  double L2_norm = 0;
  double diff_L2_norm = 0.;
  for (int w_ind = w::dmn_size() / 2; w_ind < w::dmn_size() / 2 + offset; w_ind++) {
    for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
      for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++) {
        L2_norm += std::pow(std::abs(data_.Sigma(l1, l1, k_ind, w_ind)), 2);
        diff_L2_norm += std::pow(
            std::abs(data_.Sigma(l1, l1, k_ind, w_ind) - data_.Sigma_cluster(l1, l1, k_ind, w_ind)),
            2);
      }
    }
  }

  double error_infty_norm = 0;
  offset = std::min(10, w::dmn_size() / 2);
  for (int w_ind = w::dmn_size() / 2; w_ind < w::dmn_size() / 2 + offset; w_ind++) {
    for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
      for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++) {
        error_infty_norm = std::max(error_infty_norm, abs(data_.Sigma(l1, l1, k_ind, w_ind) -
                                                          data_.Sigma_cluster(l1, l1, k_ind, w_ind)));
      }
    }
  }

  double L2_error = std::sqrt(diff_L2_norm) / double(KClusterDmn::dmn_size());
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\n\t\t |Sigma_QMC - Sigma_cg|_infty ~ " << error_infty_norm;
    std::cout << "\n\t\t |Sigma_QMC - Sigma_cg|_2 ~ " << L2_error << "\n\n";
  }
  return L2_error;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
auto CtauxClusterSolver<device_t, Parameters, Data>::local_G_k_w() const {
  if (averaged_)
    throw std::logic_error("The local data was already averaged.");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w_new(
      "G_k_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> M_k_w_new(
      "M_k_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> M_r_w_new(
      accumulator_.get_sign_times_M_r_w(), "M_r_w_new");

  M_r_w_new /= accumulator_.get_accumulated_sign();

  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);

  return G_k_w_new;
}

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_CLUSTER_SOLVER_HPP
