// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
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

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
class CtauxClusterSolver {
public:
  using DataType = Data;
  typedef parameters_type this_parameters_type;

  using rng_type = typename parameters_type::random_number_generator;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef ctaux::CtauxWalker<device_t, parameters_type, Data> walker_type;
  typedef ctaux::CtauxAccumulator<dca::linalg::CPU, parameters_type, Data> accumulator_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  using nu_nu_k_DCA_w = func::dmn_variadic<nu, nu, KClusterDmn, w>;

public:
  CtauxClusterSolver(parameters_type& parameters_ref, Data& MOMS_ref);

  template <typename Writer>
  void write(Writer& reader);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

  // Computes and returns the local value of the Green's function G(k, \omega), i.e. without
  // averaging it across processes.
  // For testing purposes.
  // Precondition: The accumulator data has not been averaged, i.e. finalize has not been called.
  auto local_G_k_w() const;

protected:
  void warm_up(walker_type& walker);

  void measure(walker_type& walker);

  void symmetrize_measurements();

  void compute_error_bars();

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
  parameters_type& parameters;
  Data& data_;
  concurrency_type& concurrency;

  double thermalization_time;
  double MC_integration_time;

  double total_time;

  rng_type rng;

  accumulator_type accumulator;

  func::function<std::complex<double>, nu_nu_k_DCA_w> Sigma_old;
  func::function<std::complex<double>, nu_nu_k_DCA_w> Sigma_new;

  int DCA_iteration;

private:
  bool averaged_;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
CtauxClusterSolver<device_t, parameters_type, Data>::CtauxClusterSolver(parameters_type& parameters_ref,
                                                                        Data& data_ref)
    : parameters(parameters_ref),
      data_(data_ref),
      concurrency(parameters.get_concurrency()),

      thermalization_time(0),
      MC_integration_time(0),

      total_time(0),

      rng(concurrency.id(), concurrency.number_of_processors(), parameters.get_seed()),

      accumulator(parameters, data_, 0),

      Sigma_old("Self-Energy-n-1-iteration"),
      Sigma_new("Self-Energy-n-0-iteration"),

      DCA_iteration(-1),
      averaged_(false) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t CT-AUX Integrator is born \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
template <typename Writer>
void CtauxClusterSolver<device_t, parameters_type, Data>::write(Writer& writer) {
  writer.open_group("CT-AUX-SOLVER-functions");

  writer.execute(Sigma_old);
  writer.execute(Sigma_new);

  accumulator.write(writer);

  writer.close_group();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::initialize(int dca_iteration) {
  DCA_iteration = dca_iteration;

  Sigma_old = data_.Sigma;

  accumulator.initialize(DCA_iteration);

  averaged_ = false;

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t CT-AUX Integrator has initialized (DCA-iteration : " << dca_iteration
              << ")\n\n";
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::integrate() {
  if (concurrency.id() == concurrency.first()) {
    std::cout << "QMC integration has started: " << dca::util::print_time() << std::endl;
  }

  walker_type walker(parameters, data_, rng, 0);

  walker.initialize();

  {
    dca::profiling::WallTime start_time;

    warm_up(walker);

    dca::profiling::WallTime mid_time;

    measure(walker);

    dca::profiling::WallTime end_time;

    dca::profiling::Duration ther_time(mid_time, start_time);
    dca::profiling::Duration meas_time(end_time, mid_time);

    dca::profiling::Duration tot_time(end_time, start_time);

    thermalization_time = ther_time.sec + 1.e-6 * ther_time.usec;
    MC_integration_time = meas_time.sec + 1.e-6 * meas_time.usec;
    total_time = tot_time.sec + 1.e-6 * tot_time.usec;
  }

  accumulator.get_error_distribution() += walker.get_error_distribution();

  if (concurrency.id() == concurrency.first()) {
    std::cout << "On-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters.get_measurements() << std::endl;

    walker.printSummary();
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
template <typename dca_info_struct_t>
double CtauxClusterSolver<device_t, parameters_type, Data>::finalize(dca_info_struct_t& dca_info_struct) {
  collect_measurements();
  symmetrize_measurements();

  // Compute new Sigma.
  compute_G_k_w_from_M_r_w();

  // FT<k_DCA,r_DCA>::execute(data_.G_k_w, data_.G_r_w);
  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(data_.G_k_w, data_.G_r_w);

  dca_info_struct.L2_Sigma_difference(DCA_iteration) = compute_S_k_w_from_G_k_w();

  for (int i = 0; i < b::dmn_size() * s::dmn_size(); i++) {
    for (int j = 0; j < KClusterDmn::dmn_size(); j++) {
      std::vector<double> x;
      for (int l = 0; l < w::dmn_size() / 4; l++)
        x.push_back(real(data_.Sigma(i, i, j, l)));

      dca_info_struct.Sigma_zero_moment(i, j, DCA_iteration) =
          math::statistics::util::mean(x);  // real(data_.Sigma(i,i,j,0));
      dca_info_struct.standard_deviation(i, j, DCA_iteration) =
          math::statistics::util::standard_deviation(x);
    }
  }

  //     if(DCA_iteration == parameters.get_dca_iterations()-1 &&
  //     parameters.additional_time_measurements())
  //       data_.G_r_t =

  if (DCA_iteration == parameters.get_dca_iterations() - 1 && parameters.get_four_point_type() != NONE)
    data_.get_G4() /= parameters.get_beta() * parameters.get_beta();

  double total = 1.e-6, integral = 0;

  for (int l = 0; l < accumulator.get_visited_expansion_order_k().size(); l++) {
    total += accumulator.get_visited_expansion_order_k()(l);
    integral += accumulator.get_visited_expansion_order_k()(l) * l;
  }

  dca_info_struct.average_expansion_order(DCA_iteration) = integral / total;

  dca_info_struct.sign(DCA_iteration) = accumulator.get_sign();

  dca_info_struct.thermalization_per_mpi_task(DCA_iteration) =
      thermalization_time / double(concurrency.number_of_processors());
  dca_info_struct.MC_integration_per_mpi_task(DCA_iteration) =
      MC_integration_time / double(concurrency.number_of_processors());

  dca_info_struct.times_per_mpi_task(DCA_iteration) =
      total_time / double(concurrency.number_of_processors());
  dca_info_struct.Gflop_per_mpi_task(DCA_iteration) =
      accumulator.get_Gflop() / double(concurrency.number_of_processors());

  dca_info_struct.Gflops_per_mpi_task(DCA_iteration) =
      dca_info_struct.Gflop_per_mpi_task(DCA_iteration) /
      dca_info_struct.times_per_mpi_task(DCA_iteration);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t CT-AUX Integrator has finalized \n" << std::endl;

  return dca_info_struct.L2_Sigma_difference(DCA_iteration);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::warm_up(walker_type& walker) {
  profiler_type profiler("thermalization", "QMCI", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t warm-up has started\n" << std::endl;

  for (int i = 0; i < parameters.get_warm_up_sweeps(); i++) {
    walker.do_sweep();
    walker.update_shell(i, parameters.get_warm_up_sweeps());
  }

  walker.is_thermalized() = true;

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t warm-up has ended\n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::measure(walker_type& walker) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t measuring has started \n" << std::endl;

  const int n_meas = parallel::util::getWorkload(parameters.get_measurements(), concurrency);

  for (int i = 0; i < n_meas; i++) {
    {
      profiler_type profiler("updating", "QMCI", __LINE__);
      walker.do_sweep();
    }

    {
      profiler_type profiler("measurements", "QMCI", __LINE__);
      accumulator.update_from(walker);
      accumulator.measure();
    }

    walker.update_shell(i, n_meas);
  }

  accumulator.finalize();

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t measuring has ended \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::compute_error_bars() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t compute-error-bars on Self-energy\t" << dca::util::print_time() << "\n\n";

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w_new(
      "G_k_w_new");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> M_r_w_new(
      "M_r_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> M_k_w_new(
      "M_k_w_new");

  const int nb_measurements = accumulator.get_number_of_measurements();
  double sign = accumulator.get_sign() / double(nb_measurements);

  for (int l = 0; l < accumulator.get_M_r_w().size(); l++)
    M_r_w_new(l) = accumulator.get_M_r_w()(l) / double(nb_measurements * sign);

  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);
  compute_S_k_w_new(G_k_w_new, Sigma_new);

  concurrency.average_and_compute_stddev(Sigma_new, data_.get_Sigma_stdv());
  concurrency.average_and_compute_stddev(G_k_w_new, data_.get_G_k_w_stdv());

  // sum G4
  if (parameters.get_four_point_type() != NONE) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t compute-error-bars on G4\t" << dca::util::print_time() << "\n\n";

    double sign = accumulator.get_sign() / double(nb_measurements);

    auto& G4 = data_.get_G4();
    for (int l = 0; l < G4.size(); l++)
      G4(l) = accumulator.get_G4()(l) / double(nb_measurements * sign);

    G4 /= parameters.get_beta() * parameters.get_beta();

    concurrency.average_and_compute_stddev(G4, data_.get_G4_stdv());
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::collect_measurements() {
  const int nb_measurements = accumulator.get_number_of_measurements();

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t Collect measurements \t" << dca::util::print_time() << "\n"
              << "\n\t\t\t QMC-time : " << total_time << " [sec]"
              << "\n\t\t\t Gflops   : " << accumulator.get_Gflop() / total_time << " [Gf]"
              << "\n\t\t\t sign     : " << accumulator.get_sign() / double(nb_measurements)
              << " \n";

  {
    profiler_type profiler("MC-time", "QMC-collectives", __LINE__);
    concurrency.sum(total_time);
  }

  {  // sum the flops
    profiler_type profiler("MC-flops", "QMC-collectives", __LINE__);
    concurrency.sum(accumulator.get_Gflop());
  }

  {  // sum the sign
    profiler_type profiler("QMC-sign", "QMC-collectives", __LINE__);
    concurrency.sum_and_average(accumulator.get_sign(), nb_measurements);
  }

  // sum M_r_w
  {
    profiler_type profiler("QMC-self-energy", "QMC-collectives", __LINE__);
    concurrency.sum_and_average(accumulator.get_M_r_w(), nb_measurements);
  }

  {
    profiler_type profiler("QMC-self-energy", "QMC-collectives", __LINE__);
    concurrency.sum_and_average(accumulator.get_M_r_w_squared(), nb_measurements);
  }

  accumulator.get_M_r_w() /= accumulator.get_sign();          // sign;
  accumulator.get_M_r_w_squared() /= accumulator.get_sign();  // sign;

  if (parameters.additional_time_measurements()) {
    profiler_type profiler("QMC-two-particle-Greens-function", "QMC-collectives", __LINE__);
    concurrency.sum_and_average(accumulator.get_G_r_t(), nb_measurements);
    concurrency.sum_and_average(accumulator.get_G_r_t_stddev(), nb_measurements);

    accumulator.get_G_r_t() /= accumulator.get_sign();
    accumulator.get_G_r_t_stddev() /= accumulator.get_sign() * std::sqrt(nb_measurements);

    concurrency.sum_and_average(accumulator.get_charge_cluster_moment(), nb_measurements);
    concurrency.sum_and_average(accumulator.get_magnetic_cluster_moment(), nb_measurements);
    concurrency.sum_and_average(accumulator.get_dwave_pp_correlator(), nb_measurements);

    accumulator.get_charge_cluster_moment() /= accumulator.get_sign();
    accumulator.get_magnetic_cluster_moment() /= accumulator.get_sign();
    accumulator.get_dwave_pp_correlator() /= accumulator.get_sign();

    data_.G_r_t = accumulator.get_G_r_t();
  }

  // sum G4
  if (parameters.get_four_point_type() != NONE) {
    {
      profiler_type profiler("QMC-two-particle-Greens-function", "QMC-collectives", __LINE__);
      concurrency.sum_and_average(accumulator.get_G4(), nb_measurements);
    }

    auto& G4 = data_.get_G4();
    for (int l = 0; l < G4.size(); l++)
      G4(l) = accumulator.get_G4()(l) / accumulator.get_sign();  // sign;
  }

  concurrency.sum(accumulator.get_visited_expansion_order_k());

  concurrency.sum(accumulator.get_error_distribution());

  averaged_ = true;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::symmetrize_measurements() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t symmetrize measurements has started \t" << dca::util::print_time() << "\n";

  symmetrize::execute(accumulator.get_M_r_w(), data_.H_symmetry);
  symmetrize::execute(accumulator.get_M_r_w_squared(), data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::compute_G_k_w_from_M_r_w() {
  // FT<RClusterDmn, KClusterDmn>::execute(accumulator.get_M_r_w(), accumulator.get_M_k_w());
  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(accumulator.get_M_r_w(),
                                                                        accumulator.get_M_k_w());

  int matrix_size = b::dmn_size() * s::dmn_size() * b::dmn_size() * s::dmn_size();
  int matrix_dim = b::dmn_size() * s::dmn_size();

  std::complex<double>* G_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_cluster_excluded_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* M_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_times_M_matrix = new std::complex<double>[matrix_size];

  // G = G0 - G0*M*G0/beta

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      memset(G_matrix, 0, sizeof(std::complex<double>) * matrix_size);
      memset(G0_times_M_matrix, 0, sizeof(std::complex<double>) * matrix_size);

      memcpy(G0_cluster_excluded_matrix, &data_.G0_k_w_cluster_excluded(0, 0, 0, 0, k_ind, w_ind),
             sizeof(std::complex<double>) * matrix_size);
      memcpy(M_matrix, &accumulator.get_M_k_w()(0, 0, 0, 0, k_ind, w_ind),
             sizeof(std::complex<double>) * matrix_size);

      // G0 * M --> G0_times_M_matrix
      dca::linalg::blas::gemm("N", "N", matrix_dim, matrix_dim, matrix_dim, 1.,
                              G0_cluster_excluded_matrix, matrix_dim, M_matrix, matrix_dim, 0.,
                              G0_times_M_matrix, matrix_dim);

      // - G0_times_M_matrix * G0 / beta --> G_matrix
      dca::linalg::blas::gemm("N", "N", matrix_dim, matrix_dim, matrix_dim,
                              -1. / parameters.get_beta(), G0_times_M_matrix, matrix_dim,
                              G0_cluster_excluded_matrix, matrix_dim, 0., G_matrix, matrix_dim);

      // G_matrix + G0_cluster_excluded_matrix --> G_matrix
      for (int l = 0; l < matrix_size; l++)
        G_matrix[l] = G_matrix[l] + G0_cluster_excluded_matrix[l];

      memcpy(&data_.G_k_w(0, 0, 0, 0, k_ind, w_ind), G_matrix,
             sizeof(std::complex<double>) * matrix_size);
    }
  }

  symmetrize::execute(data_.G_k_w, data_.H_symmetry);

  delete[] G_matrix;
  delete[] G0_cluster_excluded_matrix;
  delete[] M_matrix;
  delete[] G0_times_M_matrix;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
double CtauxClusterSolver<device_t, parameters_type, Data>::compute_S_k_w_from_G_k_w() {
  static double alpha = parameters.get_self_energy_mixing_factor();
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

  // set_non_interacting_bands_to_zero();

  symmetrize::execute(data_.Sigma, data_.H_symmetry);

  if (parameters.adjust_self_energy_for_double_counting())
    adjust_self_energy_for_double_counting();

  double L2_norm = mix_self_energy(alpha);

  return L2_norm;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::compute_G_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& M_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new) const {
  //     if(concurrency.id()==0)
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
          G_k_w_new(i, j, k_ind, w_ind) = G0_matrix(i, j) - G_matrix(i, j) / parameters.get_beta();
    }
  }

  symmetrize::execute(G_k_w_new, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::compute_S_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& S_k_w_new) {
  //     if(concurrency.id()==0)
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

  if (parameters.adjust_self_energy_for_double_counting())
    adjust_self_energy_for_double_counting();

  //     if(concurrency.id()==0)
  //       std::cout << "\n\t\t end compute-S_k_w\t" << dca::util::print_time() << "\n\n";

  symmetrize::execute(S_k_w_new, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::set_non_interacting_bands_to_zero() {
  //  for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++){
  //      for(int l2=0; l2<b::dmn_size(); l2++){
  //        for(int l1=0; l1<b::dmn_size(); l1++){
  //
  //          if( !(parameters.is_interacting_band()[l1] and
  //                parameters.is_interacting_band()[l2]))
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
  //          if( !(l1==l2 and parameters.is_interacting_band()[l1]) )
  //            data_.Sigma(l1,l2,k_ind, w_ind) = 0.;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxClusterSolver<device_t, parameters_type, Data>::adjust_self_energy_for_double_counting() {
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
  //  if(parameters.get_double_counting_method()=="constant")
  //    {
  //      std::vector<int>& interacting_bands = parameters.get_interacting_orbitals();
  //
  //      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //        for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //          for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
  //            for(int b_ind=0; b_ind<interacting_bands.size(); b_ind++)
  //              data_.Sigma(interacting_bands[b_ind], s_ind,
  //                         interacting_bands[b_ind], s_ind,
  //                         k_ind                   , w_ind) -=
  //  parameters.get_double_counting_correction();
  //    }
  //
  //  if(parameters.get_double_counting_method()=="adaptive")
  //    {
  //      std::vector<int>& interacting_bands = parameters.get_interacting_orbitals();
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

  symmetrize::execute(data_.Sigma, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
double CtauxClusterSolver<device_t, parameters_type, Data>::mix_self_energy(double alpha) {
  symmetrize::execute(data_.Sigma, data_.H_symmetry);
  symmetrize::execute(data_.Sigma_cluster, data_.H_symmetry);

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
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\t\t |Sigma_QMC - Sigma_cg|_infty ~ " << error_infty_norm;
    std::cout << "\n\t\t |Sigma_QMC - Sigma_cg|_2 ~ " << L2_error << "\n\n";
  }
  return L2_error;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
auto CtauxClusterSolver<device_t, parameters_type, Data>::local_G_k_w() const {
  if (averaged_)
    throw std::logic_error("The local data was already averaged.");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w_new(
      "G_k_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> M_k_w_new(
      "M_k_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> M_r_w_new(
      accumulator.get_M_r_w(), "M_r_w_new");

  M_r_w_new /= accumulator.get_sign();

  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);

  return G_k_w_new;
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_CLUSTER_SOLVER_HPP
