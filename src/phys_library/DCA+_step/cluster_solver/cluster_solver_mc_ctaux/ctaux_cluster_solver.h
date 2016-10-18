// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Cluster Monte Carlo integrator based on a continuous-time auxilary field expansion.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_CLUSTER_SOLVER_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_CLUSTER_SOLVER_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_template.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/statistics/util.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"

#include "comp_library/linalg/linalg.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_accumulator.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"
#include "phys_library/vertex_measurement_type.hpp"

using namespace dca;
using namespace dca::phys;

namespace DCA {

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type> {
public:
  typedef MOMS_type this_MOMS_type;
  typedef parameters_type this_parameters_type;

  using rng_type = typename parameters_type::random_number_generator;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef QMCI::MC_walker<QMCI::CT_AUX_SOLVER, device_t, parameters_type, MOMS_type> walker_type;
  typedef QMCI::MC_accumulator<QMCI::CT_AUX_SOLVER, dca::linalg::CPU, parameters_type, MOMS_type>
      accumulator_type;

  using w = func::dmn_0<frequency_domain>;
  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using r_DCA = func::dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                           CLUSTER, REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_DCA = func::dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                           CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  using nu_nu_k_DCA_w = func::dmn_variadic<nu, nu, k_DCA, w>;

public:
  cluster_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  ~cluster_solver();

  template <typename Writer>
  void write(Writer& reader);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

  // Returns the function G_k_w before the average across MPI ranks is performed.
  // For testing purposes.
  // TODO: This method duplicates parts of compute_error_bars.
  auto onNode_G_k_w();

protected:
  void warm_up(walker_type& walker);

  void measure(walker_type& walker);

  void update_shell(int i, int N, int N_k);
  void update_shell(int i, int N, int N_k, int N_s);

  void symmetrize_measurements();

  void compute_error_bars();

  // Sums/averages the quantities measured by the individual MPI ranks.
  void collect_measurements();

  void compute_G_k_w_from_M_r_w();

  double compute_S_k_w_from_G_k_w();

  void compute_G_k_w_new(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& M_k_w_new,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& G_k_w_new);

  void compute_S_k_w_new(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& G_k_w_new,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_k_w_new);

  void set_non_interacting_bands_to_zero();

  void adjust_self_energy_for_double_counting();

  double mix_self_energy(double alpha);

protected:
  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

  double thermalization_time;
  double MC_integration_time;

  double total_time;

  rng_type rng;

  accumulator_type accumulator;

  func::function<std::complex<double>, nu_nu_k_DCA_w> Sigma_old;
  func::function<std::complex<double>, nu_nu_k_DCA_w> Sigma_new;

  int DCA_iteration;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::cluster_solver(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thermalization_time(0),
      MC_integration_time(0),

      total_time(0),

      rng(concurrency.id(), concurrency.number_of_processors(), parameters.get_seed()),

      accumulator(parameters, MOMS, 0),

      Sigma_old("Self-Energy-n-1-iteration"),
      Sigma_new("Self-Energy-n-0-iteration"),

      DCA_iteration(-1) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t CT-AUX Integrator is born \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::~cluster_solver() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t CT-AUX Integrator has died \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename Writer>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::write(Writer& writer) {
  writer.open_group("CT-AUX-SOLVER-functions");

  writer.execute(Sigma_old);
  writer.execute(Sigma_new);

  accumulator.write(writer);

  writer.close_group();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::initialize(
    int dca_iteration) {
  DCA_iteration = dca_iteration;

  Sigma_old = MOMS.Sigma;

  accumulator.initialize(DCA_iteration);

  if (concurrency.id() == 0)
    std::cout << "\n\n\t CT-AUX Integrator has initialized (DCA-iteration : " << dca_iteration
              << ")\n\n";
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::integrate() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t integration has started" << std::endl;

  walker_type walker(parameters, MOMS, rng, 0);

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

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t on node integration has ended" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename dca_info_struct_t>
double cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::finalize(
    dca_info_struct_t& dca_info_struct) {
  collect_measurements();
  symmetrize_measurements();

  // Compute new Sigma.
  compute_G_k_w_from_M_r_w();

  // FT<k_DCA,r_DCA>::execute(MOMS.G_k_w, MOMS.G_r_w);
  math::transform::FunctionTransform<k_DCA, r_DCA>::execute(MOMS.G_k_w, MOMS.G_r_w);

  dca_info_struct.L2_Sigma_difference(DCA_iteration) = compute_S_k_w_from_G_k_w();

  for (int i = 0; i < b::dmn_size() * s::dmn_size(); i++) {
    for (int j = 0; j < k_DCA::dmn_size(); j++) {
      std::vector<double> x;
      for (int l = 0; l < w::dmn_size() / 4; l++)
        x.push_back(real(MOMS.Sigma(i, i, j, l)));

      dca_info_struct.Sigma_zero_moment(i, j, DCA_iteration) =
          math::statistics::util::mean(x);  // real(MOMS.Sigma(i,i,j,0));
      dca_info_struct.standard_deviation(i, j, DCA_iteration) =
          math::statistics::util::standard_deviation(x);
    }
  }

  //     if(DCA_iteration == parameters.get_DCA_iterations()-1 &&
  //     parameters.do_equal_time_measurements())
  //       MOMS.G_r_t =

  if (DCA_iteration == parameters.get_DCA_iterations() - 1 &&
      parameters.get_vertex_measurement_type() != NONE)
    MOMS.G4_k_k_w_w /= parameters.get_beta() * parameters.get_beta();

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

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::warm_up(
    walker_type& walker) {
  profiler_type profiler("thermalization", "QMCI", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t warm-up has started\n" << std::endl;

  for (int i = 0; i < parameters.get_warm_up_sweeps(); i++) {
    walker.do_sweep();

    int N_s = walker.get_configuration().size();
    int N_k = walker.get_configuration().get_number_of_interacting_HS_spins();

    update_shell(i, parameters.get_warm_up_sweeps(), N_k, N_s);
  }

  walker.is_thermalized() = true;

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t warm-up has ended\n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::measure(
    walker_type& walker) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t measuring has started \n" << std::endl;

  for (int i = 0; i < parameters.get_number_of_measurements(); i++) {
    {
      profiler_type profiler("updating", "QMCI", __LINE__);
      walker.do_sweep();
    }

    {
      profiler_type profiler("measurements", "QMCI", __LINE__);
      accumulator.update_from(walker);
      accumulator.measure();
    }

    int N_s = walker.get_configuration().size();
    int N_k = walker.get_configuration().get_number_of_interacting_HS_spins();

    update_shell(i, parameters.get_number_of_measurements(), N_k, N_s);
  }

  accumulator.finalize();

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t measuring has ended \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::update_shell(
    int i, int N, int N_k) {
  int tmp = i;

  if (concurrency.id() == concurrency.first() && N > 10 && (tmp % (N / 10)) == 0) {
    std::cout << std::scientific;
    std::cout.precision(6);

    std::cout << "\t\t\t" << double(i) / double(N) * 100. << " % completed \t ";

    std::cout << "\t <k> :" << N_k << "      ";
    std::cout << dca::util::print_time() << "\n";
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::update_shell(
    int i, int N, int N_k, int N_s) {
  int tmp = i;

  if (concurrency.id() == concurrency.first() && N > 10 && (tmp % (N / 10)) == 0) {
    std::cout << std::scientific;
    std::cout.precision(6);

    std::cout << "\t\t\t" << double(i) / double(N) * 100. << " % completed \t ";

    std::cout << "\t <k> :" << N_k << "    N : " << N_s << "      ";
    std::cout << dca::util::print_time() << "\n";
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::compute_error_bars() {
  if (concurrency.id() == 0)
    std::cout << "\n\t\t compute-error-bars on Self-energy\t" << dca::util::print_time() << "\n\n";

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> G_k_w_new("G_k_w_new");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_DCA, w>> M_r_w_new("M_r_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> M_k_w_new("M_k_w_new");

  const int nb_measurements = accumulator.get_number_of_measurements();
  double sign = accumulator.get_sign() / double(nb_measurements);

  for (int l = 0; l < accumulator.get_M_r_w().size(); l++)
    M_r_w_new(l) = accumulator.get_M_r_w()(l) / double(nb_measurements * sign);

  math::transform::FunctionTransform<r_DCA, k_DCA>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);
  compute_S_k_w_new(G_k_w_new, Sigma_new);

  concurrency.average_and_compute_stddev(Sigma_new, MOMS.Sigma_stddev, 1);
  concurrency.average_and_compute_stddev(G_k_w_new, MOMS.G_k_w_stddev, 1);

  // sum G4
  if (parameters.get_vertex_measurement_type() != NONE) {
    if (concurrency.id() == 0)
      std::cout << "\n\t\t compute-error-bars on G4\t" << dca::util::print_time() << "\n\n";

    double sign = accumulator.get_sign() / double(nb_measurements);

    for (int l = 0; l < MOMS.G4_k_k_w_w.size(); l++)
      MOMS.G4_k_k_w_w(l) = accumulator.get_G4()(l) / double(nb_measurements * sign);

    MOMS.G4_k_k_w_w /= parameters.get_beta() * parameters.get_beta();

    concurrency.average_and_compute_stddev(MOMS.G4_k_k_w_w, MOMS.G4_k_k_w_w_stddev, 1);
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::collect_measurements() {
  const int nb_measurements = accumulator.get_number_of_measurements();

  if (concurrency.id() == 0)
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
    concurrency.sum_and_average(accumulator.get_K_r_t(), nb_measurements);
  }

  {
    profiler_type profiler("QMC-self-energy", "QMC-collectives", __LINE__);
    concurrency.sum_and_average(accumulator.get_M_r_w(), nb_measurements);
  }

  {
    profiler_type profiler("QMC-self-energy", "QMC-collectives", __LINE__);
    concurrency.sum_and_average(accumulator.get_M_r_w_squared(), nb_measurements);
  }

  accumulator.get_K_r_t() /= accumulator.get_sign();          // sign;
  accumulator.get_M_r_w() /= accumulator.get_sign();          // sign;
  accumulator.get_M_r_w_squared() /= accumulator.get_sign();  // sign;

  MOMS.K_r_t = accumulator.get_K_r_t();

  if (parameters.do_equal_time_measurements()) {
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

    MOMS.G_r_t = accumulator.get_G_r_t();
  }

  // sum G4
  if (parameters.get_vertex_measurement_type() != NONE) {
    {
      profiler_type profiler("QMC-two-particle-Greens-function", "QMC-collectives", __LINE__);
      concurrency.sum_and_average(accumulator.get_G4(), nb_measurements);
    }

    for (int l = 0; l < MOMS.G4_k_k_w_w.size(); l++)
      MOMS.G4_k_k_w_w(l) = accumulator.get_G4()(l) / accumulator.get_sign();  // sign;
  }

  concurrency.sum(accumulator.get_visited_expansion_order_k());

  concurrency.sum(accumulator.get_error_distribution());
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type,
                    MOMS_type>::symmetrize_measurements() {
  if (concurrency.id() == 0)
    std::cout << "\n\t\t symmetrize measurements has started \t" << dca::util::print_time() << "\n";

  symmetrize::execute(accumulator.get_M_r_w(), MOMS.H_symmetry);
  symmetrize::execute(accumulator.get_M_r_w_squared(), MOMS.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type,
                    MOMS_type>::compute_G_k_w_from_M_r_w() {
  // FT<r_DCA, k_DCA>::execute(accumulator.get_M_r_w(), accumulator.get_M_k_w());
  math::transform::FunctionTransform<r_DCA, k_DCA>::execute(accumulator.get_M_r_w(),
                                                            accumulator.get_M_k_w());

  int matrix_size = b::dmn_size() * s::dmn_size() * b::dmn_size() * s::dmn_size();
  int matrix_dim = b::dmn_size() * s::dmn_size();

  std::complex<double>* G_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_cluster_excluded_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* M_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_times_M_matrix = new std::complex<double>[matrix_size];

  // G = G0 - G0*M*G0/beta

  for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      memset(G_matrix, 0, sizeof(std::complex<double>) * matrix_size);
      memset(G0_times_M_matrix, 0, sizeof(std::complex<double>) * matrix_size);

      memcpy(G0_cluster_excluded_matrix, &MOMS.G0_k_w_cluster_excluded(0, 0, 0, 0, k_ind, w_ind),
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

      memcpy(&MOMS.G_k_w(0, 0, 0, 0, k_ind, w_ind), G_matrix,
             sizeof(std::complex<double>) * matrix_size);
    }
  }

  symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);

  delete[] G_matrix;
  delete[] G0_cluster_excluded_matrix;
  delete[] M_matrix;
  delete[] G0_times_M_matrix;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
double cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type,
                      MOMS_type>::compute_S_k_w_from_G_k_w() {
  static double alpha = parameters.get_DCA_mixing_factor();
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

  for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      dca::linalg::matrixop::copyArrayToMatrix(matrix_dim, matrix_dim,
                                               &MOMS.G_k_w(0, 0, 0, 0, k_ind, w_ind), matrix_dim,
                                               G_inverted_matrix);
      dca::linalg::matrixop::inverse(G_inverted_matrix, ipiv, work);

      dca::linalg::matrixop::copyArrayToMatrix(
          matrix_dim, matrix_dim, &MOMS.G0_k_w_cluster_excluded(0, 0, 0, 0, k_ind, w_ind),
          matrix_dim, G0_cluster_excluded_inverted_matrix);
      dca::linalg::matrixop::inverse(G0_cluster_excluded_inverted_matrix, ipiv, work);

      for (int j = 0; j < sigma_matrix.nrCols(); ++j)
        for (int i = 0; i < sigma_matrix.nrRows(); ++i)
          sigma_matrix(i, j) = G0_cluster_excluded_inverted_matrix(i, j) - G_inverted_matrix(i, j);

      dca::linalg::matrixop::copyMatrixToArray(sigma_matrix, &MOMS.Sigma(0, 0, 0, 0, k_ind, w_ind),
                                               matrix_dim);
    }
  }

  // set_non_interacting_bands_to_zero();

  symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

  if (parameters.adjust_self_energy_for_double_counting())
    adjust_self_energy_for_double_counting();

  double L2_norm = mix_self_energy(alpha);

  return L2_norm;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::compute_G_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& M_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& G_k_w_new) {
  //     if(concurrency.id()==0)
  //       std::cout << "\n\t\t compute-G_k_w_new\t" << dca::util::print_time() << "\n\n";

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix("G_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix("G0_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> M_matrix("M_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_M_matrix("M_matrix", nu::dmn_size());

  for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = MOMS.G0_k_w_cluster_excluded(i, j, k_ind, w_ind);

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

  symmetrize::execute(G_k_w_new, MOMS.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::compute_S_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& G_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_k_w_new) {
  //     if(concurrency.id()==0)
  //       std::cout << "\n\t\t start compute-S_k_w\t" << dca::util::print_time() << "\n\n";

  int N = nu::dmn_size();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix(N);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix(N);

  for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = MOMS.G0_k_w_cluster_excluded(i, j, k_ind, w_ind);

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

  symmetrize::execute(S_k_w_new, MOMS.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type,
                    MOMS_type>::set_non_interacting_bands_to_zero() {
  //  for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
  //    for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
  //      for(int l2=0; l2<b::dmn_size(); l2++){
  //        for(int l1=0; l1<b::dmn_size(); l1++){
  //
  //          if( !(parameters.is_interacting_band()[l1] and
  //                parameters.is_interacting_band()[l2]))
  //            {
  //              MOMS.Sigma(l1,0,l2,0,k_ind, w_ind) = 0.;
  //              MOMS.Sigma(l1,1,l2,1,k_ind, w_ind) = 0.;
  //            }
  //        }
  //      }
  //    }
  //  }
  //
  //  for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //    for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
  //      for(int l2=0; l2<2*b::dmn_size(); l2++)
  //        for(int l1=0; l1<2*b::dmn_size(); l1++)
  //          if( !(l1==l2 and parameters.is_interacting_band()[l1]) )
  //            MOMS.Sigma(l1,l2,k_ind, w_ind) = 0.;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type,
                    MOMS_type>::adjust_self_energy_for_double_counting() {
  set_non_interacting_bands_to_zero();

  //  func::function<double, nu> d_0;
  //  for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
  //    for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
  //      for(int w_ind=0; w_ind<32; w_ind++)
  //        d_0(l1) += real(MOMS.Sigma(l1,l1,k_ind,w_ind));
  //
  //  d_0 /= double(32.*k_DCA::dmn_size());
  //
  //  for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
  //    for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
  //      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //        MOMS.Sigma(l1,l1,k_ind,w_ind) -= d_0(l1);
  //
  //  if(parameters.get_double_counting_method()=="constant")
  //    {
  //      std::vector<int>& interacting_bands = parameters.get_interacting_bands();
  //
  //      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //        for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
  //          for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
  //            for(int b_ind=0; b_ind<interacting_bands.size(); b_ind++)
  //              MOMS.Sigma(interacting_bands[b_ind], s_ind,
  //                         interacting_bands[b_ind], s_ind,
  //                         k_ind                   , w_ind) -=
  //  parameters.get_double_counting_correction();
  //    }
  //
  //  if(parameters.get_double_counting_method()=="adaptive")
  //    {
  //      std::vector<int>& interacting_bands = parameters.get_interacting_bands();
  //
  //      for(int b_ind=0; b_ind<interacting_bands.size(); b_ind++)
  //        for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
  //          for(int s_ind=0; s_ind<s::dmn_size(); s_ind++){
  //
  //            double value = real(MOMS.Sigma(interacting_bands[b_ind], s_ind,
  //                                           interacting_bands[b_ind], s_ind,
  //                                           k_ind                   , 0));
  //
  //            for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
  //
  //              MOMS.Sigma(interacting_bands[b_ind], s_ind,
  //                         interacting_bands[b_ind], s_ind,
  //                         k_ind                   , w_ind) -= value;
  //            }
  //          }
  //        }
  //    }

  symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
double cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::mix_self_energy(
    double alpha) {
  symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);
  symmetrize::execute(MOMS.Sigma_cluster, MOMS.H_symmetry);

  for (int l = 0; l < MOMS.Sigma.size(); l++)
    MOMS.Sigma(l) = alpha * MOMS.Sigma(l) + (1. - alpha) * MOMS.Sigma_cluster(l);

  int offset = std::min(1, w::dmn_size() / 2);

  double L2_norm = 0;
  double diff_L2_norm = 0.;
  for (int w_ind = w::dmn_size() / 2; w_ind < w::dmn_size() / 2 + offset; w_ind++) {
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
      for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++) {
        L2_norm += std::pow(std::abs(MOMS.Sigma(l1, l1, k_ind, w_ind)), 2);
        diff_L2_norm += std::pow(
            std::abs(MOMS.Sigma(l1, l1, k_ind, w_ind) - MOMS.Sigma_cluster(l1, l1, k_ind, w_ind)), 2);
      }
    }
  }

  double error_infty_norm = 0;
  offset = std::min(10, w::dmn_size() / 2);
  for (int w_ind = w::dmn_size() / 2; w_ind < w::dmn_size() / 2 + offset; w_ind++) {
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
      for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++) {
        error_infty_norm = std::max(error_infty_norm, abs(MOMS.Sigma(l1, l1, k_ind, w_ind) -
                                                          MOMS.Sigma_cluster(l1, l1, k_ind, w_ind)));
      }
    }
  }

  double L2_error = std::sqrt(diff_L2_norm) / double(k_DCA::dmn_size());
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\t\t |Sigma_QMC - Sigma_cg|_infty ~ " << error_infty_norm;
    std::cout << "\n\t\t |Sigma_QMC - Sigma_cg|_2 ~ " << L2_error << "\n\n";
  }
  return L2_error;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
auto cluster_solver<CT_AUX_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::onNode_G_k_w() {
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> G_k_w_new("G_k_w");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_DCA, w>> M_r_w_new("M_r_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> M_k_w_new("M_k_w_new");

  const int nb_measurements = accumulator.get_number_of_measurements();
  double sign = accumulator.get_sign() / double(nb_measurements);

  for (int l = 0; l < accumulator.get_M_r_w().size(); l++)
    M_r_w_new(l) = accumulator.get_M_r_w()(l) / double(nb_measurements * sign);

  math::transform::FunctionTransform<r_DCA, k_DCA>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);

  return G_k_w_new;
}

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_CLUSTER_SOLVER_H
