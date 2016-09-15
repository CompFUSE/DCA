// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Bart Ydens
//         Peter Staar (taa@zurich.ibm.com)
//
// Single-site Monte Carlo integrator based on a hybridization expansion.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_SOLVER_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_SOLVER_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_template.h"

#include <cassert>
#include <complex>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "dca/util/print_time.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg_device_types.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_accumulator.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_solver_routines.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_walker.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"

namespace DCA {

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>
    : public QMCI::ss_hybridization_solver_routines<parameters_type, MOMS_type> {
public:
  typedef MOMS_type this_MOMS_type;
  typedef parameters_type this_parameters_type;

  typedef typename parameters_type::random_number_generator rng_type;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef QMCI::ss_hybridization_solver_routines<parameters_type, MOMS_type>
      ss_hybridization_solver_routines_type;

  typedef QMCI::MC_walker<QMCI::SS_CT_HYB, dca::linalg::CPU, parameters_type, MOMS_type> walker_type;
  typedef QMCI::MC_accumulator<QMCI::SS_CT_HYB, dca::linalg::CPU, parameters_type, MOMS_type> accumulator_type;

  using w = dmn_0<frequency_domain>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using r_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  using nu_nu_k_DCA_w = dmn_variadic<nu, nu, k_DCA, w>;

  const static int MC_TYPE = SS_CT_HYB;

public:
  cluster_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  ~cluster_solver();

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

  void read(std::string filename);

  void write(std::string filename);

  template <IO::FORMAT DATA_FORMAT>
  void read(IO::reader<DATA_FORMAT>& reader);

  template <IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& reader);

  // For testing purposes:
  // TODO: Const correctness.
  FUNC_LIB::function<std::complex<double>, dmn_variadic<nu, nu, r_DCA, w>>& get_GS_r_w() {
    return accumulator.get_GS_r_w();
  }

protected:
  void warm_up(walker_type& walker);

  void measure(walker_type& walker);

  void update_shell(int i, int N, int k);

  void symmetrize_measurements();

  void compute_error_bars();

  void sum_measurements();

  void compute_G_k_w();

  double compute_S_k_w_from_G_k_w();

  void measure_Sigma();

  void compute_Sigma_new(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>>& G_r_w,
                         FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>>& GS_r_w);

  int find_w_cutoff();

  void find_tail_of_Sigma(double& S0, double& S1, int b, int s, int k);

  void adjust_self_energy_for_double_counting();

protected:
  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

  double thermalization_time;
  double MC_integration_time;

  double total_time;

  rng_type rng;
  accumulator_type accumulator;

  FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w> Sigma_old;
  FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w> Sigma_new;

  int DCA_iteration;

  FUNC_LIB::function<double, nu> mu_DC;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::cluster_solver(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref)
    : QMCI::ss_hybridization_solver_routines<parameters_type, MOMS_type>(parameters_ref, MOMS_ref),

      parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thermalization_time(0),
      MC_integration_time(0),

      total_time(0),

      rng(concurrency.id(), concurrency.number_of_processors(), parameters.get_seed()),
      accumulator(parameters, MOMS),

      Sigma_old("Self-Energy-n-1-iteration"),
      Sigma_new("Self-Energy-n-0-iteration"),

      DCA_iteration(-1) {
  concurrency << "\n\n\t SS CT-HYB Integrator is born \n\n";
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::~cluster_solver() {
  concurrency << "\n\n\t SS CT-HYB Integrator has died \n\n";
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <IO::FORMAT DATA_FORMAT>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::write(
    IO::writer<DATA_FORMAT>& writer) {
  writer.open_group("SS-HYB-SOLVER-functions");

  writer.execute(this->get_mu());
  writer.execute(this->get_mu_HALF());

  writer.execute(this->get_a0());
  writer.execute(this->get_a1());

  writer.execute(this->get_F_k_w());
  writer.execute(this->get_F_r_t());

  writer.execute(Sigma_old);
  writer.execute(Sigma_new);

  accumulator.write(writer);

  writer.close_group();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::initialize(int dca_iteration) {
  if (concurrency.id() == 0)
    std::cout << "\n\n\t SS CT-HYB Integrator has started ( DCA-iteration : " << dca_iteration
              << ")\n\n";

  DCA_iteration = dca_iteration;

  Sigma_old = MOMS.Sigma_cluster;

  ss_hybridization_solver_routines_type::initialize_functions();

  accumulator.initialize(dca_iteration);

  if (concurrency.id() == 0) {
    std::stringstream ss;
    ss.precision(6);
    ss << std::scientific;

    FUNC_LIB::function<double, nu>& mu = this->get_mu();
    FUNC_LIB::function<double, nu>& mu_half = this->get_mu_HALF();

    FUNC_LIB::function<double, nu>& a0 = this->get_a0();
    FUNC_LIB::function<double, nu>& a1 = this->get_a1();

    ss << "\n\n mu, mu_half, a0, a1\n\n";
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        ss << b_ind << "\t" << s_ind << "\t" << mu(b_ind, s_ind) << "\t" << mu_half(b_ind, s_ind)
           << "\t" << a0(b_ind, s_ind) << "\t" << a1(b_ind, s_ind) << "\n";

    ss << "\n\n";

    std::cout << ss.str();
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::integrate() {
  concurrency << "\n\t\t integration has started \n";

  walker_type walker(parameters, MOMS, rng);

  walker.initialize();

  warm_up(walker);

  measure(walker);

  concurrency << "\n\t\t on node integration has ended \n";
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename dca_info_struct_t>
double cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::finalize(
    dca_info_struct_t& dca_info_struct) {
  // Average measurements over nodes.
  sum_measurements();

  symmetrize_measurements();

  compute_G_k_w();

  math_algorithms::functional_transforms::TRANSFORM<k_DCA, r_DCA>::execute(MOMS.G_k_w, MOMS.G_r_w);

  dca_info_struct.L2_Sigma_difference(DCA_iteration) = compute_S_k_w_from_G_k_w();

  // SHOW::execute_on_bands(accumulator.get_G_r_w());
  // SHOW::execute_on_bands(accumulator.get_GS_r_w());
  // SHOW::execute_on_bands(MOMS.G_k_w);
  // SHOW::execute_on_bands(MOMS.Sigma);

  if (concurrency.id() == 0) {
    std::stringstream ss;
    ss.precision(6);
    ss << std::scientific;

    ss << "\n\n Sigma \n\n";
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++) {
      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
        double result = 0;
        for (int w_ind = 0; w_ind < 50; w_ind++)
          result += real(MOMS.Sigma(b_ind, s_ind, b_ind, s_ind, 0, w_ind)) / 50.;

        ss << b_ind << "\t" << s_ind << "\t" << result << "\n";
      }
    }
    ss << "\n\n";

    std::cout << ss.str();
  }

  for (int i = 0; i < b::dmn_size() * s::dmn_size(); i++)
    for (int j = 0; j < k_DCA::dmn_size(); j++)
      dca_info_struct.Sigma_zero_moment(i, j, DCA_iteration) = real(MOMS.Sigma(i, i, j, 0));

  double total = 1.e-6, integral = 0;
  for (int l = 0; l < accumulator.get_visited_expansion_order_k().size(); l++) {
    total += accumulator.get_visited_expansion_order_k()(l);
    integral += accumulator.get_visited_expansion_order_k()(l) * l;
  }

  dca_info_struct.average_expansion_order(DCA_iteration) = integral / total;

  dca_info_struct.sign(DCA_iteration) = accumulator.get_sign();

  concurrency << "\n\n\t SS CT-HYB Integrator has finalized \n\n";

  return dca_info_struct.L2_Sigma_difference(DCA_iteration);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::warm_up(walker_type& walker) {
  concurrency << "\n\t\t warm-up has started\n\n";

  for (int i = 0; i < parameters.get_warm_up_sweeps(); i++) {
    walker.do_sweep();

    update_shell(i, parameters.get_warm_up_sweeps(), walker.get_configuration().size());
  }

  walker.is_thermalized() = true;

  concurrency << "\n\t\t warm-up has ended\n\n";
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::measure(walker_type& walker) {
  concurrency << "\n\t\t measuring has started \n\n";

  for (int i = 0; i < parameters.get_number_of_measurements(); i++) {
    walker.do_sweep();

    accumulator.update_from(walker);

    accumulator.measure();

    update_shell(i, parameters.get_number_of_measurements(), walker.get_configuration().size());
  }

  // here we need to do a correction a la Andrey

  //  FUNC_LIB::function<double, nu> correction_to_GS;

  // for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
  // for(int b_ind=0; b_ind<b::dmn_size(); b_ind++)
  // correction_to_GS(b_ind,s_ind) = 0//parameters.get_chemical_potential_DC()
  // + mu_DC(b_ind, s_ind)
  // + walker.mu_HALF(b_ind,s_ind)
  // + walker.a0(b_ind,s_ind);

  accumulator.finalize();

  concurrency << "\n\t\t measuring has ended \n\n";
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::update_shell(int i, int N,
                                                                                   int k) {
  if (concurrency.id() == concurrency.first() && N > 10 && (i % (N / 10)) == 0) {
    std::cout << std::scientific;
    std::cout.precision(6);

    std::cout << "\t\t\t" << double(i) / double(N) * 100. << " % completed \t ";

    std::cout << "\t <k> :" << k << " \t";
    std::cout << dca::util::print_time() << "\n";
  }
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::compute_error_bars() {
  concurrency << "\n\t\t computing the error-bars \n";

  const int nb_measurements = accumulator.get_number_of_measurements();
  double sign = accumulator.get_sign() / double(nb_measurements);

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>> G_r_w("G_r_w_tmp");
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>> GS_r_w("GS_r_w_tmp");

  for (int l = 0; l < G_r_w.size(); l++)
    G_r_w(l) = accumulator.get_G_r_w()(l) / double(nb_measurements * sign);

  for (int l = 0; l < GS_r_w.size(); l++)
    GS_r_w(l) = accumulator.get_GS_r_w()(l) / double(nb_measurements * sign);

  compute_Sigma_new(G_r_w, GS_r_w);

  concurrency.average_and_compute_stddev(Sigma_new, MOMS.Sigma_stddev, 1);
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::sum_measurements() {
  concurrency << "\n\t\t sum measurements \n";

  const int nb_measurements = accumulator.get_number_of_measurements();
  
  // sum the sign
  concurrency.sum_and_average(accumulator.get_sign(), nb_measurements);

  // sum G_r_w
  concurrency.sum_and_average(accumulator.get_G_r_w(), nb_measurements);
  accumulator.get_G_r_w() /= accumulator.get_sign();

  // sum GS_r_w
  concurrency.sum_and_average(accumulator.get_GS_r_w(), nb_measurements);
  accumulator.get_GS_r_w() /= accumulator.get_sign();

  concurrency.sum(accumulator.get_visited_expansion_order_k());
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::symmetrize_measurements() {
  concurrency << "\n\t\t symmetrize measurements has started \n";

  symmetrize::execute(accumulator.get_G_r_w(), MOMS.H_symmetry);

  symmetrize::execute(accumulator.get_GS_r_w(), MOMS.H_symmetry);

  std::vector<int> flavors = parameters_type::model_type::get_flavors();
  assert(flavors.size() == b::dmn_size());

  FUNC_LIB::function<std::complex<double>, b> f_val;
  FUNC_LIB::function<std::complex<double>, b> f_tot;

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++) {
      f_val = 0;
      f_tot = 0;

      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
        f_tot(flavors[b_ind]) += 1;
        f_val(flavors[b_ind]) += accumulator.get_G_r_w()(b_ind, s_ind, b_ind, s_ind, 0, w_ind);
      }

      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        accumulator.get_G_r_w()(b_ind, s_ind, b_ind, s_ind, 0, w_ind) =
            f_val(flavors[b_ind]) / f_tot(flavors[b_ind]);
    }
  }

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++) {
      f_val = 0;
      f_tot = 0;

      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
        f_tot(flavors[b_ind]) += 1;
        f_val(flavors[b_ind]) += accumulator.get_GS_r_w()(b_ind, s_ind, b_ind, s_ind, 0, w_ind);
      }

      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        accumulator.get_GS_r_w()(b_ind, s_ind, b_ind, s_ind, 0, w_ind) =
            f_val(flavors[b_ind]) / f_tot(flavors[b_ind]);
    }
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::compute_G_k_w() {
  math_algorithms::functional_transforms::TRANSFORM<r_DCA, k_DCA>::execute(accumulator.get_G_r_w(),
                                                                           MOMS.G_k_w);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
double cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::compute_S_k_w_from_G_k_w() {
  double alpha = DCA_iteration > 0 ? parameters.get_DCA_mixing_factor() : 1;

  double L2_difference_norm = 1.e-6;
  double L2_Sigma_norm = 1.e-6;

  int w_cutoff = find_w_cutoff();

  // measure_Sigma();

  compute_Sigma_new(accumulator.get_G_r_w(), accumulator.get_GS_r_w());

  symmetrize::execute(Sigma_new, MOMS.H_symmetry);

  for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
    if (ss_hybridization_solver_routines_type::is_interacting_band(b_ind)) {
      for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++) {
        for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
          double Sigma_0, Sigma_1;
          find_tail_of_Sigma(Sigma_0, Sigma_1, b_ind, s_ind, k_ind);
          for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
            std::complex<double> new_sigma;
            if (w_ind >= w::dmn_size() / 2 - w_cutoff && w_ind < w::dmn_size() / 2 + w_cutoff)
              new_sigma = Sigma_new(b_ind, s_ind, b_ind, s_ind, k_ind, w_ind);
            else
              new_sigma =
                  std::complex<double>(Sigma_0, Sigma_1 / w::parameter_type::get_elements()[w_ind]);

            std::complex<double> old_sigma = Sigma_old(b_ind, s_ind, b_ind, s_ind, k_ind, w_ind);

            if (w::dmn_size() / 2 - 16 < w_ind and w_ind < w::dmn_size() / 2 + 16) {
              L2_Sigma_norm += imag(new_sigma) * imag(new_sigma);
              L2_difference_norm += imag(old_sigma - new_sigma) * imag(old_sigma - new_sigma);
            }

            MOMS.Sigma(b_ind, s_ind, b_ind, s_ind, k_ind, w_ind) =
                alpha * (new_sigma) + (1 - alpha) * old_sigma;
          }
        }
      }
    }
  }

  symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t |Sigma_old-Sigma_new| : " << L2_difference_norm / L2_Sigma_norm << std::endl;

  return L2_difference_norm / L2_Sigma_norm;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::compute_Sigma_new(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>>& G_r_w,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>>& GS_r_w) {
  Sigma_new = 0;

  FUNC_LIB::function<double, nu>& mu_HALF = ss_hybridization_solver_routines_type::get_mu_HALF();

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++)
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        if (ss_hybridization_solver_routines_type::is_interacting_band(b_ind))
          Sigma_new(b_ind, s_ind, b_ind, s_ind, 0, w_ind) =
              (GS_r_w(b_ind, s_ind, b_ind, s_ind, 0, w_ind) /
               G_r_w(b_ind, s_ind, b_ind, s_ind, 0, w_ind));

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++)
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        if (ss_hybridization_solver_routines_type::is_interacting_band(b_ind))
          Sigma_new(b_ind, s_ind, b_ind, s_ind, 0, w_ind) -= (mu_HALF(b_ind, s_ind));
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
int cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::find_w_cutoff() {
  ;
  return std::max(
      1.0, std::min(parameters.get_Sigma_tail_cutoff() * parameters.get_beta() / (2.0 * M_PI) - 0.5,
                    1.0 * (w::dmn_size() / 2)));
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<SS_CT_HYB, device_t, parameters_type, MOMS_type>::find_tail_of_Sigma(
    double& S0, double& S1, int b, int s, int k) {
  int w_cutoff = find_w_cutoff();
  S0 = 0.0;
  S1 = 0.0;

  S0 = real(Sigma_new(b, s, b, s, k, w::dmn_size() / 2 + w_cutoff - 1));
  S1 = imag(Sigma_new(b, s, b, s, k, w::dmn_size() / 2 + w_cutoff - 1)) *
       w::parameter_type::get_elements()[w::dmn_size() / 2 + w_cutoff - 1];
}

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_SOLVER_H
