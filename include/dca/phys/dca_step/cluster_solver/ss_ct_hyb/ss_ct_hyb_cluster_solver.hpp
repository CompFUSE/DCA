// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// Single-site Monte Carlo integrator based on a hybridization expansion.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_CLUSTER_SOLVER_HPP

#include <cassert>
#include <complex>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_walker.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_hybridization_solver_routines.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/util/plot.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
class SsCtHybClusterSolver : public cthyb::ss_hybridization_solver_routines<parameters_type, Data> {
public:
  using DataType = Data;
  typedef parameters_type this_parameters_type;

  typedef typename parameters_type::random_number_generator rng_type;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef cthyb::ss_hybridization_solver_routines<parameters_type, Data> ss_hybridization_solver_routines_type;

  typedef cthyb::SsCtHybWalker<dca::linalg::CPU, parameters_type, Data> walker_type;
  typedef cthyb::SsCtHybAccumulator<dca::linalg::CPU, parameters_type, Data> accumulator_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  using nu_nu_k_DCA_w = func::dmn_variadic<nu, nu, KClusterDmn, w>;

  const static int MC_TYPE = SS_CT_HYB;

public:
  SsCtHybClusterSolver(parameters_type& parameters_ref, Data& MOMS_ref);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

  template <typename Writer>
  void write(Writer& writer);

  // Computes and returns the local value of the Green's function G(k, \omega), i.e. without
  // averaging it across processes.
  // For testing purposes.
  // Precondition: The accumulator data has not been averaged, i.e. finalize has not been called.
  auto local_G_k_w() const;

  // Computes and returns the local value of the product of the Green's function G(r, \omega) and
  // the self-energy \Sigma(r, \omega), i.e. without averaging it across processes.
  // For testing purposes.
  // Precondition: The accumulator data has not been averaged, i.e. finalize has not been called.
  auto local_GS_r_w() const;

protected:
  void warm_up(walker_type& walker);

  void measure(walker_type& walker);

  void symmetrize_measurements();

  void compute_error_bars();

  // Sums/averages the quantities measured by the individual MPI ranks.
  void collect_measurements();

  double compute_S_k_w_from_G_k_w();

  void compute_Sigma_new(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>>& G_r_w,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>>& GS_r_w);

  int find_w_cutoff();

  void find_tail_of_Sigma(double& S0, double& S1, int b, int s, int k);

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

  func::function<double, nu> mu_DC;

private:
  bool averaged_;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
SsCtHybClusterSolver<device_t, parameters_type, Data>::SsCtHybClusterSolver(
    parameters_type& parameters_ref, Data& data_ref)
    : cthyb::ss_hybridization_solver_routines<parameters_type, Data>(parameters_ref, data_ref),

      parameters(parameters_ref),
      data_(data_ref),
      concurrency(parameters.get_concurrency()),

      thermalization_time(0),
      MC_integration_time(0),

      total_time(0),

      rng(concurrency.id(), concurrency.number_of_processors(), parameters.get_seed()),
      accumulator(parameters, data_),

      Sigma_old("Self-Energy-n-1-iteration"),
      Sigma_new("Self-Energy-n-0-iteration"),

      DCA_iteration(-1),
      averaged_(false) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t SS CT-HYB Integrator is born \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
template <typename Writer>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::write(Writer& writer) {
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

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::initialize(int dca_iteration) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t SS CT-HYB Integrator has started ( DCA-iteration : " << dca_iteration
              << ")\n\n";

  DCA_iteration = dca_iteration;

  Sigma_old = data_.Sigma_cluster;

  ss_hybridization_solver_routines_type::initialize_functions();

  accumulator.initialize(dca_iteration);

  averaged_ = false;

  if (concurrency.id() == concurrency.first()) {
    std::stringstream ss;
    ss.precision(6);
    ss << std::scientific;

    func::function<double, nu>& mu = this->get_mu();
    func::function<double, nu>& mu_half = this->get_mu_HALF();

    func::function<double, nu>& a0 = this->get_a0();
    func::function<double, nu>& a1 = this->get_a1();

    ss << "\n\n mu, mu_half, a0, a1\n\n";
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        ss << b_ind << "\t" << s_ind << "\t" << mu(b_ind, s_ind) << "\t" << mu_half(b_ind, s_ind)
           << "\t" << a0(b_ind, s_ind) << "\t" << a1(b_ind, s_ind) << "\n";

    ss << "\n\n";

    std::cout << ss.str();
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::integrate() {
  if (concurrency.id() == concurrency.first()) {
    std::cout << "QMC integration has started: " << dca::util::print_time() << std::endl;
  }

  walker_type walker(parameters, data_, rng);

  walker.initialize();

  warm_up(walker);

  measure(walker);

  if (concurrency.id() == concurrency.first()) {
    std::cout << "On-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters.get_measurements() << std::endl;

    walker.printSummary();
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
template <typename dca_info_struct_t>
double SsCtHybClusterSolver<device_t, parameters_type, Data>::finalize(
    dca_info_struct_t& dca_info_struct) {
  collect_measurements();
  symmetrize_measurements();

  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(accumulator.get_G_r_w(),
                                                                        data_.G_k_w);

  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(data_.G_k_w, data_.G_r_w);

  dca_info_struct.L2_Sigma_difference(DCA_iteration) = compute_S_k_w_from_G_k_w();

  // util::Plot::plotBandsLines(accumulator.get_G_r_w());
  // util::Plot::plotBandsLines(accumulator.get_GS_r_w());
  // util::Plot::plotBandsLines(data_.G_k_w);
  // util::Plot::plotBandsLines(data_.Sigma);

  if (concurrency.id() == concurrency.first()) {
    std::stringstream ss;
    ss.precision(6);
    ss << std::scientific;

    ss << "\n\n Sigma \n\n";
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++) {
      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
        double result = 0;
        for (int w_ind = 0; w_ind < 50; w_ind++)
          result += real(data_.Sigma(b_ind, s_ind, b_ind, s_ind, 0, w_ind)) / 50.;

        ss << b_ind << "\t" << s_ind << "\t" << result << "\n";
      }
    }
    ss << "\n\n";

    std::cout << ss.str();
  }

  for (int i = 0; i < b::dmn_size() * s::dmn_size(); i++)
    for (int j = 0; j < KClusterDmn::dmn_size(); j++)
      dca_info_struct.Sigma_zero_moment(i, j, DCA_iteration) = real(data_.Sigma(i, i, j, 0));

  double total = 1.e-6, integral = 0;
  for (int l = 0; l < accumulator.get_visited_expansion_order_k().size(); l++) {
    total += accumulator.get_visited_expansion_order_k()(l);
    integral += accumulator.get_visited_expansion_order_k()(l) * l;
  }

  dca_info_struct.average_expansion_order(DCA_iteration) = integral / total;

  dca_info_struct.sign(DCA_iteration) = accumulator.get_average_sign();

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t SS CT-HYB Integrator has finalized \n" << std::endl;

  return dca_info_struct.L2_Sigma_difference(DCA_iteration);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::warm_up(walker_type& walker) {
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
void SsCtHybClusterSolver<device_t, parameters_type, Data>::measure(walker_type& walker) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t measuring has started \n" << std::endl;

  const int n_meas = dca::parallel::util::getWorkload(parameters.get_measurements(), concurrency);

  for (int i = 0; i < n_meas; i++) {
    walker.do_sweep();

    accumulator.update_from(walker);

    accumulator.measure();

    walker.update_shell(i, n_meas);
  }

  // here we need to do a correction a la Andrey

  //  func::function<double, nu> correction_to_GS;

  // for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
  // for(int b_ind=0; b_ind<b::dmn_size(); b_ind++)
  // correction_to_GS(b_ind,s_ind) = 0//parameters.get_chemical_potential_DC()
  // + mu_DC(b_ind, s_ind)
  // + walker.mu_HALF(b_ind,s_ind)
  // + walker.a0(b_ind,s_ind);

  accumulator.finalize();

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t measuring has ended \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::compute_error_bars() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t computing the error-bars" << std::endl;

  const int nb_measurements = accumulator.get_number_of_measurements();
  double sign = accumulator.get_accumulated_sign() / double(nb_measurements);

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> G_r_w(
      "G_r_w_tmp");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> GS_r_w(
      "GS_r_w_tmp");

  for (int l = 0; l < G_r_w.size(); l++)
    G_r_w(l) = accumulator.get_G_r_w()(l) / double(nb_measurements * sign);

  for (int l = 0; l < GS_r_w.size(); l++)
    GS_r_w(l) = accumulator.get_GS_r_w()(l) / double(nb_measurements * sign);

  compute_Sigma_new(G_r_w, GS_r_w);

  concurrency.average_and_compute_stddev(Sigma_new, data_.get_Sigma_stdv());
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::collect_measurements() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t Collect measurements" << std::endl;

  const int nb_measurements = accumulator.get_number_of_measurements();

  // sum the sign
  double accumulated_sign = accumulator.get_accumulated_sign();
  concurrency.sum_and_average(accumulated_sign);

  // sum G_r_w
  concurrency.sum_and_average(accumulator.get_G_r_w());
  accumulator.get_G_r_w() /= accumulated_sign;

  // sum GS_r_w
  concurrency.sum_and_average(accumulator.get_GS_r_w());
  accumulator.get_GS_r_w() /= accumulated_sign;

  concurrency.sum(accumulator.get_visited_expansion_order_k());
  averaged_ = true;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::symmetrize_measurements() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t symmetrize measurements has started" << std::endl;

  symmetrize::execute(accumulator.get_G_r_w(), data_.H_symmetry);

  symmetrize::execute(accumulator.get_GS_r_w(), data_.H_symmetry);

  std::vector<int> flavors = parameters_type::model_type::get_flavors();
  assert(flavors.size() == b::dmn_size());

  func::function<std::complex<double>, b> f_val;
  func::function<std::complex<double>, b> f_tot;

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

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
double SsCtHybClusterSolver<device_t, parameters_type, Data>::compute_S_k_w_from_G_k_w() {
  double alpha = DCA_iteration > 0 ? parameters.get_self_energy_mixing_factor() : 1;

  double L2_difference_norm = 1.e-6;
  double L2_Sigma_norm = 1.e-6;

  int w_cutoff = find_w_cutoff();

  compute_Sigma_new(accumulator.get_G_r_w(), accumulator.get_GS_r_w());

  symmetrize::execute(Sigma_new, data_.H_symmetry);

  for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
    if (ss_hybridization_solver_routines_type::is_interacting_band(b_ind)) {
      for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++) {
        for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
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

            data_.Sigma(b_ind, s_ind, b_ind, s_ind, k_ind, w_ind) =
                alpha * (new_sigma) + (1 - alpha) * old_sigma;
          }
        }
      }
    }
  }

  symmetrize::execute(data_.Sigma, data_.H_symmetry);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t |Sigma_old-Sigma_new| : " << L2_difference_norm / L2_Sigma_norm << std::endl;

  return L2_difference_norm / L2_Sigma_norm;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::compute_Sigma_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>>& G_r_w,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>>& GS_r_w) {
  Sigma_new = 0;

  func::function<double, nu>& mu_HALF = ss_hybridization_solver_routines_type::get_mu_HALF();

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

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
int SsCtHybClusterSolver<device_t, parameters_type, Data>::find_w_cutoff() {
  return std::max(
      1.0,
      std::min(parameters.get_self_energy_tail_cutoff() * parameters.get_beta() / (2.0 * M_PI) - 0.5,
               1.0 * (w::dmn_size() / 2)));
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void SsCtHybClusterSolver<device_t, parameters_type, Data>::find_tail_of_Sigma(double& S0, double& S1,
                                                                               int b, int s, int k) {
  int w_cutoff = find_w_cutoff();
  S0 = 0.0;
  S1 = 0.0;

  S0 = real(Sigma_new(b, s, b, s, k, w::dmn_size() / 2 + w_cutoff - 1));
  S1 = imag(Sigma_new(b, s, b, s, k, w::dmn_size() / 2 + w_cutoff - 1)) *
       w::parameter_type::get_elements()[w::dmn_size() / 2 + w_cutoff - 1];
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
auto SsCtHybClusterSolver<device_t, parameters_type, Data>::local_G_k_w() const {
  if (averaged_)
    throw std::logic_error("The local data was already averaged.");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w;
  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(accumulator.get_G_r_w(),
                                                                        G_k_w);
  G_k_w /= accumulator.get_sign();

  return G_k_w;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
auto SsCtHybClusterSolver<device_t, parameters_type, Data>::local_GS_r_w() const {
  if (averaged_)
    throw std::logic_error("The local data was already averaged.");

  auto GS_r_w = accumulator.get_GS_r_w();
  GS_r_w /= accumulator.get_accumulated_sign();

  return GS_r_w;
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_CLUSTER_SOLVER_HPP
