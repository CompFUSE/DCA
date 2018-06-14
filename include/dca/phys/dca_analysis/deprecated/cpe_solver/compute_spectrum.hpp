// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the spectrum using a CPE analytic continuation method.
//
// TODO: If the CPE solver gets reactivated, need to
//       - write, read and broadcast G_k_t in DcaData,
//       - initialize (and write) frequency_domain_imag_axis in Parameters.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_COMPUTE_SPECTRUM_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_COMPUTE_SPECTRUM_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/statistics/gaussian_probability.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/dca_analysis/cpe_solver/continuous_pole_expansion.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_imag_axis.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/util/plot.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <class parameters_type, class basis_function_t>
class compute_spectrum {
public:
  using concurrency_type = typename parameters_type::concurrency_type;
  using random_number_generator = typename parameters_type::random_number_generator;

  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;
  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;
  using w_IMAG = func::dmn_0<domains::frequency_domain_imag_axis>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = func::dmn_variadic<nu, nu>;

  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  compute_spectrum(parameters_type& parameters);

  template <typename MOMS_imag_type, typename MOMS_real_type>
  void write(std::string file_name, MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real);
  template <typename Writer>
  void write(Writer& writer);

  template <typename MOMS_imag_type, typename MOMS_real_type>
  void execute(MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real);

private:
  template <typename MOMS_imag_type, typename MOMS_real_type>
  void execute_without_error_bars(MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real);

  template <typename MOMS_imag_type, typename MOMS_real_type>
  void execute_with_error_bars(MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real);

  template <typename MOMS_imag_type, typename MOMS_real_type>
  void test_A_w_versus_G_t(MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real);

  void generate_f_original(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w);

  void generate_new_sample(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w_mean,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w_stddev,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w_sample);

  template <typename k_dmn_t, typename w_imag_dmn_t, typename w_real_dmn_t>
  void perform_analytic_continuation(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_imag_dmn_t>>& S_k_w_imag,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_real_dmn_t>>& S_k_w_real);

  template <typename k_dmn_t, typename w_dmn_t>
  void compute_G_k_w_on_cluster(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G0_k_w,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& Sigma,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G_k_w);

  template <typename k_host_dmn_t, typename k_cluster_dmn_t, typename w_dmn_t>
  void compute_G_k_w_on_lattice(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_host_dmn_t>>& H_k,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_cluster_dmn_t, w_dmn_t>>& Sigma,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_cluster_dmn_t, w_dmn_t>>& G_k_w);

  template <typename k_dmn_t, typename w_dmn_t>
  void compute_A_w(
      func::function<double, w_dmn_t>& A_w,
      func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G_k_w);

  template <typename k_dmn_t>
  void compute_G_k_beta_over_2(func::function<double, func::dmn_variadic<b, s, k_dmn_t>>& G_k_beta_over_2,
                               func::function<double, func::dmn_variadic<nu, nu, k_dmn_t, t>>& G_k_t);

  template <typename k_dmn_t>
  void accumulate_integrated_A_k_div_cosh(
      func::function<double, func::dmn_variadic<b, s, k_dmn_t>>& integrated_A_k_div_cosh,
      func::function<double, func::dmn_variadic<b, s, k_dmn_t>>& integrated_A_k_div_cosh_stddev,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_REAL>>& G_k_w);

  template <typename scalar_type, typename dmn_type>
  void accumulate(func::function<scalar_type, dmn_type>& f,
                  func::function<scalar_type, dmn_type>& f_average,
                  func::function<scalar_type, dmn_type>& f_square);

  template <typename scalar_type, typename dmn_type>
  void accumulate(func::function<std::complex<scalar_type>, dmn_type>& f,
                  func::function<std::complex<scalar_type>, dmn_type>& f_average,
                  func::function<std::complex<scalar_type>, dmn_type>& f_square);

  template <typename k_dmn_t, typename w_dmn_t>
  void accumulate_A_w(
      func::function<double, w_dmn_t>& A_w, func::function<double, w_dmn_t>& A_w_stddev,
      func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w,
      func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w_stddev,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G_k_w);

  template <typename k_dmn_t, typename w_dmn_t>
  void accumulate_f_K_w(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_sample,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_mean,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_square);

  template <typename scalar_type, typename dmn_type>
  void compute_mean_and_stddev(int nb_sample, func::function<scalar_type, dmn_type>& f_average,
                               func::function<scalar_type, dmn_type>& f_square);

  template <typename scalar_type, typename dmn_type>
  void compute_mean_and_stddev(int nb_sample,
                               func::function<std::complex<scalar_type>, dmn_type>& f_average,
                               func::function<std::complex<scalar_type>, dmn_type>& f_square);

  template <typename w_dmn_t>
  void compute_mean_and_stddev(int nb_samples, func::function<double, w_dmn_t>& A_w,
                               func::function<double, w_dmn_t>& A_w_stddev);

  template <typename w_dmn_t>
  void compute_mean_and_stddev(
      int nb_samples, func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w,
      func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w_stddev);

  template <typename k_dmn_t, typename w_dmn_t>
  void compute_mean_and_stddev(
      int nb_samples,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_mean,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_stddev);

  template <typename MOMS_real_type>
  void print_check_sums(MOMS_real_type& MOMS_real);

  parameters_type& parameters;
  concurrency_type& concurrency;

  random_number_generator rng;

  continuous_pole_expansion<parameters_type, basis_function_t, k_DCA, w_REAL, WEIGHTED_GRADIENT_METHOD> cpe_obj;

  func::function<double, func::dmn_variadic<b, s, k_DCA>> G_k_beta_over_2;
  func::function<double, func::dmn_variadic<b, s, k_DCA>> integrated_A_k_div_cosh;
  func::function<double, func::dmn_variadic<b, s, k_DCA>> integrated_A_k_div_cosh_stddev;

  func::function<double, func::dmn_variadic<b, s, k_DCA>> G0_k_beta_over_2;
  func::function<double, func::dmn_variadic<b, s, k_DCA>> integrated_A0_k_div_cosh;
  func::function<double, func::dmn_variadic<b, s, k_DCA>> integrated_A0_k_div_cosh_stddev;

  func::function<double, func::dmn_variadic<nu, nu, k_DCA>> error_function;
  func::function<double, func::dmn_variadic<nu, nu, k_DCA>> error_function_stddev;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> S_approx;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> S_approx_stddev;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w_IMAG>> f_orig;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w_IMAG>> f_approx;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w_IMAG>> f_approx_stddev;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w_IMAG>> f_measured;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w_IMAG>> f_measured_stddev;
};

template <class parameters_type, class basis_function_t>
compute_spectrum<parameters_type, basis_function_t>::compute_spectrum(parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      rng(concurrency.id(), concurrency.number_of_processors(), parameters.get_seed()),

      cpe_obj(parameters, concurrency),

      G_k_beta_over_2("G_k_beta_over_2"),
      integrated_A_k_div_cosh("integrated_A_k_div_cosh"),
      integrated_A_k_div_cosh_stddev("integrated_A_k_div_cosh-stddev"),

      G0_k_beta_over_2("G0_k_beta_over_2"),
      integrated_A0_k_div_cosh("integrated_A0_k_div_cosh"),
      integrated_A0_k_div_cosh_stddev("integrated_A0_k_div_cosh-stddev"),

      error_function("L2-CPE-error"),
      error_function_stddev("L2-CPE-error-stddev"),

      S_approx("Sigma-approx"),
      S_approx_stddev("Sigma-approx-stddev"),

      f_orig("f-original"),

      f_approx("f-approx"),
      f_approx_stddev("f-approx-stddev"),

      f_measured("f-measured"),
      f_measured_stddev("f-measured-stddev") {}

template <class parameters_type, class basis_function_t>
template <typename MOMS_imag_type, typename MOMS_real_type>
void compute_spectrum<parameters_type, basis_function_t>::write(std::string file_name,
                                                                MOMS_imag_type& /*MOMS_imag*/,
                                                                MOMS_real_type& MOMS_real) {
  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  const std::string& output_format = parameters.get_output_format();

  if (output_format == "JSON") {
    dca::io::JSONWriter writer;
    writer.open_file(file_name);

    parameters.write(writer);
    MOMS_real.write(writer);
    cpe_obj.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else if (output_format == "HDF5") {
    dca::io::HDF5Writer writer;
    writer.open_file(file_name);

    parameters.write(writer);
    MOMS_real.write(writer);
    cpe_obj.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else
    throw std::logic_error(__FUNCTION__);
}

template <class parameters_type, class basis_function_t>
template <typename Writer>
void compute_spectrum<parameters_type, basis_function_t>::write(Writer& writer) {
  writer.open_group("CPE-functions");

  writer.execute(error_function);
  writer.execute(error_function_stddev);

  writer.execute(S_approx);
  writer.execute(S_approx_stddev);

  writer.execute(f_orig);

  writer.execute(f_approx);
  writer.execute(f_approx_stddev);

  writer.execute(f_measured);
  writer.execute(f_measured_stddev);

  writer.close_group();

  writer.open_group("A-versus-G-functions");

  writer.execute(G_k_beta_over_2);
  writer.execute(integrated_A_k_div_cosh);

  writer.execute(G0_k_beta_over_2);
  writer.execute(integrated_A0_k_div_cosh);

  writer.close_group();
}

template <class parameters_type, class basis_function_t>
template <typename MOMS_imag_type, typename MOMS_real_type>
void compute_spectrum<parameters_type, basis_function_t>::execute(MOMS_imag_type& MOMS_imag,
                                                                  MOMS_real_type& MOMS_real) {
  if (parameters.simulate_gaussian_noise())
    execute_with_error_bars(MOMS_imag, MOMS_real);
  else
    execute_without_error_bars(MOMS_imag, MOMS_real);
}

template <class parameters_type, class basis_function_t>
template <typename MOMS_imag_type, typename MOMS_real_type>
void compute_spectrum<parameters_type, basis_function_t>::execute_without_error_bars(
    MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t start analytic-continuation without error-bars (time = "
              << dca::util::print_time() << ")\n";

  if (parameters.compute_free_spectrum()) {
    MOMS_real.Sigma = 0;
    compute_G_k_w_on_lattice(MOMS_imag.H_HOST, MOMS_real.Sigma, MOMS_real.G0_k_w);

    compute_A_w(MOMS_real.A0_w, MOMS_real.A0_nu_w, MOMS_real.G0_k_w);

    util::Plot::plotLinesPoints(MOMS_real.A0_w);
    util::Plot::plotBandsLinesPoints(MOMS_real.A0_nu_w);

    compute_G_k_beta_over_2(G0_k_beta_over_2, MOMS_imag.G0_k_t);

    accumulate_integrated_A_k_div_cosh(integrated_A0_k_div_cosh, integrated_A0_k_div_cosh_stddev,
                                       MOMS_real.G0_k_w);
  }

  if (parameters.compute_cluster_spectrum() or parameters.compute_lattice_spectrum())
    perform_analytic_continuation(MOMS_imag.Sigma, MOMS_real.Sigma);

  if (parameters.compute_lattice_spectrum())
    compute_G_k_w_on_lattice(MOMS_imag.H_HOST, MOMS_real.Sigma, MOMS_real.G_k_w);

  if (parameters.compute_cluster_spectrum())
    compute_G_k_w_on_cluster(MOMS_real.G0_k_w, MOMS_real.Sigma, MOMS_real.G_k_w);

  if (parameters.compute_cluster_spectrum() or parameters.compute_lattice_spectrum()) {
    {
      compute_G_k_beta_over_2(G_k_beta_over_2, MOMS_imag.G_k_t);

      accumulate_integrated_A_k_div_cosh(integrated_A_k_div_cosh, integrated_A_k_div_cosh_stddev,
                                         MOMS_real.G_k_w);

      accumulate_A_w(MOMS_real.A_w, MOMS_real.A_w_stddev, MOMS_real.A_nu_w, MOMS_real.A_nu_w_stddev,
                     MOMS_real.G_k_w);

      //      accumulate_f_K_w(MOMS_real.Sigma, MOMS_real.Sigma, MOMS_real.Sigma_stddev);
      //      accumulate_f_K_w(MOMS_real.G_k_w, MOMS_real.G_k_w, MOMS_real.G_k_w_stddev);
    }

    {
      util::Plot::plotLinesPoints(MOMS_real.A_w);
      util::Plot::plotBandsLinesPoints(MOMS_real.A_nu_w);

      print_check_sums(MOMS_real);
    }

    {
      integrated_A_k_div_cosh_stddev = 0.;

      MOMS_real.A_w_stddev = 0.;
      MOMS_real.A_nu_w_stddev = 0.;

      MOMS_real.Sigma_stddev = 0.;
      MOMS_real.G_k_w_stddev = 0.;

      f_approx_stddev = 0.;
      f_measured_stddev = 0.;
    }

    if (concurrency.id() == concurrency.last()) {
      std::cout.precision(6);
      std::cout << std::scientific;

      std::cout
          << "\n\n\t beta/2*G(beta/2) versus 1./(2*T)\\int_{-inf}^{inf} dw A(w)/cosh(w/(2*T)) \n\n";

      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
          std::cout << "\t" << G_k_beta_over_2(0, 0, k_ind) << "\t"
                    << integrated_A_k_div_cosh(0, 0, k_ind) << "\t"
                    << integrated_A_k_div_cosh_stddev(0, 0, k_ind) << "\n";
    }
  }
}

template <class parameters_type, class basis_function_t>
template <typename MOMS_imag_type, typename MOMS_real_type>
void compute_spectrum<parameters_type, basis_function_t>::execute_with_error_bars(
    MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real) {
  if (parameters.compute_free_spectrum()) {
    MOMS_real.Sigma = 0;
    compute_G_k_w_on_lattice(MOMS_imag.H_HOST, MOMS_real.Sigma, MOMS_real.G0_k_w);

    compute_A_w(MOMS_real.A0_w, MOMS_real.A0_nu_w, MOMS_real.G0_k_w);

    util::Plot::plotLinesPoints(MOMS_real.A0_w);
    util::Plot::plotBandsLinesPoints(MOMS_real.A0_nu_w);
    // util::Plot::plotBandsLinesPoints(MOMS_real.G0_k_w);

    compute_G_k_beta_over_2(G0_k_beta_over_2, MOMS_imag.G0_k_t);
  }

  {
    generate_f_original(MOMS_imag.Sigma);

    compute_G_k_beta_over_2(G_k_beta_over_2, MOMS_imag.G_k_t);
  }

  int nb_samples = parameters.get_nr_of_CPE_samples();
  {
    double magnitude = parameters.get_simulated_CPE_stddev();

    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>> S_K_wi;
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w_REAL>> S_K_wr;

    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w_REAL>> G_K_wr;

    if (concurrency.id() == concurrency.first())
      std::cout << "\n\n";

    for (int l = 0; l < nb_samples; l++) {
      if (concurrency.id() == concurrency.first())
        std::cout << "\t start analytic-continuation on sample = " << l
                  << " (time = " << dca::util::print_time() << ")\n";

      {  // generate a new sample that is equal for each MPI-task!

        if (concurrency.id() == concurrency.first()) {
          // MOMS_imag.Sigma_stddev = magnitude;

          // uniform error
          for (int l = 0; l < MOMS_imag.Sigma_stddev.size(); l++) {
            std::complex<double> val = MOMS_imag.Sigma(l);

            MOMS_imag.Sigma_stddev(l).real(magnitude * (real(val) * real(val)) /
                                           (real(val) * real(val) + 1.e-6));
            MOMS_imag.Sigma_stddev(l).imag(magnitude * (imag(val) * imag(val)) /
                                           (imag(val) * imag(val) + 1.e-6));
          }

          generate_new_sample(MOMS_imag.Sigma, MOMS_imag.Sigma_stddev, S_K_wi);
        }
        else {
          S_K_wi = 0;
        }

        concurrency.sum(S_K_wi);
      }

      {  // do the analytic continuation ...
        if (parameters.compute_cluster_spectrum() or parameters.compute_lattice_spectrum())
          perform_analytic_continuation(S_K_wi, S_K_wr);

        if (parameters.compute_lattice_spectrum())
          compute_G_k_w_on_lattice(MOMS_imag.H_HOST, S_K_wr, G_K_wr);

        if (parameters.compute_cluster_spectrum())
          compute_G_k_w_on_cluster(MOMS_real.G0_k_w, S_K_wr, G_K_wr);
      }

      {  // accumulate the function and the squared-function ...

        accumulate_integrated_A_k_div_cosh(integrated_A_k_div_cosh, integrated_A_k_div_cosh_stddev,
                                           G_K_wr);

        accumulate_A_w(MOMS_real.A_w, MOMS_real.A_w_stddev, MOMS_real.A_nu_w,
                       MOMS_real.A_nu_w_stddev, G_K_wr);

        accumulate_f_K_w(S_K_wr, MOMS_real.Sigma, MOMS_real.Sigma_stddev);
        accumulate_f_K_w(G_K_wr, MOMS_real.G_k_w, MOMS_real.G_k_w_stddev);
      }

      if (concurrency.id() == concurrency.first())
        std::cout << "\n";
    }
  }

  {
    compute_mean_and_stddev(nb_samples, error_function, error_function_stddev);

    compute_mean_and_stddev(nb_samples, integrated_A_k_div_cosh, integrated_A_k_div_cosh_stddev);

    compute_mean_and_stddev(nb_samples, MOMS_real.A_w, MOMS_real.A_w_stddev);
    compute_mean_and_stddev(nb_samples, MOMS_real.A_nu_w, MOMS_real.A_nu_w_stddev);

    compute_mean_and_stddev(nb_samples, MOMS_real.Sigma, MOMS_real.Sigma_stddev);
    compute_mean_and_stddev(nb_samples, MOMS_real.G_k_w, MOMS_real.G_k_w_stddev);

    compute_mean_and_stddev(nb_samples, f_approx, f_approx_stddev);
    compute_mean_and_stddev(nb_samples, f_measured, f_measured_stddev);
  }

  {
    util::Plot::plotErrorBars(MOMS_real.A_w, MOMS_real.A_w_stddev);

    print_check_sums(MOMS_real);

    if (concurrency.id() == concurrency.last()) {
      std::cout.precision(6);
      std::cout << std::scientific;

      std::cout
          << "\n\n\t beta/2*G(beta/2) versus 1./(2*T)\\int_{-inf}^{inf} dw A(w)/cosh(w/(2*T)) \n\n";

      for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
        for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
          std::cout << "\t" << G_k_beta_over_2(0, 0, k_ind) << "\t"
                    << integrated_A_k_div_cosh(0, 0, k_ind) << "\t"
                    << integrated_A_k_div_cosh_stddev(0, 0, k_ind) << "\n";
    }

    // test_A_w_versus_G_t(MOMS_imag, MOMS_real);
  }
}

template <class parameters_type, class basis_function_t>
void compute_spectrum<parameters_type, basis_function_t>::generate_f_original(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w) {
  for (int w_0 = 0; w_0 < w::dmn_size(); w_0++) {
    for (int w_1 = 0; w_1 < w_IMAG::dmn_size(); w_1++) {
      if (std::abs(w::get_elements()[w_0] - w_IMAG::get_elements()[w_1]) < 1.e-6) {
        for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
          for (int nu_j = 0; nu_j < nu::dmn_size(); nu_j++)
            for (int nu_i = 0; nu_i < nu::dmn_size(); nu_i++)
              f_orig(nu_i, nu_j, k_ind, w_1) = S_K_w(nu_i, nu_j, k_ind, w_0);
      }
    }
  }

  util::Plot::plotLinesPoints(f_orig);
}

template <class parameters_type, class basis_function_t>
void compute_spectrum<parameters_type, basis_function_t>::generate_new_sample(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w_mean,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w_stddev,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w_sample) {
  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
      for (int b_j = 0; b_j < b::dmn_size(); b_j++) {
        for (int b_i = 0; b_i < b::dmn_size(); b_i++) {
          if ((b_i == b_j) and parameters.is_an_interacting_band(b_i)) {
            double error_re = rng();
            double error_im = rng();

            std::complex<double> mean = S_K_w_mean(b_i, 0, b_j, 0, k_ind, w_ind);
            std::complex<double> stddev = S_K_w_stddev(b_i, 0, b_j, 0, k_ind, w_ind);

            std::complex<double> error(
                math::statistics::gauss::argTailProbability(error_re, 0, real(stddev)),
                math::statistics::gauss::argTailProbability(error_im, 0, imag(stddev)));

            assert(error_re == error_re);
            assert(error_im == error_im);

            S_K_w_sample(b_i, 0, b_j, 0, k_ind, w_ind) = mean + error;
            S_K_w_sample(b_i, 1, b_j, 1, k_ind, w_ind) = mean + error;
          }
          else {
            S_K_w_sample(b_i, 0, b_j, 0, k_ind, w_ind) = 0.;
            S_K_w_sample(b_i, 1, b_j, 1, k_ind, w_ind) = 0.;
          }
        }
      }
    }
  }

  { symmetrize::execute(S_K_w_sample); }
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t, typename w_imag_dmn_t, typename w_real_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::perform_analytic_continuation(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_imag_dmn_t>>& S_k_w_imag,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_real_dmn_t>>& S_k_w_real) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\t\t start CPE (time = " << dca::util::print_time() << ") --> ";

  // cpe_obj.execute_st(S_k_w_imag, S_k_w_real);
  {
    cpe_obj.execute_mt(S_k_w_imag, S_k_w_real);

    accumulate(cpe_obj.get_error_function(), error_function, error_function_stddev);

    accumulate_f_K_w(cpe_obj.get_S_approx(), S_approx, S_approx_stddev);

    accumulate_f_K_w(cpe_obj.get_f_approx(), f_approx, f_approx_stddev);
    accumulate_f_K_w(cpe_obj.get_f_measured(), f_measured, f_measured_stddev);
  }

  if (concurrency.id() == concurrency.first())
    std::cout << " (time = " << dca::util::print_time() << ")\n";
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t, typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::compute_G_k_w_on_cluster(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G0_k_w,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& Sigma,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G_k_w) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\t\t start AC on G_K_w (time = " << dca::util::print_time() << ") --> ";

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix("G-matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> S_matrix("S-matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix("G0-matrix", nu::dmn_size());

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> S_times_G0_matrix("SxG0_matrix",
                                                                                nu::dmn_size());
  // Allocate the work space for inverse only once.
  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          S_matrix(i, j) = Sigma(i, j, k_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = G0_k_w(i, j, k_ind, w_ind);

      dca::linalg::matrixop::gemm(S_matrix, G0_matrix, S_times_G0_matrix);

      for (int i = 0; i < nu::dmn_size(); i++)
        S_times_G0_matrix(i, i) = 1. - S_times_G0_matrix(i, i);

      dca::linalg::matrixop::inverse(S_times_G0_matrix, ipiv, work);

      dca::linalg::matrixop::gemm(S_times_G0_matrix, G0_matrix, G_matrix);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_k_w(i, j, k_ind, w_ind) = G_matrix(i, j);
    }
  }

  if (concurrency.id() == concurrency.first())
    std::cout << " (time = " << dca::util::print_time() << ")\n";
}

template <class parameters_type, class basis_function_t>
template <typename k_host_dmn_t, typename k_cluster_dmn_t, typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::compute_G_k_w_on_lattice(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_host_dmn_t>>& H_k,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_cluster_dmn_t, w_dmn_t>>& Sigma,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_cluster_dmn_t, w_dmn_t>>& G_k_w) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\t\t start TIM (time = " << dca::util::print_time() << ") --> ";

  clustermapping::coarsegraining_sp<parameters_type, k_DCA> coarsegraining_sp_obj(parameters);

  coarsegraining_sp_obj.compute_G_K_w_with_TIM(H_k, Sigma, G_k_w);

  if (concurrency.id() == concurrency.first())
    std::cout << " (time = " << dca::util::print_time() << ")\n";
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t, typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::compute_A_w(
    func::function<double, w_dmn_t>& A_w,
    func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G_k_w) {
  A_w = 0;
  A_nu_w = 0;

  for (int b = 0; b < b::dmn_size(); b++) {
    for (int s = 0; s < s::dmn_size(); s++) {
      for (int w_ind = 0; w_ind < w_REAL::dmn_size(); w_ind++) {
        for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
          double term =
              -1. / M_PI * imag(G_k_w(b, s, b, s, k_ind, w_ind)) / double(k_DCA::dmn_size());

          A_w(w_ind) += term / s::dmn_size();
          A_nu_w(b, s, w_ind) += term;
        }
      }
    }
  }
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::compute_G_k_beta_over_2(
    func::function<double, func::dmn_variadic<b, s, k_dmn_t>>& G_k_beta_over_2,
    func::function<double, func::dmn_variadic<nu, nu, k_dmn_t, t>>& G_k_t) {
  int t_ind = 3. * t::dmn_size() / 4.;
  assert(std::abs(t::get_elements()[t_ind] - parameters.get_beta() / 2.) < 1.e-6);

  for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
      for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++)
        G_k_beta_over_2(b_ind, s_ind, k_ind) =
            -parameters.get_beta() * G_k_t(b_ind, s_ind, b_ind, s_ind, k_ind, t_ind);
}

template <class parameters_type, class basis_function_t>
template <typename scalar_type, typename dmn_type>
void compute_spectrum<parameters_type, basis_function_t>::accumulate(
    func::function<scalar_type, dmn_type>& f, func::function<scalar_type, dmn_type>& f_average,
    func::function<scalar_type, dmn_type>& f_square) {
  for (int i = 0; i < f.size(); i++) {
    f_average(i) += f(i);
    f_square(i) += f(i) * f(i);
  }
}

template <class parameters_type, class basis_function_t>
template <typename scalar_type, typename dmn_type>
void compute_spectrum<parameters_type, basis_function_t>::accumulate(
    func::function<std::complex<scalar_type>, dmn_type>& f,
    func::function<std::complex<scalar_type>, dmn_type>& f_average,
    func::function<std::complex<scalar_type>, dmn_type>& f_square) {
  for (int i = 0; i < f.size(); i++) {
    f_average(i) += f(i);

    real(f_square(i)) += real(f(i)) * real(f(i));
    imag(f_square(i)) += imag(f(i)) * imag(f(i));
  }
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::accumulate_integrated_A_k_div_cosh(
    func::function<double, func::dmn_variadic<b, s, k_dmn_t>>& int_A_k_div_cosh,
    func::function<double, func::dmn_variadic<b, s, k_dmn_t>>& int_A_k_div_cosh_stddev,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_REAL>>& G_k_w) {
  double twoT = 2. / parameters.get_beta();

  double result = 0;
  for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
    for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++) {
      for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
        result = 0;
        for (int w_ind = 0; w_ind < w_REAL::dmn_size(); w_ind++)
          result += -1. / (twoT * M_PI) * imag(G_k_w(b_ind, s_ind, b_ind, s_ind, k_ind, w_ind)) /
                    (std::cosh(w_REAL::get_elements()[w_ind] / twoT)) *
                    basis_function_t::volume(w_ind);

        int_A_k_div_cosh(b_ind, s_ind, k_ind) += result;
        int_A_k_div_cosh_stddev(b_ind, s_ind, k_ind) += result * result;
      }
    }
  }
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t, typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::accumulate_A_w(
    func::function<double, w_dmn_t>& A_w, func::function<double, w_dmn_t>& A_w_stddev,
    func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w,
    func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w_stddev,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& G_k_w) {
  for (int w_ind = 0; w_ind < w_REAL::dmn_size(); w_ind++) {
    double A_w_result = 0;

    for (int b = 0; b < b::dmn_size(); b++) {
      for (int s = 0; s < s::dmn_size(); s++) {
        double A_nu_w_result = 0;

        for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
          A_w_result +=
              -1. / M_PI * imag(G_k_w(b, s, b, s, k_ind, w_ind)) / double(k_DCA::dmn_size());
          A_nu_w_result +=
              -1. / M_PI * imag(G_k_w(b, s, b, s, k_ind, w_ind)) / double(k_DCA::dmn_size());
        }

        A_nu_w(b, s, w_ind) += A_nu_w_result;
        A_nu_w_stddev(b, s, w_ind) += A_nu_w_result * A_nu_w_result;
      }
    }

    A_w_result /= s::dmn_size();

    A_w(w_ind) += A_w_result;
    A_w_stddev(w_ind) += A_w_result * A_w_result;
  }
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t, typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::accumulate_f_K_w(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_sample,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_mean,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_square) {
  for (int i = 0; i < f_K_w_sample.size(); i++) {
    f_K_w_mean(i) += f_K_w_sample(i);

    f_K_w_square(i).real(real(f_K_w_square(i)) + real(f_K_w_sample(i)) * real(f_K_w_sample(i)));
    f_K_w_square(i).imag(imag(f_K_w_square(i)) + imag(f_K_w_sample(i)) * imag(f_K_w_sample(i)));
  }
}

template <class parameters_type, class basis_function_t>
template <typename scalar_type, typename dmn_type>
void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(
    int nb_samples, func::function<scalar_type, dmn_type>& f_average,
    func::function<scalar_type, dmn_type>& f_square) {
  f_average /= scalar_type(nb_samples);
  f_square /= scalar_type(nb_samples);

  for (int i = 0; i < f_average.size(); i++)
    f_square(i) = sqrt(abs(f_square(i) - f_average(i) * f_average(i)));
}

template <class parameters_type, class basis_function_t>
template <typename scalar_type, typename dmn_type>
void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(
    int nb_samples, func::function<std::complex<scalar_type>, dmn_type>& f_average,
    func::function<std::complex<scalar_type>, dmn_type>& f_square) {
  f_average /= scalar_type(nb_samples);
  f_square /= scalar_type(nb_samples);

  for (int i = 0; i < f_average.size(); i++) {
    real(f_square(i)) = sqrt(abs(real(f_square(i)) - real(f_average(i)) * real(f_average(i))));
    imag(f_square(i)) = sqrt(abs(imag(f_square(i)) - imag(f_average(i)) * imag(f_average(i))));
  }
}

template <class parameters_type, class basis_function_t>
template <typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(
    int nb_samples, func::function<double, w_dmn_t>& A_w,
    func::function<double, w_dmn_t>& A_w_stddev) {
  A_w /= double(nb_samples);
  A_w_stddev /= double(nb_samples);

  for (int i = 0; i < A_w.size(); i++)
    A_w_stddev(i) = sqrt(std::abs(A_w_stddev(i) - A_w(i) * A_w(i)));
}

template <class parameters_type, class basis_function_t>
template <typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(
    int nb_samples, func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w,
    func::function<double, func::dmn_variadic<b, s, w_dmn_t>>& A_nu_w_stddev) {
  A_nu_w /= double(nb_samples);
  A_nu_w_stddev /= double(nb_samples);

  for (int i = 0; i < A_nu_w.size(); i++)
    A_nu_w_stddev(i) = sqrt(std::abs(A_nu_w_stddev(i) - A_nu_w(i) * A_nu_w(i)));
}

template <class parameters_type, class basis_function_t>
template <typename k_dmn_t, typename w_dmn_t>
void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(
    int nb_samples,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_K_w_stddev) {
  f_K_w /= double(nb_samples);
  f_K_w_stddev /= double(nb_samples);

  for (int i = 0; i < f_K_w.size(); i++) {
    f_K_w_stddev(i).real(std::sqrt(std::abs(real(f_K_w_stddev(i)) - real(f_K_w(i)) * real(f_K_w(i)))));
    f_K_w_stddev(i).imag(std::sqrt(std::abs(imag(f_K_w_stddev(i)) - imag(f_K_w(i)) * imag(f_K_w(i)))));
  }
}

template <class parameters_type, class basis_function_t>
template <typename MOMS_real_type>
void compute_spectrum<parameters_type, basis_function_t>::print_check_sums(MOMS_real_type& MOMS_real) {
  if (concurrency.id() == concurrency.first()) {
    double result = 0;

    {
      std::cout << "\n\n\t integrated G0 and G : \n\n";
      for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
        math::util::print(k_DCA::get_elements()[k_ind]);

        result = 0;
        for (int w_ind = 0; w_ind < w_REAL::dmn_size() - 1; w_ind++) {
          double f_w = imag(MOMS_real.G0_k_w(0, 0, 0, w_ind + 1) + MOMS_real.G0_k_w(0, 0, 0, w_ind));
          double delta_w = (w_REAL::get_elements()[w_ind + 1] - w_REAL::get_elements()[w_ind]);

          result += -1. / M_PI * delta_w * f_w / 2;
        }
        std::cout << result << "\t";

        result = 0;
        for (int w_ind = 0; w_ind < w_REAL::dmn_size() - 1; w_ind++) {
          double f_w = imag(MOMS_real.G_k_w(0, 0, 0, w_ind + 1) + MOMS_real.G_k_w(0, 0, 0, w_ind));
          double delta_w = (w_REAL::get_elements()[w_ind + 1] - w_REAL::get_elements()[w_ind]);

          result += -1. / M_PI * delta_w * f_w / 2;
        }
        std::cout << result << "\n";
      }
    }

    std::cout << "\n\n";

    {
      std::cout << "integrated A0 and A: \n\n\t";

      result = 0;
      for (int w_ind = 0; w_ind < w_REAL::dmn_size() - 1; w_ind++)
        result += (w_REAL::get_elements()[w_ind + 1] - w_REAL::get_elements()[w_ind]) *
                  (MOMS_real.A0_w(w_ind) + MOMS_real.A0_w(w_ind + 1)) / 2;

      std::cout << result << "\t";

      result = 0;
      for (int w_ind = 0; w_ind < w_REAL::dmn_size() - 1; w_ind++)
        result += (w_REAL::get_elements()[w_ind + 1] - w_REAL::get_elements()[w_ind]) *
                  (MOMS_real.A_w(w_ind) + MOMS_real.A_w(w_ind + 1)) / 2;

      std::cout << result << "\n";
    }
  }
}

template <class parameters_type, class basis_function_t>
template <typename MOMS_imag_type, typename MOMS_real_type>
void compute_spectrum<parameters_type, basis_function_t>::test_A_w_versus_G_t(
    MOMS_imag_type& MOMS_imag, MOMS_real_type& MOMS_real) {
  double twoT = 2. / parameters.get_beta();

  double result = 0;
  for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
    result = 0;
    for (int w_ind = 0; w_ind < w_REAL::dmn_size(); w_ind++)
      result += -1. / M_PI * imag(MOMS_real.G_k_w(0, 0, k_ind, w_ind)) /
                (std::cosh(w_REAL::get_elements()[w_ind] / twoT)) *
                basis_function_t::volume(w_ind) / twoT;

    integrated_A_k_div_cosh(k_ind) = result;

    result = 0;
    for (int w_ind = 0; w_ind < w_REAL::dmn_size(); w_ind++)
      result += -1. / M_PI * imag(MOMS_real.G0_k_w(0, 0, k_ind, w_ind)) /
                (std::cosh(w_REAL::get_elements()[w_ind] / twoT)) *
                basis_function_t::volume(w_ind) / twoT;

    integrated_A0_k_div_cosh(k_ind) = result;
  }

  for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++) {
    G_k_beta_over_2(k_ind) =
        -parameters.get_beta() * MOMS_imag.G_k_t(0, 0, k_ind, 3 * t::dmn_size() / 4);
    G0_k_beta_over_2(k_ind) =
        -parameters.get_beta() * MOMS_imag.G0_k_t(0, 0, k_ind, 3 * t::dmn_size() / 4);
  }

  if (concurrency.id() == concurrency.last()) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::cout << "\n\n\t G0(beta/2) versus \\int_{-inf}^{inf} dw A0(w)/cosh(w) \n\n";

    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      std::cout << "\t" << G0_k_beta_over_2(k_ind) << "\t" << integrated_A0_k_div_cosh(k_ind) << "\n";

    std::cout << "\n\n\t G(beta/2) versus \\int_{-inf}^{inf} dw A(w)/cosh(w) \n\n";

    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      std::cout << "\t" << G_k_beta_over_2(k_ind) << "\t" << integrated_A_k_div_cosh(k_ind) << "\n";

    std::cout << "\n";
  }
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_COMPUTE_SPECTRUM_HPP
