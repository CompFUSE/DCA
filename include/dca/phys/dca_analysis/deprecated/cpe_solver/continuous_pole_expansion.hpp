// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implemenents the continuous pole expansion method for solving the analytic continution
// problem. It is templated on the minimization algorithm.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_CONTINUOUS_POLE_EXPANSION_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_CONTINUOUS_POLE_EXPANSION_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/parallel/util/get_bounds.hpp"
#include "dca/parallel/util/threading_data.hpp"
#include "dca/phys/dca_analysis/cpe_solver/cpe_data.hpp"
#include "dca/phys/dca_analysis/cpe_solver/minimization_method.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_imag_axis.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

// Empty class template.
template <class parameters_type, class basis_function_t, class k_dmn_t, class w_dmn_t,
          MinimizationMethod minimization_method_t>
class continuous_pole_expansion {};

// Specialization for a CPE analytic continuation using a weighted gradient method as the
// minimization algorithm.
template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
class continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                                WEIGHTED_GRADIENT_METHOD> {
public:
  enum gradient_minimization_state {
    BELOW_MAX_ERROR,
    FOUND_THE_MINIMUM,
    FOUND_A_MINIMUM,
    END_OF_PATH
  };
  using gradient_minimization_state_t = gradient_minimization_state;

  using scalartype = double;

  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;
  using Threading = typename parameters_type::ThreadingType;

  using w = func::dmn_0<domains::frequency_domain>;
  using w_IMAG = func::dmn_0<domains::frequency_domain_imag_axis>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  using nu_nu_k_dmn = func::dmn_variadic<nu, nu, k_dmn_t>;
  using nu_nu_k_dmn_w = func::dmn_variadic<nu, nu, k_dmn_t, w>;
  using nu_nu_k_dmn_w_IMAG = func::dmn_variadic<nu, nu, k_dmn_t, w_IMAG>;
  using nu_nu_k_DCA_w_IMAG = func::dmn_variadic<nu, nu, k_DCA, w_IMAG>;

  using alpha_dmn_t = func::dmn_0<basis_function_t>;
  using nu_nu_k_dmn_alpha_dmn = func::dmn_variadic<nu, nu, k_dmn_t, alpha_dmn_t>;

  using CPE_data_type = CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>;

  continuous_pole_expansion(parameters_type& parameters, concurrency_type& concurrency,
                            bool fixed_zero_moment = false, double zero_moment = 0,
                            bool fixed_first_moment = false, double first_moment = 1);

  func::function<scalartype, nu_nu_k_dmn>& get_error_function() {
    return error_function;
  }

  func::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG>& get_f_approx() {
    return f_approx;
  }
  func::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG>& get_f_measured() {
    return f_measured;
  }

  func::function<std::complex<scalartype>, nu_nu_k_dmn_w>& get_S_approx() {
    return S_approx;
  }

  void execute_st(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& f_source,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_target);

  void execute_mt(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& f_source,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_target);

  static void* threaded_analytical_continuation(void* void_ptr);

  void write(std::string filename);

private:
  void initialize();

  static void read_function_values(int nu_ind, int k_ind, CPE_data_type& CPE_data_obj);

  static void perform_continuous_pole_expansion_threaded_1(CPE_data_type& CPE_data_obj);
  static void perform_continuous_pole_expansion_threaded_2(CPE_data_type& CPE_data_obj);

  static void write_function_values(int nu_ind, int k_ind, CPE_data_type& CPE_data_obj);

  static void compute_gradient_alphas(CPE_data_type& CPE_data_obj);

  static void find_new_Sigma_0(CPE_data_type& CPE_data_obj);

  static void compute_gradient_Sigma_0(CPE_data_type& CPE_data_obj);

  static gradient_minimization_state find_new_alpha_1(CPE_data_type& CPE_data_obj);
  static int find_new_alpha_2(CPE_data_type& CPE_data_obj);

  static scalartype find_minimum_1(int index, CPE_data_type& CPE_data_obj);
  static scalartype find_minimum_2(int index, CPE_data_type& CPE_data_obj);

  static void compute_new_alpha(double lambda, CPE_data_type& CPE_data_obj);

  static double evaluate_Lambda_norm(CPE_data_type& CPE_data_obj);

  static void project_from_real_axis_to_imaginary_axis(
      dca::linalg::Vector<scalartype, dca::linalg::CPU>& real_values,
      dca::linalg::Vector<scalartype, dca::linalg::CPU>& imag_values,
      dca::linalg::Matrix<scalartype, dca::linalg::CPU>& matrix);

  static void project_from_imaginary_axis_to_real_axis(
      dca::linalg::Vector<scalartype, dca::linalg::CPU>& imag_values,
      dca::linalg::Vector<scalartype, dca::linalg::CPU>& real_values,
      dca::linalg::Matrix<scalartype, dca::linalg::CPU>& matrix);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

public:
  func::function<scalartype, nu_nu_k_dmn> Sigma_0_moment;
  func::function<scalartype, nu_nu_k_dmn_alpha_dmn> alpha_function;
  func::function<scalartype, nu_nu_k_dmn> error_function;

  func::function<std::complex<scalartype>, nu_nu_k_dmn_w> S_approx;

  func::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG> f_approx;
  func::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG> f_measured;

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> A_matrix;

  dca::linalg::Matrix<scalartype, dca::linalg::CPU> A_matrix_re;
  dca::linalg::Matrix<scalartype, dca::linalg::CPU> A_matrix_im;
};

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::
    continuous_pole_expansion(parameters_type& parameters_ref, concurrency_type& concurrency_ref,
                              bool /*fixed_zeroth_moment*/, double /*zeroth_moment_val*/,
                              bool /*fixed_first_moment*/, double /*first_moment_val*/)
    : parameters(parameters_ref),
      concurrency(concurrency_ref),

      alpha_function("alpha"),
      error_function("L2-CPE-error"),

      S_approx("Sigma-approx"),

      f_approx("f-approx"),
      f_measured("f-measured"),

      A_matrix("A_matrix"),

      A_matrix_re("A_matrix_re"),
      A_matrix_im("A_matrix_im") {
  basis_function_t::initialize(parameters);
  initialize();
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::write(std::string file_name) {
  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  const std::string& output_format = parameters.get_output_format();

  if (output_format == "JSON") {
    dca::io::JSONWriter writer;
    writer.open_file(file_name);

    parameters.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else if (output_format == "HDF5") {
    dca::io::HDF5Writer writer;
    writer.open_file(file_name);

    parameters.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else
    throw std::logic_error(__FUNCTION__);
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::initialize() {
  {
    A_matrix.resizeNoCopy(std::pair<int, int>(w_IMAG::dmn_size(), alpha_dmn_t::dmn_size()));

    A_matrix_re.resizeNoCopy(std::pair<int, int>(w_IMAG::dmn_size(), alpha_dmn_t::dmn_size()));
    A_matrix_im.resizeNoCopy(std::pair<int, int>(w_IMAG::dmn_size(), alpha_dmn_t::dmn_size()));

    for (int wn_ind = 0; wn_ind < w_IMAG::dmn_size(); wn_ind++) {
      std::complex<double> z(0., w_IMAG::get_elements()[wn_ind]);

      for (int alpha_ind = 0; alpha_ind < alpha_dmn_t::dmn_size(); alpha_ind++) {
        A_matrix(wn_ind, alpha_ind) = basis_function_t::phi(alpha_ind, z);

        A_matrix_re(wn_ind, alpha_ind) = real(A_matrix(wn_ind, alpha_ind));
        A_matrix_im(wn_ind, alpha_ind) = imag(A_matrix(wn_ind, alpha_ind));
      }
    }
  }
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::execute_st(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& f_source,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_target) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  // initialize f-measured
  for (int w_ind = 0; w_ind < w_IMAG::dmn_size(); w_ind++)
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++)
      for (int nu_j = 0; nu_j < nu::dmn_size(); nu_j++)
        for (int nu_i = 0; nu_i < nu::dmn_size(); nu_i++)
          f_measured(nu_i, nu_j, k_ind, w_ind) =
              f_source(nu_i, nu_j, k_ind, w::dmn_size() / 2 + w_ind);

  CPE_data_type CPE_data_obj;

  CPE_data_obj.initialize(parameters, f_target, *this);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
    for (int nu_ind = 0; nu_ind < nu::dmn_size(); nu_ind++) {
      read_function_values(nu_ind, k_ind, CPE_data_obj);

      // find the alpha
      perform_continuous_pole_expansion_threaded(CPE_data_obj);
      write_function_values(nu_ind, k_ind, CPE_data_obj);
    }
  }

  symmetrize::execute(f_target);
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::execute_mt(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& f_source,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_target) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  {
    S_approx = 0.;

    f_approx = 0;
    f_target = 0;

    f_measured = 0;

    Sigma_0_moment = 0;
    alpha_function = 0;
    error_function = 0;
  }

  symmetrize::execute(f_source);

  {  // initialize f-measured
    for (int w_ind = 0; w_ind < w_IMAG::dmn_size(); w_ind++)
      for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++)
        for (int nu_j = 0; nu_j < nu::dmn_size(); nu_j++)
          for (int nu_i = 0; nu_i < nu::dmn_size(); nu_i++)
            f_measured(nu_i, nu_j, k_ind, w_ind) =
                f_source(nu_i, nu_j, k_ind, w::dmn_size() / 2 + w_ind);
  }

  {
    func::dmn_variadic<b, s, k_dmn_t> b_s_k_dmn;
    std::pair<int, int> bounds = concurrency.get_bounds(b_s_k_dmn);

    int nr_tasks = bounds.second - bounds.first;

    if (nr_tasks > 0) {
      int nr_threads = std::min(8, nr_tasks);

      std::vector<CPE_data_type> CPE_data_vector(nr_threads);

      for (int l = 0; l < nr_threads; l++)
        CPE_data_vector[l].initialize(l, bounds, parameters, f_target, *this);

      Threading parallelization_obj;

      parallelization_obj.execute(nr_threads, threaded_analytical_continuation,
                                  (void*)&CPE_data_vector);
    }

    // sum
    concurrency.sum(Sigma_0_moment);
    concurrency.sum(alpha_function);
    concurrency.sum(error_function);

    concurrency.sum(S_approx);
    concurrency.sum(f_approx);
    concurrency.sum(f_target);
  }

  symmetrize::execute(f_target);
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void* continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                                WEIGHTED_GRADIENT_METHOD>::threaded_analytical_continuation(void* void_ptr) {
  dca::parallel::ThreadingData* data_ptr = static_cast<dca::parallel::ThreadingData*>(void_ptr);
  std::vector<CPE_data_type>* CPE_data_vec_ptr =
      static_cast<std::vector<CPE_data_type>*>(data_ptr->arg);

  int id = data_ptr->id;
  int nr_threads = data_ptr->num_threads;

  std::vector<CPE_data_type>& CPE_data_vec = *(CPE_data_vec_ptr);

  std::pair<int, int> MPI_bounds = CPE_data_vec[id].bounds;

  func::dmn_variadic<b, s> nu_dmn;
  func::dmn_variadic<b, s, k_dmn_t> b_s_k_dmn;

  std::pair<int, int> bounds = dca::parallel::util::getBounds(id, nr_threads, MPI_bounds);

  int coor[3];
  for (int l = bounds.first; l < bounds.second; l++) {
    b_s_k_dmn.linind_2_subind(l, coor);

    int nu_ind = nu_dmn(coor[0], coor[1]);
    int k_ind = coor[2];

    read_function_values(nu_ind, k_ind, CPE_data_vec[id]);
    if (CPE_data_vec[id].is_interacting_band[nu_ind]) {  // find the alpha
      perform_continuous_pole_expansion_threaded_1(CPE_data_vec[id]);
    }
    else {
      for (int n_ind = 0; n_ind < alpha_dmn_t::dmn_size(); n_ind++)
        CPE_data_vec[id].alpha_vec_d[n_ind] = 0.;
    }

    write_function_values(nu_ind, k_ind, CPE_data_vec[id]);
  }
  return 0;
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::read_function_values(int nu_ind, int k_ind,
                                                                               CPE_data_type& CPE_data_obj) {
  func::function<std::complex<scalartype>, func::dmn_variadic<nu, nu, k_dmn_t, w_IMAG>>& f_measured_func =
      *(CPE_data_obj.f_measured_ptr);

  {  // read in the values on the imaginary axis.
    for (int w_ind = 0; w_ind < w_IMAG::dmn_size(); w_ind++) {
      std::complex<double> value = f_measured_func(nu_ind, nu_ind, k_ind, w_ind);

      CPE_data_obj.F_wn_vec[w_ind] = value;

      CPE_data_obj.F_wn_re_vec[w_ind] = real(value);
      CPE_data_obj.F_wn_im_vec[w_ind] = imag(value);
    }
  }
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::write_function_values(int nu_ind, int k_ind,
                                                                                CPE_data_type&
                                                                                    CPE_data_obj) {
  {
    func::function<scalartype, func::dmn_variadic<nu, nu, k_dmn_t>>& Sigma_0_func =
        *(CPE_data_obj.Sigma_0_moment_ptr);
    func::function<scalartype, func::dmn_variadic<nu, nu, k_dmn_t, alpha_dmn_t>>& alpha_func =
        *(CPE_data_obj.alpha_function_ptr);
    func::function<scalartype, func::dmn_variadic<nu, nu, k_dmn_t>>& error_func =
        *(CPE_data_obj.error_function_ptr);

    Sigma_0_func(nu_ind, nu_ind, k_ind) = CPE_data_obj.Sigma_0;

    for (int n_ind = 0; n_ind < alpha_dmn_t::dmn_size(); n_ind++)
      alpha_func(nu_ind, nu_ind, k_ind, n_ind) = CPE_data_obj.alpha_vec_d[n_ind];

    error_func(nu_ind, nu_ind, k_ind) = evaluate_Lambda_norm(CPE_data_obj);
  }

  {  // write in the values on the imaginary axis.
    dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>& A_matrix =
        *(CPE_data_obj.A_matrix_ptr);

    func::function<std::complex<scalartype>, nu_nu_k_DCA_w_IMAG>& f_approx_func =
        *(CPE_data_obj.f_approx_ptr);

    for (int m_ind = 0; m_ind < w_IMAG::dmn_size(); m_ind++) {
      std::complex<double> value = CPE_data_obj.Sigma_0;

      for (int n_ind = 0; n_ind < alpha_dmn_t::dmn_size(); n_ind++)
        value += A_matrix(m_ind, n_ind) * CPE_data_obj.alpha_vec_d[n_ind];

      f_approx_func(nu_ind, nu_ind, k_ind, m_ind) = value;
    }
  }

  {
    func::function<std::complex<scalartype>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& S_approx_func =
        *(CPE_data_obj.S_approx_ptr);

    for (int w_ind = w::dmn_size() / 2; w_ind < w::dmn_size(); w_ind++) {
      std::complex<scalartype> value = CPE_data_obj.Sigma_0;

      double w_re = 0;
      double w_im = w::get_elements()[w_ind];

      std::complex<double> z(w_re, w_im);

      for (int n_ind = 0; n_ind < alpha_dmn_t::dmn_size(); n_ind++)
        value += basis_function_t::phi(n_ind, z) * CPE_data_obj.alpha_vec_d[n_ind];

      S_approx_func(nu_ind, nu_ind, k_ind, w_ind) = value;
    }
  }

  {
    func::function<std::complex<scalartype>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_target_func =
        *(CPE_data_obj.f_target_ptr);

    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
      std::complex<scalartype> value = CPE_data_obj.Sigma_0;

      double w_re = w_dmn_t::get_elements()[w_ind];
      double w_im = CPE_data_obj.real_axis_off_set;

      std::complex<double> z(w_re, w_im);

      for (int n_ind = 0; n_ind < alpha_dmn_t::dmn_size(); n_ind++)
        value += basis_function_t::phi(n_ind, z) * CPE_data_obj.alpha_vec_d[n_ind];

      f_target_func(nu_ind, nu_ind, k_ind, w_ind) = value;
    }
  }
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<
    parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
    WEIGHTED_GRADIENT_METHOD>::perform_continuous_pole_expansion_threaded_1(CPE_data_type& CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  CPE_data_obj.Sigma_0 = 0;
  CPE_data_obj.grad_Sigma_0 = 0;

  for (int n = 0; n < alpha_dmn_t::dmn_size(); n++)
    CPE_data_obj.alpha_vec_d[n] = 1. / double(alpha_dmn_t::dmn_size());

  int MAX_ITERATIONS = CPE_data_obj.max_iterations;
  int CPE_iteration = 0;

  gradient_minimization_state result = END_OF_PATH;

  while (MAX_ITERATIONS > CPE_iteration) {
    find_new_Sigma_0(CPE_data_obj);

    result = find_new_alpha_1(CPE_data_obj);

    if (result == BELOW_MAX_ERROR or result == FOUND_THE_MINIMUM)
      break;

    CPE_iteration++;
  }
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<
    parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
    WEIGHTED_GRADIENT_METHOD>::perform_continuous_pole_expansion_threaded_2(CPE_data_type& CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  CPE_data_obj.Sigma_0 = 0;
  CPE_data_obj.grad_Sigma_0 = 0;

  for (int n = 0; n < alpha_dmn_t::dmn_size(); n++)
    CPE_data_obj.alpha_vec_d[n] = 1. / double(alpha_dmn_t::dmn_size());

  int MAX_ITERATIONS = CPE_data_obj.max_iterations;
  int CPE_iteration = 0;

  while (MAX_ITERATIONS > CPE_iteration) {
    find_new_Sigma_0(CPE_data_obj);

    int index = find_new_alpha_2(CPE_data_obj);

    if (index == 0)
      break;

    CPE_iteration++;
  }
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::compute_gradient_Sigma_0(CPE_data_type&
                                                                                       CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  project_from_real_axis_to_imaginary_axis(CPE_data_obj.alpha_vec_d, CPE_data_obj.f_wn_re_vec,
                                           *(CPE_data_obj.A_matrix_re_ptr));

  for (int n = w_IMAG::dmn_size() / 2; n < w_IMAG::dmn_size(); n++)
    CPE_data_obj.f_wn_re_vec[n] =
        CPE_data_obj.F_wn_re_vec[n] - (CPE_data_obj.f_wn_re_vec[n] + CPE_data_obj.Sigma_0);

  CPE_data_obj.grad_Sigma_0 = 0.;
  for (int n = w_IMAG::dmn_size() / 2; n < w_IMAG::dmn_size(); n++)
    CPE_data_obj.grad_Sigma_0 += -2. * CPE_data_obj.f_wn_re_vec[n];
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::find_new_Sigma_0(CPE_data_type& CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  for (int i = 0; i < 10; i++) {
    compute_gradient_Sigma_0(CPE_data_obj);
    CPE_data_obj.Sigma_0 =
        CPE_data_obj.Sigma_0 - CPE_data_obj.grad_Sigma_0 / scalartype(2. * w_IMAG::dmn_size() / 2.);
  }
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::compute_gradient_alphas(CPE_data_type&
                                                                                      CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  //     grad_F = -2*A1*(real(spoints) - transpose(A1)*a) ...
  //               - 2*A2*(imag(spoints) - transpose(A2)*a);

  {    // compute the gradient
    {  // -2*A1*(real(spoints) - transpose(A1)*a)
      project_from_real_axis_to_imaginary_axis(CPE_data_obj.alpha_vec_d, CPE_data_obj.f_wn_re_vec,
                                               *(CPE_data_obj.A_matrix_re_ptr));

      for (int n = 0; n < w_IMAG::dmn_size(); n++)
        CPE_data_obj.f_wn_re_vec[n] =
            CPE_data_obj.F_wn_re_vec[n] - (CPE_data_obj.f_wn_re_vec[n] + CPE_data_obj.Sigma_0);

      project_from_imaginary_axis_to_real_axis(CPE_data_obj.f_wn_re_vec, CPE_data_obj.gradient_re,
                                               *(CPE_data_obj.A_matrix_re_ptr));

      for (int n = 0; n < alpha_dmn_t::dmn_size(); n++)
        CPE_data_obj.gradient_re[n] *= -2;
    }

    {  // -2*A2*(imag(spoints) - transpose(A2)*a)
      project_from_real_axis_to_imaginary_axis(CPE_data_obj.alpha_vec_d, CPE_data_obj.f_wn_im_vec,
                                               *(CPE_data_obj.A_matrix_im_ptr));

      for (int n = 0; n < w_IMAG::dmn_size(); n++)
        CPE_data_obj.f_wn_im_vec[n] = CPE_data_obj.F_wn_im_vec[n] - CPE_data_obj.f_wn_im_vec[n];

      project_from_imaginary_axis_to_real_axis(CPE_data_obj.f_wn_im_vec, CPE_data_obj.gradient_im,
                                               *(CPE_data_obj.A_matrix_im_ptr));

      for (int n = 0; n < alpha_dmn_t::dmn_size(); n++)
        CPE_data_obj.gradient_im[n] *= -2;
    }
  }

  // if(normalize)
  // normalize
  for (int n = 0; n < alpha_dmn_t::dmn_size(); n++)
    CPE_data_obj.gradient[n] = -(CPE_data_obj.gradient_re[n] + CPE_data_obj.gradient_im[n]);

  double norm_gradient = 0;
  for (int n = 0; n < alpha_dmn_t::dmn_size(); n++)
    norm_gradient += CPE_data_obj.gradient[n] * CPE_data_obj.gradient[n];

  for (int n = 0; n < alpha_dmn_t::dmn_size(); n++)
    CPE_data_obj.gradient[n] /= std::sqrt(norm_gradient);
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
typename continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::gradient_minimization_state continuous_pole_expansion<
    parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
    WEIGHTED_GRADIENT_METHOD>::find_new_alpha_1(CPE_data_type& CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  gradient_minimization_state result = END_OF_PATH;

  compute_gradient_alphas(CPE_data_obj);

  double LAMBDA_MAX = 1;

  // follow the path along the gradient
  int index = 0;
  for (index = 0; index < CPE_data_obj.y.size(); index++) {
    double lambda = 0.;
    double delta_lambda = LAMBDA_MAX / 10.;

    CPE_data_obj.x[index] = lambda + index * delta_lambda;

    for (int n = 0; n < alpha_dmn_t::dmn_size(); n++) {
      scalartype value =
          CPE_data_obj.alpha_vec_d[n] + CPE_data_obj.x[index] * CPE_data_obj.gradient[n];

      CPE_data_obj.alpha_vec_z[n].real(value < 0. ? 0. : value);
      CPE_data_obj.alpha_vec_z[n].imag(0.);
    }

    CPE_data_obj.y[index] = evaluate_Lambda_norm(CPE_data_obj);

    if (CPE_data_obj.y[index] < CPE_data_obj.max_error)
      return BELOW_MAX_ERROR;

    if (index > 1 && CPE_data_obj.y[index - 1] < CPE_data_obj.y[index]) {
      result = FOUND_A_MINIMUM;
      break;
    }
  }

  double lambda = 0;

  {  // find the lambda at the minimum
    switch (result) {
      case FOUND_A_MINIMUM: {
        assert(index >= 2 and index < CPE_data_obj.x.size());
        lambda = find_minimum_1(index, CPE_data_obj);

        if (lambda < 0) {
          return FOUND_THE_MINIMUM;
        }
        else {
          assert(CPE_data_obj.x[index - 2] - 1.e-6 < lambda and
                 lambda < CPE_data_obj.x[index - 0] + 1.e-6);
        }
      } break;

      case END_OF_PATH: {
        assert(index == CPE_data_obj.x.size());

        lambda = CPE_data_obj.x[index - 1];
      } break;

      default:
        throw std::logic_error(__FUNCTION__);
    }
  }

  compute_new_alpha(lambda, CPE_data_obj);
  return FOUND_A_MINIMUM;
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::compute_new_alpha(double lambda,
                                                                            CPE_data_type& CPE_data_obj) {
  {  // compute alpha at the minimum
    for (int n = 0; n < alpha_dmn_t::dmn_size(); n++) {
      scalartype value = CPE_data_obj.alpha_vec_d[n] + lambda * CPE_data_obj.gradient[n];

      CPE_data_obj.alpha_vec_z[n].real(value < 0. ? 0. : value);
      CPE_data_obj.alpha_vec_z[n].imag(0.);
    }
  }

  {  // smooth alpha out
    int factor = CPE_data_obj.smoothing_factor;

    for (int n = 0; n < alpha_dmn_t::dmn_size(); n++) {
      double N = 0.;
      CPE_data_obj.alpha_vec_d[n] = 0;

      for (int l = 0 - factor; l <= 0 + factor; l++) {
        if (n + l > -1 && n + l < alpha_dmn_t::dmn_size()) {
          CPE_data_obj.alpha_vec_d[n] += real(CPE_data_obj.alpha_vec_z[n + l]);
          N += 1.;
        }
      }

      CPE_data_obj.alpha_vec_d[n] /= N;
    }
  }
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
int continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                              WEIGHTED_GRADIENT_METHOD>::find_new_alpha_2(CPE_data_type& CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  compute_gradient_alphas(CPE_data_obj);

  double LAMBDA_MAX = 1;

  double lambda = 0.;
  double delta_lambda = LAMBDA_MAX / 10.;

  int index = 0;
  for (int tmp = 0; tmp < 100; tmp++) {
    CPE_data_obj.x[tmp] = lambda + tmp * delta_lambda;

    for (int n = 0; n < alpha_dmn_t::dmn_size(); n++) {
      scalartype value = CPE_data_obj.alpha_vec_d[n] + CPE_data_obj.x[tmp] * CPE_data_obj.gradient[n];

      real(CPE_data_obj.alpha_vec_z[n]) = value < 0. ? 0. : value;
      imag(CPE_data_obj.alpha_vec_z[n]) = 0.;
    }

    CPE_data_obj.y[tmp] = evaluate_Lambda_norm(CPE_data_obj);

    if (tmp > 1 && CPE_data_obj.y[tmp - 1] < CPE_data_obj.y[tmp])
      break;
    else
      index++;
  }

  if (index == 100) {
    lambda = CPE_data_obj.x[99];
  }
  else {
    dca::linalg::Matrix<scalartype, dca::linalg::CPU> V(3);

    V(0, 0) = CPE_data_obj.x[index - 2] * CPE_data_obj.x[index - 2];
    V(0, 1) = CPE_data_obj.x[index - 2];
    V(0, 2) = 1.;
    V(1, 0) = CPE_data_obj.x[index - 1] * CPE_data_obj.x[index - 1];
    V(1, 1) = CPE_data_obj.x[index - 1];
    V(1, 2) = 1.;
    V(2, 0) = CPE_data_obj.x[index - 0] * CPE_data_obj.x[index - 0];
    V(2, 1) = CPE_data_obj.x[index - 0];
    V(2, 2) = 1.;

    dca::linalg::matrixop::inverse(V);

    scalartype a = V(0, 0) * CPE_data_obj.y[index - 2] + V(0, 1) * CPE_data_obj.y[index - 1] +
                   V(0, 2) * CPE_data_obj.y[index - 0];
    scalartype b = V(1, 0) * CPE_data_obj.y[index - 2] + V(1, 1) * CPE_data_obj.y[index - 1] +
                   V(1, 2) * CPE_data_obj.y[index - 0];

    lambda = -b / (2 * a);
    lambda = lambda > 0 ? lambda : 0;
  }

  {
    for (int n = 0; n < alpha_dmn_t::dmn_size(); n++) {
      scalartype value = CPE_data_obj.alpha_vec_d[n];

      real(CPE_data_obj.alpha_vec_z[n]) = value < 0. ? 0. : value;
      imag(CPE_data_obj.alpha_vec_z[n]) = 0.;
    }

    scalartype L2_norm_0 = evaluate_Lambda_norm(CPE_data_obj);

    for (int n = 0; n < alpha_dmn_t::dmn_size(); n++) {
      scalartype value = CPE_data_obj.alpha_vec_d[n] + lambda * CPE_data_obj.gradient[n];

      real(CPE_data_obj.alpha_vec_z[n]) = value < 0. ? 0. : value;
      imag(CPE_data_obj.alpha_vec_z[n]) = 0.;
    }

    scalartype L2_norm_lambda = evaluate_Lambda_norm(CPE_data_obj);

    if (L2_norm_lambda < L2_norm_0) {
      int factor = CPE_data_obj.smoothing_factor;

      for (int n = 0; n < alpha_dmn_t::dmn_size(); n++) {
        double N = 0.;
        CPE_data_obj.alpha_vec_d[n] = 0;

        for (int l = 0 - factor; l <= 0 + factor; l++) {
          if (n + l > -1 && n + l < alpha_dmn_t::dmn_size()) {
            CPE_data_obj.alpha_vec_d[n] += real(CPE_data_obj.alpha_vec_z[n + l]);
            N += 1.;
          }
        }

        CPE_data_obj.alpha_vec_d[n] /= N;
      }
    }
    else {
      index = 0;
    }
  }

  return index;
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
double continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                                 WEIGHTED_GRADIENT_METHOD>::find_minimum_1(int index,
                                                                           CPE_data_type& CPE_data_obj) {
  double lambda = 0;

  scalartype x0 = CPE_data_obj.x[index - 2];
  scalartype x1 = CPE_data_obj.x[index - 1];
  scalartype x2 = CPE_data_obj.x[index - 0];

  scalartype y0 = CPE_data_obj.y[index - 2];
  scalartype y1 = CPE_data_obj.y[index - 1];
  scalartype y2 = CPE_data_obj.y[index - 0];

  scalartype a =
      (x2 * (-y0 + y1) + x1 * (y0 - y2) + x0 * (-y1 + y2)) / ((x0 - x1) * (x0 - x2) * (x1 - x2));
  scalartype b =
      (std::pow(x2, 2) * (y0 - y1) + std::pow(x0, 2) * (y1 - y2) + std::pow(x1, 2) * (-y0 + y2)) /
      ((x0 - x1) * (x0 - x2) * (x1 - x2));
  lambda = -b / (2 * a);

  return lambda;
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
double continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                                 WEIGHTED_GRADIENT_METHOD>::find_minimum_2(int index,
                                                                           CPE_data_type& CPE_data_obj) {
  double lambda = 0;

  dca::linalg::Matrix<scalartype, dca::linalg::CPU> V(3);

  V(0, 0) = CPE_data_obj.x[index - 2] * CPE_data_obj.x[index - 2];
  V(0, 1) = CPE_data_obj.x[index - 2];
  V(0, 2) = 1.;
  V(1, 0) = CPE_data_obj.x[index - 1] * CPE_data_obj.x[index - 1];
  V(1, 1) = CPE_data_obj.x[index - 1];
  V(1, 2) = 1.;
  V(2, 0) = CPE_data_obj.x[index - 0] * CPE_data_obj.x[index - 0];
  V(2, 1) = CPE_data_obj.x[index - 0];
  V(2, 2) = 1.;

  dca::linalg::matrixop::inverse(V);

  scalartype a = V(0, 0) * CPE_data_obj.y[index - 2] + V(0, 1) * CPE_data_obj.y[index - 1] +
                 V(0, 2) * CPE_data_obj.y[index - 0];
  scalartype b = V(1, 0) * CPE_data_obj.y[index - 2] + V(1, 1) * CPE_data_obj.y[index - 1] +
                 V(1, 2) * CPE_data_obj.y[index - 0];
  lambda = -b / (2 * a);

  return lambda;
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
double continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                                 WEIGHTED_GRADIENT_METHOD>::evaluate_Lambda_norm(CPE_data_type&
                                                                                     CPE_data_obj) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  std::complex<scalartype> MIN_ONE(-1, 0);
  std::complex<scalartype> PLUS_ONE(1, 0);

  for (int m = 0; m < w_IMAG::dmn_size(); m++)
    CPE_data_obj.f_wn_vec[m] = CPE_data_obj.F_wn_vec[m] - CPE_data_obj.Sigma_0;

  dca::linalg::matrixop::gemv('N', PLUS_ONE, *(CPE_data_obj.A_matrix_ptr), CPE_data_obj.alpha_vec_z,
                              MIN_ONE, CPE_data_obj.f_wn_vec);

  scalartype L2_norm = 0;
  for (int l = 0; l < w_IMAG::dmn_size(); l++)
    L2_norm += std::norm(CPE_data_obj.f_wn_vec[l]);

  return std::sqrt(L2_norm / (w_IMAG::dmn_size()));
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::
    project_from_real_axis_to_imaginary_axis(
        dca::linalg::Vector<scalartype, dca::linalg::CPU>& real_values,
        dca::linalg::Vector<scalartype, dca::linalg::CPU>& imag_values,
        dca::linalg::Matrix<scalartype, dca::linalg::CPU>& matrix) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  dca::linalg::matrixop::gemv('N', matrix, real_values, imag_values);
}

template <class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t,
                               WEIGHTED_GRADIENT_METHOD>::
    project_from_imaginary_axis_to_real_axis(
        dca::linalg::Vector<scalartype, dca::linalg::CPU>& imag_values,
        dca::linalg::Vector<scalartype, dca::linalg::CPU>& real_values,
        dca::linalg::Matrix<scalartype, dca::linalg::CPU>& matrix) {
  profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

  dca::linalg::matrixop::gemv('T', matrix, imag_values, real_values);
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_CONTINUOUS_POLE_EXPANSION_HPP
