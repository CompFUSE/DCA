// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements a data structure for the CPE analytic continuation. It is used for the
// multithreaded implementation.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_CPE_DATA_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_CPE_DATA_HPP

#include <complex>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_imag_axis.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
class CPE_data {
public:
  using alpha_dmn_t = func::dmn_0<basis_function_t>;

  using w = func::dmn_0<domains::frequency_domain>;
  using w_IMAG = func::dmn_0<domains::frequency_domain_imag_axis>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  CPE_data();

  template <typename parameters_type, typename function_target_t, typename CPE_solver_type>
  void initialize(int l, std::pair<int, int>& bounds_ref, parameters_type& parameters,
                  function_target_t& f_target, CPE_solver_type& CPE_solver);

  int id;

  std::vector<bool> is_interacting_band;

  std::pair<int, int> bounds;

  int max_iterations;
  int smoothing_factor;

  double max_error;
  double real_axis_off_set;

  func::function<std::complex<scalartype>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>* f_target_ptr;

  func::function<scalartype, func::dmn_variadic<nu, nu, k_dmn_t>>* Sigma_0_moment_ptr;
  func::function<scalartype, func::dmn_variadic<nu, nu, k_dmn_t, alpha_dmn_t>>* alpha_function_ptr;
  func::function<scalartype, func::dmn_variadic<nu, nu, k_dmn_t>>* error_function_ptr;

  func::function<std::complex<scalartype>, func::dmn_variadic<nu, nu, k_dmn_t, w>>* S_approx_ptr;

  func::function<std::complex<scalartype>, func::dmn_variadic<nu, nu, k_dmn_t, w_IMAG>>* f_approx_ptr;
  func::function<std::complex<scalartype>, func::dmn_variadic<nu, nu, k_dmn_t, w_IMAG>>* f_measured_ptr;

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU>* A_matrix_ptr;

  dca::linalg::Matrix<scalartype, dca::linalg::CPU>* A_matrix_re_ptr;
  dca::linalg::Matrix<scalartype, dca::linalg::CPU>* A_matrix_im_ptr;

  scalartype Sigma_0;
  scalartype grad_Sigma_0;

  dca::linalg::Vector<scalartype, dca::linalg::CPU> alpha_vec_d;
  dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU> alpha_vec_z;

  dca::linalg::Vector<scalartype, dca::linalg::CPU> gradient;
  dca::linalg::Vector<scalartype, dca::linalg::CPU> gradient_re;
  dca::linalg::Vector<scalartype, dca::linalg::CPU> gradient_im;

  dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU> F_wn_vec;
  dca::linalg::Vector<scalartype, dca::linalg::CPU> F_wn_re_vec;
  dca::linalg::Vector<scalartype, dca::linalg::CPU> F_wn_im_vec;

  dca::linalg::Vector<std::complex<scalartype>, dca::linalg::CPU> f_wn_vec;
  dca::linalg::Vector<scalartype, dca::linalg::CPU> f_wn_re_vec;
  dca::linalg::Vector<scalartype, dca::linalg::CPU> f_wn_im_vec;

  dca::linalg::Vector<scalartype, dca::linalg::CPU> x;
  dca::linalg::Vector<scalartype, dca::linalg::CPU> y;
};

template <typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>::CPE_data()
    : id(0),

      is_interacting_band(0),

      max_iterations(100),
      smoothing_factor(0),

      max_error(0.0),
      real_axis_off_set(0.1),

      f_target_ptr(nullptr),

      Sigma_0_moment_ptr(nullptr),
      alpha_function_ptr(nullptr),
      error_function_ptr(nullptr),

      S_approx_ptr(nullptr),

      f_approx_ptr(nullptr),
      f_measured_ptr(nullptr),

      A_matrix_ptr(nullptr),
      A_matrix_re_ptr(nullptr),
      A_matrix_im_ptr(nullptr),

      Sigma_0(0),
      grad_Sigma_0(0),

      alpha_vec_d("alpha_vec_d", alpha_dmn_t::dmn_size()),
      alpha_vec_z("alpha_vec_z", alpha_dmn_t::dmn_size()),

      gradient("gradient", alpha_dmn_t::dmn_size()),
      gradient_re("gradient_re", alpha_dmn_t::dmn_size()),
      gradient_im("gradient_im", alpha_dmn_t::dmn_size()),

      // measured
      F_wn_vec("F_wn_vec", w_IMAG::dmn_size()),
      F_wn_re_vec("F_wn_re_vec", w_IMAG::dmn_size()),
      F_wn_im_vec("F_wn_im_vec", w_IMAG::dmn_size()),

      // approx
      f_wn_vec("f_wn_vec", w_IMAG::dmn_size()),
      f_wn_re_vec("f_wn_re_vec", w_IMAG::dmn_size()),
      f_wn_im_vec("f_wn_im_vec", w_IMAG::dmn_size()),

      x("x", 128),
      y("y", 128) {}

template <typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
template <typename parameters_type, typename function_target_t, typename CPE_solver_type>
void CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>::initialize(
    int l, std::pair<int, int>& bounds_ref, parameters_type& parameters,
    function_target_t& f_target, CPE_solver_type& CPE_solver) {
  id = l;

  for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++)
    for (int s_ind = 0; s_ind < b::dmn_size(); s_ind++)
      is_interacting_band.push_back(parameters.is_an_interacting_band(b_ind));

  bounds = bounds_ref;

  max_iterations = parameters.get_max_CPE_iterations();
  smoothing_factor = parameters.get_CPE_smoothing_factor();

  max_error = parameters.get_max_CPE_error();
  real_axis_off_set = parameters.get_real_frequencies_off_set();

  f_target_ptr = &f_target;

  Sigma_0_moment_ptr = &(CPE_solver.Sigma_0_moment);
  alpha_function_ptr = &(CPE_solver.alpha_function);
  error_function_ptr = &(CPE_solver.error_function);

  S_approx_ptr = &(CPE_solver.S_approx);

  f_approx_ptr = &(CPE_solver.f_approx);
  f_measured_ptr = &(CPE_solver.f_measured);

  A_matrix_ptr = &(CPE_solver.A_matrix);

  A_matrix_re_ptr = &(CPE_solver.A_matrix_re);
  A_matrix_im_ptr = &(CPE_solver.A_matrix_im);
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_CPE_DATA_HPP
