//-*-C++-*-

#ifndef CONTINUOUS_POLE_EXPANSION_DATA_H
#define CONTINUOUS_POLE_EXPANSION_DATA_H

namespace DCA
{
  /*!
   *  \ingroup CPE
   *
   *  \author  Peter Staar
   *  \brief   This class implements a data-structure for the CPE analytic continuation. 
   *           It is very usefull for the multithreaded implementation.
   *  \version 1.0
   */
  template<typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  class CPE_data
  {
#include "type_definitions.h"

    typedef dmn_0<basis_function_t> alpha_dmn_t;

  public:

    CPE_data();
    CPE_data(const CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>& CPE_data_ref);

    ~CPE_data();

    template<typename parameters_type, typename function_target_t, typename CPE_solver_type>
    void initialize(int l, 
		    std::pair<int, int>& bounds_ref,
		    parameters_type&     parameters,
                    function_target_t&   f_target,
                    CPE_solver_type&     CPE_solver);

  public:

    int id;

    std::vector<bool> is_interacting_band;

    std::pair<int, int> bounds;

    int max_iterations;
    int smoothing_factor;
    
    double max_error;
    double real_axis_off_set;
    

    FUNC_LIB::function<std::complex<scalartype>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >* f_target_ptr;

    FUNC_LIB::function<scalartype, dmn_3<nu,nu,k_dmn_t>             >* Sigma_0_moment_ptr;
    FUNC_LIB::function<scalartype, dmn_4<nu,nu,k_dmn_t,alpha_dmn_t> >* alpha_function_ptr;
    FUNC_LIB::function<scalartype, dmn_3<nu,nu,k_dmn_t>             >* error_function_ptr;

    FUNC_LIB::function<std::complex<scalartype>, dmn_4<nu,nu,k_dmn_t,w> >*      S_approx_ptr;

    FUNC_LIB::function<std::complex<scalartype>, dmn_4<nu,nu,k_dmn_t,w_IMAG> >* f_approx_ptr;
    FUNC_LIB::function<std::complex<scalartype>, dmn_4<nu,nu,k_dmn_t,w_IMAG> >* f_measured_ptr;

    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>* A_matrix_ptr;

    LIN_ALG::matrix<             scalartype , LIN_ALG::CPU>* A_matrix_re_ptr;
    LIN_ALG::matrix<             scalartype , LIN_ALG::CPU>* A_matrix_im_ptr;

    scalartype      Sigma_0;
    scalartype grad_Sigma_0;

    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> alpha_vec_d;
    LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU> alpha_vec_z;

    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> gradient;
    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> gradient_re;
    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> gradient_im;

    LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU> F_wn_vec;
    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> F_wn_re_vec;
    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> F_wn_im_vec;

    LIN_ALG::vector<std::complex<scalartype>, LIN_ALG::CPU> f_wn_vec;
    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> f_wn_re_vec;
    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> f_wn_im_vec;

    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> x;
    LIN_ALG::vector<             scalartype , LIN_ALG::CPU> y;
  };

  template<typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>::CPE_data():
    id(0),

    is_interacting_band(0),

    max_iterations(100),
    smoothing_factor(0),

    max_error        (0.0),
    real_axis_off_set(0.1),

    f_target_ptr(NULL),

    Sigma_0_moment_ptr(NULL),
    alpha_function_ptr(NULL),
    error_function_ptr(NULL),

    S_approx_ptr  (NULL),

    f_approx_ptr  (NULL),
    f_measured_ptr(NULL),

    A_matrix_ptr(NULL),
    A_matrix_re_ptr(NULL),
    A_matrix_im_ptr(NULL),

    Sigma_0(0),
    grad_Sigma_0(0),

    alpha_vec_d("alpha_vec_d", alpha_dmn_t::dmn_size()),
    alpha_vec_z("alpha_vec_z", alpha_dmn_t::dmn_size()),

    gradient   ("gradient"   , alpha_dmn_t::dmn_size()),
    gradient_re("gradient_re", alpha_dmn_t::dmn_size()),
    gradient_im("gradient_im", alpha_dmn_t::dmn_size()),

    // measured
    F_wn_vec   ("F_wn_vec"   , w_IMAG::dmn_size()),
    F_wn_re_vec("F_wn_re_vec", w_IMAG::dmn_size()),
    F_wn_im_vec("F_wn_im_vec", w_IMAG::dmn_size()),

    // approx
    f_wn_vec   ("f_wn_vec"   , w_IMAG::dmn_size()),
    f_wn_re_vec("f_wn_re_vec", w_IMAG::dmn_size()),
    f_wn_im_vec("f_wn_im_vec", w_IMAG::dmn_size()),

    x("x", 128),
    y("y", 128)
  {
  }

  template<typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>::CPE_data(const CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>& /*CPE_data_ref*/):
    id(0),

    is_interacting_band(0),

    max_iterations(100),
    smoothing_factor(0),

    max_error        (0.0),
    real_axis_off_set(0.1),

    f_target_ptr(NULL),

    Sigma_0_moment_ptr(NULL),
    alpha_function_ptr(NULL),
    error_function_ptr(NULL),

    f_approx_ptr  (NULL),
    f_measured_ptr(NULL),

    A_matrix_ptr(NULL),
    A_matrix_re_ptr(NULL),
    A_matrix_im_ptr(NULL),

    Sigma_0(0),
    grad_Sigma_0(0),

    alpha_vec_d("alpha_vec_d", alpha_dmn_t::dmn_size()),
    alpha_vec_z("alpha_vec_z", alpha_dmn_t::dmn_size()),

    gradient   ("gradient"   , alpha_dmn_t::dmn_size()),
    gradient_re("gradient_re", alpha_dmn_t::dmn_size()),
    gradient_im("gradient_im", alpha_dmn_t::dmn_size()),

    // measured
    F_wn_vec   ("F_wn_vec"   , w_IMAG::dmn_size()),
    F_wn_re_vec("F_wn_re_vec", w_IMAG::dmn_size()),
    F_wn_im_vec("F_wn_im_vec", w_IMAG::dmn_size()),

    // approx
    f_wn_vec   ("f_wn_vec"   , w_IMAG::dmn_size()),
    f_wn_re_vec("f_wn_re_vec", w_IMAG::dmn_size()),
    f_wn_im_vec("f_wn_im_vec", w_IMAG::dmn_size()),

    x("x", 128),
    y("y", 128)
  {}

  template<typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>::~CPE_data()
  {}

  template<typename scalartype, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  template<typename parameters_type, typename function_target_t, typename CPE_solver_type>
  void CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t>::initialize(int l, 
									    std::pair<int, int>& bounds_ref,
									    parameters_type&     parameters,
                                                                            function_target_t&   f_target,
                                                                            CPE_solver_type&     CPE_solver)
  {
    id = l;

    for(int b_ind=0; b_ind<b::dmn_size(); b_ind++)
      for(int s_ind=0; s_ind<b::dmn_size(); s_ind++)
	is_interacting_band.push_back(parameters.is_an_interacting_band(b_ind));

    bounds = bounds_ref;

    max_iterations   = parameters.get_max_CPE_iterations();
    smoothing_factor = parameters.get_CPE_smoothing_factor();

    max_error         = parameters.get_max_CPE_error();
    real_axis_off_set = parameters.get_real_frequencies_off_set();

    f_target_ptr = &f_target;

    Sigma_0_moment_ptr = &(CPE_solver.Sigma_0_moment);
    alpha_function_ptr = &(CPE_solver.alpha_function);
    error_function_ptr = &(CPE_solver.error_function);

    S_approx_ptr   = &(CPE_solver.S_approx);

    f_approx_ptr   = &(CPE_solver.f_approx);
    f_measured_ptr = &(CPE_solver.f_measured);

    A_matrix_ptr = &(CPE_solver.A_matrix);

    A_matrix_re_ptr = &(CPE_solver.A_matrix_re);
    A_matrix_im_ptr = &(CPE_solver.A_matrix_im);
  }

}

#endif
