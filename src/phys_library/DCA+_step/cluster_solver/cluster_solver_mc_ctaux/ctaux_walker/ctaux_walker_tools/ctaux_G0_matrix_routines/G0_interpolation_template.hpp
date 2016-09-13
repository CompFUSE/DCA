// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the interpolation of \f$G^{0}\f$ towards the \f$G^{0}\f$-matrix.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G0_MATRIX_ROUTINES_G0_INTERPOLATION_TEMPLATE_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G0_MATRIX_ROUTINES_G0_INTERPOLATION_TEMPLATE_HPP

#include "dca/linalg/matrix.hpp"

#include "comp_library/function_library/include_function_library.h"
#include "math_library/interpolation_library/akima_interpolation.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"
#include "phys_library/domains/time_and_frequency/time_domain_left_oriented.h"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <typename parameters_type>
class G0_INTERPOLATION_TEMPLATE {
public:
  using t = dmn_0<time_domain>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index

  using r_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     REAL_SPACE, BRILLOUIN_ZONE>>;
  using r_dmn_t = r_DCA;
  typedef typename r_dmn_t::parameter_type r_cluster_type;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

  typedef dmn_0<time_domain_left_oriented> shifted_t;
  typedef dmn_4<nu, nu, r_dmn_t, shifted_t> nu_nu_r_dmn_t_shifted_t;

  typedef dmn_0<dmn<4, int>> akima_dmn_t;
  typedef dmn_5<akima_dmn_t, nu, nu, r_dmn_t, shifted_t> akima_nu_nu_r_dmn_t_shifted_t;

public:
  G0_INTERPOLATION_TEMPLATE(int id, parameters_type& parameters);

  template <class MOMS_type>
  void initialize(MOMS_type& MOMS);

protected:
  template <class MOMS_type>
  void initialize_linear_coefficients(MOMS_type& MOMS);

  template <class MOMS_type>
  void initialize_akima_coefficients(MOMS_type& MOMS);

protected:
  int thread_id;

  parameters_type& parameters;
  concurrency_type& concurrency;

  nu_nu_r_dmn_t_shifted_t nu_nu_r_dmn_t_t_shifted_dmn;

  dca::linalg::Matrix<double, dca::linalg::CPU> r1_minus_r0;

  FUNC_LIB::function<double, nu_nu_r_dmn_t_shifted_t> G0_r_t_shifted;
  FUNC_LIB::function<double, nu_nu_r_dmn_t_shifted_t> grad_G0_r_t_shifted;

  FUNC_LIB::function<double, akima_nu_nu_r_dmn_t_shifted_t> akima_coefficients;

  int N_t, linind, t_ind;
  double beta, N_div_beta, new_tau, scaled_tau, delta_tau, f_0, grad;
};

template <typename parameters_type>
G0_INTERPOLATION_TEMPLATE<parameters_type>::G0_INTERPOLATION_TEMPLATE(int id,
                                                                      parameters_type& parameters_ref)
    : thread_id(id),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      r1_minus_r0(r_dmn_t::dmn_size()) {
  beta = parameters.get_beta();

  // N_t        = parameters.get_number_of_positive_times()+1;
  N_t = parameters.get_sp_time_intervals() + 1;
  N_div_beta = parameters.get_sp_time_intervals() / beta;

  for (int r1_ind = 0; r1_ind < r_dmn_t::dmn_size(); r1_ind++)
    for (int r0_ind = 0; r0_ind < r_dmn_t::dmn_size(); r0_ind++)
      r1_minus_r0(r0_ind, r1_ind) = r_cluster_type::subtract(r0_ind, r1_ind);
}

/*!
 *  \brief  Set the functions 'G0_r_t_shifted' and 'grad_G0_r_t_shifted'
 */
template <typename parameters_type>
template <class MOMS_type>
void G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize(MOMS_type& MOMS) {
  initialize_linear_coefficients(MOMS);

  initialize_akima_coefficients(MOMS);
}

template <typename parameters_type>
template <class MOMS_type>
void G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize_linear_coefficients(MOMS_type& MOMS) {
  for (int t_ind = 0; t_ind < t::dmn_size() / 2 - 1; t_ind++) {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < b::dmn_size() * s::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < b::dmn_size() * s::dmn_size(); nu0_ind++) {
          G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind) =
              MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);
          grad_G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind) =
              (MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind + 1) -
               MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind));
        }
      }
    }
  }

  for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size() - 1; t_ind++) {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < b::dmn_size() * s::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < b::dmn_size() * s::dmn_size(); nu0_ind++) {
          G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind - 1) =
              MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);
          grad_G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind - 1) =
              (MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind + 1) -
               MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind));
        }
      }
    }
  }
}

template <typename parameters_type>
template <class MOMS_type>
void G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize_akima_coefficients(MOMS_type& MOMS) {
  int size = t::dmn_size() / 2;

  math_algorithms::interpolation::akima_interpolation<double> ai_obj(size);

  double* x = new double[size];
  double* y = new double[size];

  for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
    x[t_ind] = t_ind;

  {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < b::dmn_size() * s::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < b::dmn_size() * s::dmn_size(); nu0_ind++) {
          for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
            y[t_ind] = MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);

          ai_obj.initialize(x, y);

          for (int t_ind = 0; t_ind < t::dmn_size() / 2 - 1; t_ind++)
            for (int l = 0; l < 4; l++)
              akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind) = ai_obj.get_alpha(l, t_ind);
        }
      }
    }
  }

  {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < b::dmn_size() * s::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < b::dmn_size() * s::dmn_size(); nu0_ind++) {
          for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size(); t_ind++)
            y[t_ind - t::dmn_size() / 2] =
                MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);

          ai_obj.initialize(x, y);

          for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size() - 1; t_ind++)
            for (int l = 0; l < 4; l++)
              akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind - 1) =
                  ai_obj.get_alpha(l, t_ind - t::dmn_size() / 2);
        }
      }
    }
  }

  delete[] x;
  delete[] y;

  /*
    {
    cout << "\n\n\n";
    for(int t_ind=0; t_ind<shifted_t::dmn_size(); t_ind++)
    cout << t_ind << "\t" << G0_r_t_shifted(0, 0, 0, t_ind) << endl;
    cout << "\n\n\n";

    for(int t_ind=0; t_ind<shifted_t::dmn_size(); t_ind++)
    {
    int linind    = 4*nu_nu_r_dmn_t_t_shifted_dmn(0,0,0,t_ind);
    double* a_ptr = &akima_coefficents(linind);

    for(double x=0; x<1.05; x+=0.1)
    cout << t_ind+x << "\t" << (a_ptr[0] + x*(a_ptr[1] + x*(a_ptr[2] + x*a_ptr[3]))) << endl;
    }


    sleep(10);
    throw std::logic_error(__FUNCTION__);
    }
  */
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G0_MATRIX_ROUTINES_G0_INTERPOLATION_TEMPLATE_HPP
