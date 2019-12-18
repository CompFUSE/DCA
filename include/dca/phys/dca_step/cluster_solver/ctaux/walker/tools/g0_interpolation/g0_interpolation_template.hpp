// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the interpolation of \f$G^{0}\f$ towards the \f$G^{0}\f$-matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_TEMPLATE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_TEMPLATE_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/interpolation/akima_interpolation.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <typename Parameters, typename Real>
class G0_INTERPOLATION_TEMPLATE {
public:
  using t = func::dmn_0<domains::time_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;

  using r_dmn_t = RClusterDmn;
  typedef typename r_dmn_t::parameter_type r_cluster_type;

  typedef typename Parameters::concurrency_type concurrency_type;
  typedef typename Parameters::profiler_type profiler_t;

  typedef func::dmn_0<domains::time_domain_left_oriented> shifted_t;
  typedef func::dmn_variadic<nu, nu, r_dmn_t, shifted_t> nu_nu_r_dmn_t_shifted_t;

  typedef func::dmn_0<func::dmn<4, int>> akima_dmn_t;
  typedef func::dmn_variadic<akima_dmn_t, nu, nu, r_dmn_t, shifted_t> akima_nu_nu_r_dmn_t_shifted_t;

public:
  G0_INTERPOLATION_TEMPLATE(int id, Parameters& parameters);

  template <class MOMS_type>
  void initialize(MOMS_type& MOMS);

protected:
  template <class MOMS_type>
  void initialize_linear_coefficients(MOMS_type& MOMS);

  template <class MOMS_type>
  void initialize_akima_coefficients(MOMS_type& MOMS);

protected:
  int thread_id;

  Parameters& parameters;
  concurrency_type& concurrency;

  nu_nu_r_dmn_t_shifted_t nu_nu_r_dmn_t_t_shifted_dmn;

  dca::linalg::Matrix<Real, dca::linalg::CPU> r1_minus_r0;

  func::function<Real, nu_nu_r_dmn_t_shifted_t> G0_r_t_shifted;
  func::function<Real, nu_nu_r_dmn_t_shifted_t> grad_G0_r_t_shifted;

  func::function<Real, akima_nu_nu_r_dmn_t_shifted_t> akima_coefficients;

  int N_t, linind, t_ind;
  Real beta, N_div_beta, new_tau, scaled_tau, delta_tau, f_0, grad;
};

template <typename Parameters, typename Real>
G0_INTERPOLATION_TEMPLATE<Parameters, Real>::G0_INTERPOLATION_TEMPLATE(int id,
                                                                       Parameters& parameters_ref)
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
template <typename Parameters, typename Real>
template <class MOMS_type>
void G0_INTERPOLATION_TEMPLATE<Parameters, Real>::initialize(MOMS_type& MOMS) {
  initialize_linear_coefficients(MOMS);

  initialize_akima_coefficients(MOMS);
}

template <typename Parameters, typename Real>
template <class MOMS_type>
void G0_INTERPOLATION_TEMPLATE<Parameters, Real>::initialize_linear_coefficients(MOMS_type& MOMS) {
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

template <typename Parameters, typename Real>
template <class MOMS_type>
void G0_INTERPOLATION_TEMPLATE<Parameters, Real>::initialize_akima_coefficients(MOMS_type& MOMS) {
  int size = t::dmn_size() / 2;

  math::interpolation::akima_interpolation<Real> ai_obj(size);

  Real* x = new Real[size];
  Real* y = new Real[size];

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
    Real* a_ptr = &akima_coefficents(linind);

    for(Real x=0; x<1.05; x+=0.1)
    cout << t_ind+x << "\t" << (a_ptr[0] + x*(a_ptr[1] + x*(a_ptr[2] + x*a_ptr[3]))) << endl;
    }


    sleep(10);
    throw std::logic_error(__FUNCTION__);
    }
  */
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_TEMPLATE_HPP
