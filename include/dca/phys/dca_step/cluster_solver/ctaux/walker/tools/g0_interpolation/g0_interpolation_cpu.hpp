// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.eth.ch)
//
// This class organizes the interpolation of \f$G^{0}\f$ towards the \f$G^{0}\f$-matrix.
// Template specialization for CPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_CPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_CPU_HPP

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_base.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_tmpl.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <typename Parameters>
class G0Interpolation<dca::linalg::CPU, Parameters> : public G0InterpolationBase<Parameters> {
public:
  using vertex_singleton_type = vertex_singleton;
  using shifted_t = func::dmn_0<domains::time_domain_left_oriented>;

  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;

  using r_dmn_t = RClusterDmn;
  using r_cluster_type = typename r_dmn_t::parameter_type;

  typedef typename Parameters::concurrency_type concurrency_type;
  typedef typename Parameters::profiler_type profiler_t;

public:
  G0Interpolation(int id, Parameters& parameters);

  using G0InterpolationBase<Parameters>::initialize;

  template <class configuration_type>
  void build_G0_matrix(configuration_type& configuration,
                       dca::linalg::Matrix<double, dca::linalg::CPU>& G0, e_spin_states_type spin);

  void build_G0_matrix(std::vector<vertex_singleton_type>& configuration,
                       dca::linalg::Matrix<double, dca::linalg::CPU>& G0);

  template <class configuration_type>
  void update_G0_matrix(configuration_type& configuration,
                        dca::linalg::Matrix<double, dca::linalg::CPU>& G0, e_spin_states_type spin);

  std::size_t deviceFingerprint() const {
    return 0;
  }

private:
  double interpolate(int nu_0, int nu_1, int delta_r, double delta_time);

  double interpolate_on_diagonal(int nu_i);

  double interpolate_akima(int nu_0, int nu_1, int delta_r, double tau);

private:
  using G0InterpolationBase<Parameters>::parameters;
  using G0InterpolationBase<Parameters>::concurrency;

  using G0InterpolationBase<Parameters>::nu_nu_r_dmn_t_t_shifted_dmn;

  using G0InterpolationBase<Parameters>::r1_minus_r0;

  using G0InterpolationBase<Parameters>::G0_r_t_shifted;
  using G0InterpolationBase<Parameters>::grad_G0_r_t_shifted;

  using G0InterpolationBase<Parameters>::akima_coefficients;

  using G0InterpolationBase<Parameters>::N_t;
  using G0InterpolationBase<Parameters>::linind;
  using G0InterpolationBase<Parameters>::t_ind;

  using G0InterpolationBase<Parameters>::beta;
  using G0InterpolationBase<Parameters>::N_div_beta;
  using G0InterpolationBase<Parameters>::new_tau;
  using G0InterpolationBase<Parameters>::scaled_tau;
  using G0InterpolationBase<Parameters>::delta_tau;
  using G0InterpolationBase<Parameters>::f_0;
  using G0InterpolationBase<Parameters>::grad;

  int thread_id;
};

template <typename Parameters>
G0Interpolation<dca::linalg::CPU, Parameters>::G0Interpolation(int id, Parameters& parameters_ref)
    : G0InterpolationBase<Parameters>(id, parameters_ref),

      thread_id(id) {}

template <typename Parameters>
template <class configuration_type>
void G0Interpolation<dca::linalg::CPU, Parameters>::build_G0_matrix(
    configuration_type& configuration, dca::linalg::Matrix<double, dca::linalg::CPU>& G0_e_spin,
    e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  G0_e_spin.resizeNoCopy(configuration_size);

  for (int j = 0; j < configuration_size; j++) {
    vertex_singleton_type& v_j = configuration_e_spin[j];

    {  // i < j
      for (int i = 0; i < j; i++) {
        vertex_singleton_type& v_i = configuration_e_spin[i];

        G0_e_spin(i, j) = interpolate_akima(v_i.get_spin_orbital(), v_j.get_spin_orbital(),
                                            r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
                                            (v_i.get_tau() - v_j.get_tau()));
      }
    }

    {  // i == j
      G0_e_spin(j, j) = interpolate_on_diagonal(v_j.get_spin_orbital());
    }

    {  // j > i
      for (int i = j + 1; i < configuration_size; i++) {
        vertex_singleton_type& v_i = configuration_e_spin[i];

        G0_e_spin(i, j) = interpolate_akima(v_i.get_spin_orbital(), v_j.get_spin_orbital(),
                                            r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
                                            (v_i.get_tau() - v_j.get_tau()));
      }
    }
  }
}

template <typename Parameters>
void G0Interpolation<dca::linalg::CPU, Parameters>::build_G0_matrix(
    std::vector<vertex_singleton_type>& configuration,
    dca::linalg::Matrix<double, dca::linalg::CPU>& G0) {
  // profiler_t profiler(concurrency, "G0-matrix (build)", "CT-AUX", __LINE__);

  int configuration_size = configuration.size();

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  G0.resizeNoCopy(configuration_size);

  for (int j = 0; j < configuration_size; j++) {
    vertex_singleton_type& v_j = configuration[j];

    {  // i < j
      for (int i = 0; i < j; i++) {
        vertex_singleton_type& v_i = configuration[i];

        G0(i, j) = interpolate_akima(v_i.get_spin_orbital(), v_j.get_spin_orbital(),
                                     r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
                                     (v_i.get_tau() - v_j.get_tau()));
      }
    }

    {  // i == j
      G0(j, j) = interpolate_on_diagonal(v_j.get_spin_orbital());
    }

    {  // j > i
      for (int i = j + 1; i < configuration_size; i++) {
        vertex_singleton_type& v_i = configuration[i];

        G0(i, j) = interpolate_akima(v_i.get_spin_orbital(), v_j.get_spin_orbital(),
                                     r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
                                     (v_i.get_tau() - v_j.get_tau()));
      }
    }
  }
}

template <typename Parameters>
template <class configuration_type>
void G0Interpolation<dca::linalg::CPU, Parameters>::update_G0_matrix(
    configuration_type& configuration, dca::linalg::Matrix<double, dca::linalg::CPU>& G0,
    e_spin_states_type e_spin) {
  // profiler_t profiler("G0-matrix (update)", "CT-AUX", __LINE__);

  std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  G0.resize(configuration_size);

  int first_shuffled_index = configuration.get_first_shuffled_spin_index(e_spin);

  for (int j = 0; j < first_shuffled_index; j++) {
    double* G0_ptr = G0.ptr(0, j);

    vertex_singleton_type& v_j = configuration_e_spin[j];

    for (int i = first_shuffled_index; i < configuration_size; i++) {
      assert(i >= first_shuffled_index || j >= first_shuffled_index);

      vertex_singleton_type& v_i = configuration_e_spin[i];

      G0_ptr[i] = interpolate_akima(v_i.get_spin_orbital(), v_j.get_spin_orbital(),
                                    r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
                                    (v_i.get_tau() - v_j.get_tau()));
    }
  }

  for (int j = first_shuffled_index; j < configuration_size; j++) {
    double* G0_ptr = G0.ptr(0, j);

    vertex_singleton_type& v_j = configuration_e_spin[j];

    // i<j
    for (int i = 0; i < j; i++) {
      assert(i >= first_shuffled_index || j >= first_shuffled_index);

      vertex_singleton_type& v_i = configuration_e_spin[i];

      G0_ptr[i] = interpolate_akima(v_i.get_spin_orbital(), v_j.get_spin_orbital(),
                                    r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
                                    (v_i.get_tau() - v_j.get_tau()));
    }

    {  // i == j
      G0_ptr[j] = interpolate_on_diagonal(v_j.get_spin_orbital());
    }

    // i>j
    for (int i = j + 1; i < configuration_size; i++) {
      assert(i >= first_shuffled_index || j >= first_shuffled_index);

      vertex_singleton_type& v_i = configuration_e_spin[i];

      G0_ptr[i] = interpolate_akima(v_i.get_spin_orbital(), v_j.get_spin_orbital(),
                                    r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
                                    (v_i.get_tau() - v_j.get_tau()));
    }
  }

  /*
    for(int j=first_shuffled_index; j<configuration_size; j++){

    double* G0_ptr = G0.ptr(0,j);

    vertex_singleton_type& v_j = configuration_e_spin[j];

    for(int i=0; i<configuration_size; i++){
    assert(i >= first_shuffled_index || j >= first_shuffled_index);

    vertex_singleton_type& v_i = configuration_e_spin[i];

    G0_ptr[i] = interpolate(v_i.get_spin_orbital(),v_j.get_spin_orbital(),
    r1_minus_r0(v_j.get_r_site(), v_i.get_r_site()),
    (v_i.get_tau()-v_j.get_tau()));
    }
    }
  */
}

template <typename Parameters>
inline double G0Interpolation<dca::linalg::CPU, Parameters>::interpolate(int nu_0, int nu_1,
                                                                         int delta_r, double tau) {
  // make sure that new_tau is positive !!
  new_tau = tau + beta;

  scaled_tau = new_tau * N_div_beta;

  t_ind = scaled_tau;
  assert(shifted_t::get_elements()[t_ind] <= tau &&
         tau < shifted_t::get_elements()[t_ind] + 1. / N_div_beta);

  delta_tau = scaled_tau - t_ind;
  assert(delta_tau > -1.e-16 && delta_tau <= 1 + 1.e-16);

  linind = nu_nu_r_dmn_t_t_shifted_dmn(nu_0, nu_1, delta_r, t_ind);

  f_0 = G0_r_t_shifted(linind);
  grad = grad_G0_r_t_shifted(linind);

  return -(f_0 + grad * delta_tau);
}

template <typename Parameters>
inline double G0Interpolation<dca::linalg::CPU, Parameters>::interpolate_on_diagonal(int nu_i) {
  const static int t_0_index = shifted_t::dmn_size() / 2;
  const static int r_0_index = r_cluster_type::origin_index();

  return -G0_r_t_shifted(nu_i, nu_i, r_0_index, t_0_index);
}

template <typename Parameters>
inline double G0Interpolation<dca::linalg::CPU, Parameters>::interpolate_akima(int nu_0, int nu_1,
                                                                               int delta_r,
                                                                               double tau) {
  // make sure that new_tau is positive !!
  new_tau = tau + beta;

  scaled_tau = new_tau * N_div_beta;

  t_ind = scaled_tau;
  assert(shifted_t::get_elements()[t_ind] <= tau &&
         tau < shifted_t::get_elements()[t_ind] + 1. / N_div_beta);

  delta_tau = scaled_tau - t_ind;
  assert(delta_tau > -1.e-16 && delta_tau <= 1 + 1.e-16);

  linind = 4 * nu_nu_r_dmn_t_t_shifted_dmn(nu_0, nu_1, delta_r, t_ind);

  double* a_ptr = &akima_coefficients(linind);

  return -(a_ptr[0] + delta_tau * (a_ptr[1] + delta_tau * (a_ptr[2] + delta_tau * a_ptr[3])));
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_CPU_HPP
