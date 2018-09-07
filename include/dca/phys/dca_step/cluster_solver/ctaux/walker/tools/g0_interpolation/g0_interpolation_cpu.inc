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
// Template specialization for CPU.

template <typename parameters_type>
class G0_INTERPOLATION<dca::linalg::CPU, parameters_type>
    : public G0_INTERPOLATION_TEMPLATE<parameters_type> {
public:
  using vertex_singleton_type = vertex_singleton;
  using shifted_t = func::dmn_0<domains::time_domain_left_oriented>;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;

  using r_dmn_t = RClusterDmn;
  using r_cluster_type = typename r_dmn_t::parameter_type;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  G0_INTERPOLATION(int id, parameters_type& parameters);

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize;

  template <class configuration_type>
  void build_G0_matrix(configuration_type& configuration,
                       dca::linalg::Matrix<double, dca::linalg::CPU>& G0, e_spin_states_type spin);

  void build_G0_matrix(std::vector<vertex_singleton_type>& configuration,
                       dca::linalg::Matrix<double, dca::linalg::CPU>& G0);

  template <class configuration_type>
  void update_G0_matrix(configuration_type& configuration,
                        dca::linalg::Matrix<double, dca::linalg::CPU>& G0, e_spin_states_type spin);

  std::size_t deviceFingerprint() const {return 0;}

private:
  double interpolate(int nu_0, int nu_1, int delta_r, double delta_time);

  double interpolate_on_diagonal(int nu_i);

  double interpolate_akima(int nu_0, int nu_1, int delta_r, double tau);

private:
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::parameters;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::concurrency;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::nu_nu_r_dmn_t_t_shifted_dmn;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::r1_minus_r0;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::G0_r_t_shifted;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::grad_G0_r_t_shifted;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::akima_coefficients;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::N_t;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::linind;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::t_ind;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::beta;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::N_div_beta;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::new_tau;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::scaled_tau;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::delta_tau;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::f_0;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::grad;

  int thread_id;
};

template <typename parameters_type>
G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::G0_INTERPOLATION(int id,
                                                                      parameters_type& parameters_ref)
    : G0_INTERPOLATION_TEMPLATE<parameters_type>(id, parameters_ref),

      thread_id(id) {}

template <typename parameters_type>
template <class configuration_type>
void G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::build_G0_matrix(
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

template <typename parameters_type>
void G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::build_G0_matrix(
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

template <typename parameters_type>
template <class configuration_type>
void G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::update_G0_matrix(
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

template <typename parameters_type>
inline double G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::interpolate(int nu_0, int nu_1,
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

  linind = nu_nu_r_dmn_t_t_shifted_dmn(nu_0, nu_1, delta_r, t_ind);

  f_0 = G0_r_t_shifted(linind);
  grad = grad_G0_r_t_shifted(linind);

  return -(f_0 + grad * delta_tau);
}

template <typename parameters_type>
inline double G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::interpolate_on_diagonal(int nu_i) {
  const static int t_0_index = shifted_t::dmn_size() / 2;
  const static int r_0_index = r_cluster_type::origin_index();

  return -G0_r_t_shifted(nu_i, nu_i, r_0_index, t_0_index);
}

template <typename parameters_type>
inline double G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::interpolate_akima(int nu_0,
                                                                                     int nu_1,
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

/*
  template<typename parameters_type>
  inline double G0_INTERPOLATION<dca::linalg::CPU, parameters_type>::interpolate(int nu_0, int nu_1,
  int
  delta_r, double tau)
  {
  new_tau = tau+beta;

  scaled_tau = new_tau*N_div_beta;

  t_ind = scaled_tau;
  assert(shifted_t::get_elements()[t_ind]<=tau &&
  tau<shifted_t::get_elements()[t_ind]+1./N_div_beta);

  delta_tau = scaled_tau-t_ind;
  assert(delta_tau > -1.e-6 && delta_tau <= 1);

  linind = nu_nu_r_dmn_t_t_shifted_dmn(nu_0,nu_1,delta_r,t_ind);

  f_0  =      G0_r_t_shifted(linind);
  grad = grad_G0_r_t_shifted(linind);

  return -(f_0 + grad*delta_tau);
  }
*/
