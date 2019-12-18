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
// Template specialization for GPU.

template <typename Parameters, typename Real>
class G0_INTERPOLATION<dca::linalg::GPU, Parameters, Real>
    : public G0_INTERPOLATION_TEMPLATE<Parameters, Real> {
public:
  using vertex_singleton_type = vertex_singleton;
  using shifted_t = func::dmn_0<domains::time_domain_left_oriented>;
  using b = func::dmn_0<domains::electron_band_domain>;

  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using r_dmn_t = RClusterDmn;

  typedef typename Parameters::concurrency_type concurrency_type;
  typedef typename Parameters::profiler_type profiler_t;

public:
  G0_INTERPOLATION(int id, Parameters& parameters);

  template <class MOMS_type>
  void initialize(MOMS_type& MOMS);

  template <class configuration_type>
  void build_G0_matrix(configuration_type& configuration,
                       dca::linalg::Matrix<Real, dca::linalg::GPU>& G0, e_spin_states_type spin);

  template <class configuration_type>
  void update_G0_matrix(configuration_type& configuration,
                        dca::linalg::Matrix<Real, dca::linalg::GPU>& G0, e_spin_states_type spin);

  int deviceFingerprint() const {
    return G0_r_t_GPU.deviceFingerprint() + akima_coefficients_GPU.deviceFingerprint() +
           r1_minus_r0.deviceFingerprint();
  }

private:
  int thread_id;
  int stream_id;

  using G0_INTERPOLATION_TEMPLATE<Parameters, Real>::parameters;
  using G0_INTERPOLATION_TEMPLATE<Parameters, Real>::concurrency;

  using G0_INTERPOLATION_TEMPLATE<Parameters, Real>::G0_r_t_shifted;
  using G0_INTERPOLATION_TEMPLATE<Parameters, Real>::grad_G0_r_t_shifted;

  using G0_INTERPOLATION_TEMPLATE<Parameters, Real>::akima_coefficients;

  using G0_INTERPOLATION_TEMPLATE<Parameters, Real>::r1_minus_r0;

  dca::linalg::Matrix<Real, dca::linalg::GPU> r1_min_r0_GPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> G0_r_t_CPU;
  dca::linalg::Matrix<Real, dca::linalg::GPU> G0_r_t_GPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> grad_G0_r_t_CPU;
  dca::linalg::Matrix<Real, dca::linalg::GPU> grad_G0_r_t_GPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> akima_coefficients_CPU;
  dca::linalg::Matrix<Real, dca::linalg::GPU> akima_coefficients_GPU;

  int Nb, Nr, Nt;

  dca::linalg::Vector<int, dca::linalg::CPU> b_ind;
  dca::linalg::Vector<int, dca::linalg::CPU> r_ind;
  dca::linalg::Vector<Real, dca::linalg::CPU> tau;

  dca::linalg::Vector<int, dca::linalg::GPU> b_ind_GPU;
  dca::linalg::Vector<int, dca::linalg::GPU> r_ind_GPU;
  dca::linalg::Vector<Real, dca::linalg::GPU> tau_GPU;

  using G0_INTERPOLATION_TEMPLATE<Parameters, Real>::beta;
};

template <typename Parameters, typename Real>
G0_INTERPOLATION<dca::linalg::GPU, Parameters, Real>::G0_INTERPOLATION(int id,
                                                                       Parameters& parameters_ref)
    : G0_INTERPOLATION_TEMPLATE<Parameters, Real>(id, parameters_ref),

      thread_id(id),
      stream_id(0),

      r1_min_r0_GPU(r1_minus_r0),

      G0_r_t_CPU(std::pair<int, int>(shifted_t::dmn_size(),
                                     b::dmn_size() * b::dmn_size() * r_dmn_t::dmn_size())),
      G0_r_t_GPU(std::pair<int, int>(shifted_t::dmn_size(),
                                     b::dmn_size() * b::dmn_size() * r_dmn_t::dmn_size())),

      grad_G0_r_t_CPU(std::pair<int, int>(shifted_t::dmn_size(),
                                          b::dmn_size() * b::dmn_size() * r_dmn_t::dmn_size())),
      grad_G0_r_t_GPU(std::pair<int, int>(shifted_t::dmn_size(),
                                          b::dmn_size() * b::dmn_size() * r_dmn_t::dmn_size())),

      akima_coefficients_CPU(std::pair<int, int>(
          4 * shifted_t::dmn_size(), b::dmn_size() * b::dmn_size() * r_dmn_t::dmn_size())),
      akima_coefficients_GPU(std::pair<int, int>(
          4 * shifted_t::dmn_size(), b::dmn_size() * b::dmn_size() * r_dmn_t::dmn_size())),

      Nb(b::dmn_size()),
      Nr(r_dmn_t::dmn_size()),
      Nt(shifted_t::dmn_size()),

      b_ind("b_ind G0_INTERPOLATION<dca::linalg::GPU>", 4096),
      r_ind("r_ind G0_INTERPOLATION<dca::linalg::GPU>", 4096),
      tau("tau   G0_INTERPOLATION<dca::linalg::GPU>", 4096),

      b_ind_GPU("b_ind_GPU G0_INTERPOLATION<dca::linalg::GPU>", 4096),
      r_ind_GPU("r_ind_GPU G0_INTERPOLATION<dca::linalg::GPU>", 4096),
      tau_GPU("tau_GPU   G0_INTERPOLATION<dca::linalg::GPU>", 4096) {}

/*!
 *  \brief  Set the functions 'G0_r_t_shifted' and 'grad_G0_r_t_shifted'
 */
template <typename Parameters, typename Real>
template <class MOMS_type>
void G0_INTERPOLATION<dca::linalg::GPU, Parameters, Real>::initialize(MOMS_type& MOMS) {
  G0_INTERPOLATION_TEMPLATE<Parameters, Real>::initialize(MOMS);

  for (int t_ind = 0; t_ind < shifted_t::dmn_size(); t_ind++) {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < b::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < b::dmn_size(); nu0_ind++) {
          G0_r_t_CPU(t_ind, nu0_ind + Nb * (nu1_ind + Nb * r_ind)) =
              G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind);
          grad_G0_r_t_CPU(t_ind, nu0_ind + Nb * (nu1_ind + Nb * r_ind)) =
              grad_G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind);
        }
      }
    }
  }

  G0_r_t_GPU = G0_r_t_CPU;
  grad_G0_r_t_GPU = grad_G0_r_t_CPU;

  for (int t_ind = 0; t_ind < shifted_t::dmn_size(); t_ind++)
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
      for (int nu1_ind = 0; nu1_ind < b::dmn_size(); nu1_ind++)
        for (int nu0_ind = 0; nu0_ind < b::dmn_size(); nu0_ind++)
          for (int l_ind = 0; l_ind < 4; l_ind++)
            akima_coefficients_CPU(l_ind + 4 * t_ind, nu0_ind + Nb * (nu1_ind + Nb * r_ind)) =
                akima_coefficients(l_ind, nu0_ind, nu1_ind, r_ind, t_ind);

  akima_coefficients_GPU = akima_coefficients_CPU;
}

template <typename Parameters, typename Real>
template <class configuration_type>
void G0_INTERPOLATION<dca::linalg::GPU, Parameters, Real>::build_G0_matrix(
    configuration_type& configuration, dca::linalg::Matrix<Real, dca::linalg::GPU>& G0_e_spin,
    e_spin_states_type e_spin) {
  // profiler_t profiler(concurrency, "G0-matrix (build)", "CT-AUX", __LINE__);

  std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  G0_e_spin.resizeNoCopy(configuration_size);

  b_ind.resize(configuration_size);
  r_ind.resize(configuration_size);
  tau.resize(configuration_size);

  for (int l = 0; l < configuration_size; ++l) {
    b_ind[l] = configuration_e_spin[l].get_band();
    r_ind[l] = configuration_e_spin[l].get_r_site();
    tau[l] = configuration_e_spin[l].get_tau();
  }

  b_ind_GPU = b_ind;
  r_ind_GPU = r_ind;
  tau_GPU = tau;

  int first_shuffled_index = 0;  // configuration.get_first_shuffled_spin_index(e_spin);
  g0kernels::akima_interpolation_on_GPU(
      Nb, Nr, Nt, beta, first_shuffled_index, configuration_size, b_ind_GPU.ptr(), r_ind_GPU.ptr(),
      tau_GPU.ptr(), G0_e_spin.ptr(), G0_e_spin.size(), G0_e_spin.capacity(), r1_min_r0_GPU.ptr(),
      r1_min_r0_GPU.size(), r1_min_r0_GPU.capacity(), akima_coefficients_GPU.ptr(),
      akima_coefficients_GPU.size(), akima_coefficients_GPU.capacity());
}

template <typename Parameters, typename Real>
template <class configuration_type>
void G0_INTERPOLATION<dca::linalg::GPU, Parameters, Real>::update_G0_matrix(
    configuration_type& configuration, dca::linalg::Matrix<Real, dca::linalg::GPU>& G0_e_spin,
    e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  G0_e_spin.resize(configuration_size);

  int first_shuffled_index = configuration.get_first_shuffled_spin_index(e_spin);

  b_ind.resize(configuration_size);
  r_ind.resize(configuration_size);
  tau.resize(configuration_size);

  for (int l = 0; l < configuration_size; ++l) {
    b_ind[l] = configuration_e_spin[l].get_band();
    r_ind[l] = configuration_e_spin[l].get_r_site();
    tau[l] = configuration_e_spin[l].get_tau();
  }

  auto stream = linalg::util::getStream(thread_id, stream_id);
  b_ind_GPU.setAsync(b_ind, stream);
  r_ind_GPU.setAsync(r_ind, stream);
  tau_GPU.setAsync(tau, stream);

  g0kernels::akima_interpolation_on_GPU(
      Nb, Nr, Nt, beta, first_shuffled_index, configuration_size, b_ind_GPU.ptr(), r_ind_GPU.ptr(),
      tau_GPU.ptr(), G0_e_spin.ptr(), G0_e_spin.size(), G0_e_spin.capacity(), r1_min_r0_GPU.ptr(),
      r1_min_r0_GPU.size(), r1_min_r0_GPU.capacity(), akima_coefficients_GPU.ptr(),
      akima_coefficients_GPU.size(), akima_coefficients_GPU.capacity(), thread_id, stream_id);
}
