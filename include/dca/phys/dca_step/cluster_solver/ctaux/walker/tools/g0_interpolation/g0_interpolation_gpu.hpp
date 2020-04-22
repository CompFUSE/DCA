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
// Template specialization for GPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_GPU_HPP

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_tmpl.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_base.hpp"

#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/multi_vector.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_kernels.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <typename Parameters, typename Real>
class G0Interpolation<dca::linalg::GPU, Parameters, Real>
    : public G0InterpolationBase<Parameters, Real> {
public:
  using vertex_singleton_type = vertex_singleton;
  using shifted_t = func::dmn_0<domains::time_domain_left_oriented>;
  using b = func::dmn_0<domains::electron_band_domain>;

  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using r_dmn_t = RClusterDmn;

  typedef typename Parameters::concurrency_type concurrency_type;
  typedef typename Parameters::profiler_type profiler_t;

  using Base = G0InterpolationBase<Parameters, Real>;

public:
  G0Interpolation(int id, const Parameters& parameters);

  template <class MOMS_type>
  void initialize(MOMS_type& MOMS);

  template <class Configuration>
  void build_G0_matrix(Configuration& configuration,
                       dca::linalg::Matrix<Real, dca::linalg::GPU>& G0, e_spin_states_type spin);

  template <class Configuration>
  void update_G0_matrix(Configuration& configuration,
                        dca::linalg::Matrix<Real, dca::linalg::GPU>& G0, e_spin_states_type spin);

  int deviceFingerprint() const {
    return G0_r_t_GPU.deviceFingerprint() + akima_coefficients_GPU.deviceFingerprint() +
           r1_minus_r0.deviceFingerprint();
  }

private:
  template <class Configuration>
  void uploadConfiguration(const Configuration& configuration);

  int thread_id;
  int stream_id;

  using Base::parameters;
  using Base::concurrency;

  using Base::G0_r_t_shifted;
  using Base::grad_G0_r_t_shifted;

  using Base::akima_coefficients;

  using Base::r1_minus_r0;

  dca::linalg::Matrix<Real, dca::linalg::GPU> r1_min_r0_GPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> G0_r_t_CPU;
  dca::linalg::Matrix<Real, dca::linalg::GPU> G0_r_t_GPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> grad_G0_r_t_CPU;
  dca::linalg::Matrix<Real, dca::linalg::GPU> grad_G0_r_t_GPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> akima_coefficients_CPU;
  dca::linalg::Matrix<Real, dca::linalg::GPU> akima_coefficients_GPU;

  int Nb, Nr, Nt;

  // Store indices of band and site, and the imaginary time.
  linalg::MultiVector<linalg::CPU, int, int, Real> g0_labels_cpu_;
  linalg::MultiVector<linalg::GPU, int, int, Real> g0_labels_gpu_;

  using Base::beta;

  linalg::util::CudaEvent config_copied_;
};

template <typename Parameters, typename Real>
G0Interpolation<dca::linalg::GPU, Parameters, Real>::G0Interpolation(int id,
                                                                     const Parameters& parameters_ref)
    : Base(id, parameters_ref),

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
      Nt(shifted_t::dmn_size()) {}

/*!
 *  \brief  Set the functions 'G0_r_t_shifted' and 'grad_G0_r_t_shifted'
 */
template <typename Parameters, typename Real>
template <class MOMS_type>
void G0Interpolation<dca::linalg::GPU, Parameters, Real>::initialize(MOMS_type& MOMS) {
  Base::initialize(MOMS);

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
template <class Configuration>
void G0Interpolation<dca::linalg::GPU, Parameters, Real>::build_G0_matrix(
    Configuration& configuration, dca::linalg::Matrix<Real, dca::linalg::GPU>& G0_e_spin,
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

  uploadConfiguration(configuration_e_spin);

  const auto b_ind_gpu = g0_labels_gpu_.template get<0>();
  const auto r_ind_gpu = g0_labels_gpu_.template get<1>();
  const auto tau_gpu = g0_labels_gpu_.template get<2>();

  int first_shuffled_index = 0;  // configuration.get_first_shuffled_spin_index(e_spin);
  g0kernels::akima_interpolation_on_GPU(
      Nb, Nr, Nt, beta, first_shuffled_index, configuration_size, b_ind_gpu, r_ind_gpu, tau_gpu,
      G0_e_spin.ptr(), G0_e_spin.size(), G0_e_spin.capacity(), r1_min_r0_GPU.ptr(),
      r1_min_r0_GPU.size(), r1_min_r0_GPU.capacity(), akima_coefficients_GPU.ptr(),
      akima_coefficients_GPU.size(), akima_coefficients_GPU.capacity(), thread_id, stream_id);
}

template <typename Parameters, typename Real>
template <class Configuration>
void G0Interpolation<dca::linalg::GPU, Parameters, Real>::update_G0_matrix(
    Configuration& configuration, dca::linalg::Matrix<Real, dca::linalg::GPU>& G0_e_spin,
    e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  G0_e_spin.resize(configuration_size);

  uploadConfiguration(configuration_e_spin);

  int first_shuffled_index = configuration.get_first_shuffled_spin_index(e_spin);

  const auto b_ind_gpu = g0_labels_gpu_.template get<0>();
  const auto r_ind_gpu = g0_labels_gpu_.template get<1>();
  const auto tau_gpu = g0_labels_gpu_.template get<2>();

  g0kernels::akima_interpolation_on_GPU(
      Nb, Nr, Nt, beta, first_shuffled_index, configuration_size, b_ind_gpu, r_ind_gpu, tau_gpu,
      G0_e_spin.ptr(), G0_e_spin.size(), G0_e_spin.capacity(), r1_min_r0_GPU.ptr(),
      r1_min_r0_GPU.size(), r1_min_r0_GPU.capacity(), akima_coefficients_GPU.ptr(),
      akima_coefficients_GPU.size(), akima_coefficients_GPU.capacity(), thread_id, stream_id);
}

template <typename Parameters, typename Real>
template <class Configuration>
void G0Interpolation<dca::linalg::GPU, Parameters, Real>::uploadConfiguration(
    const Configuration& configuration) {
  const int configuration_size = configuration.size();

  config_copied_.block();
  g0_labels_cpu_.resizeNoCopy(configuration_size);
  g0_labels_gpu_.resizeNoCopy(configuration_size);

  auto band = g0_labels_cpu_.template get<0>();
  auto site = g0_labels_cpu_.template get<1>();
  auto tau = g0_labels_cpu_.template get<2>();

  for (int l = 0; l < configuration_size; ++l) {
    band[l] = configuration[l].get_band();
    site[l] = configuration[l].get_r_site();
    tau[l] = configuration[l].get_tau();
  }

  const auto& stream = linalg::util::getStream(thread_id, stream_id);
  g0_labels_gpu_.setAsync(g0_labels_cpu_, stream);
  config_copied_.record(stream);
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_GPU_HPP
