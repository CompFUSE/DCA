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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G0_MATRIX_ROUTINES_CTAUX_G0_MATRIX_ROUTINES_GPU_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G0_MATRIX_ROUTINES_CTAUX_G0_MATRIX_ROUTINES_GPU_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/ctaux_G0_matrix_routines_TEM.h"

#include <cassert>
#include <utility>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/linalg/matrix.hpp"

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_vertex_singleton.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/G0_interpolation_template.hpp"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/time_and_frequency/time_domain_left_oriented.h"

namespace DCA {
namespace QMCI {
namespace G0_INTERPOLATION_KERNELS {
// DCA::QMCI::G0_INTERPOLATION_KERNELS::

void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, double beta, int Nc, int Nv, int* b, int* r,
                                double* t, double* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, double* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                double* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs);

void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, double beta, int Nc, int Nv, int* b, int* r,
                                double* t, double* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, double* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                double* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs, int thread_id, int stream_id);

}  // G0_INTERPOLATION_KERNELS

template <typename parameters_type>
class G0_INTERPOLATION<dca::linalg::GPU, parameters_type>
    : public G0_INTERPOLATION_TEMPLATE<parameters_type> {
public:
  using vertex_singleton_type = vertex_singleton;
  using shifted_t = func::dmn_0<time_domain_left_oriented>;
  using b = func::dmn_0<electron_band_domain>;

  using r_DCA = func::dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                           CLUSTER, REAL_SPACE, BRILLOUIN_ZONE>>;
  using r_dmn_t = r_DCA;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  G0_INTERPOLATION(int id, parameters_type& parameters);

  template <class MOMS_type>
  void initialize(MOMS_type& MOMS);

  template <class configuration_type>
  void build_G0_matrix(configuration_type& configuration,
                       dca::linalg::Matrix<double, dca::linalg::GPU>& G0, e_spin_states_type spin);

  template <class configuration_type>
  void update_G0_matrix(configuration_type& configuration,
                        dca::linalg::Matrix<double, dca::linalg::GPU>& G0, e_spin_states_type spin);

private:
  void build_G0_matrix(std::vector<vertex_singleton_type>& configuration,
                       dca::linalg::Matrix<double, dca::linalg::GPU>& G0);

  double interpolate(int nu_0, int nu_1, int delta_r, double delta_time);

private:
  int thread_id;
  int stream_id;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::parameters;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::concurrency;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::G0_r_t_shifted;
  using G0_INTERPOLATION_TEMPLATE<parameters_type>::grad_G0_r_t_shifted;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::akima_coefficients;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::r1_minus_r0;

  dca::linalg::Matrix<double, dca::linalg::GPU> r1_min_r0_GPU;

  dca::linalg::Matrix<double, dca::linalg::CPU> G0_r_t_CPU;
  dca::linalg::Matrix<double, dca::linalg::GPU> G0_r_t_GPU;

  dca::linalg::Matrix<double, dca::linalg::CPU> grad_G0_r_t_CPU;
  dca::linalg::Matrix<double, dca::linalg::GPU> grad_G0_r_t_GPU;

  dca::linalg::Matrix<double, dca::linalg::CPU> akima_coefficients_CPU;
  dca::linalg::Matrix<double, dca::linalg::GPU> akima_coefficients_GPU;

  int Nb, Nr, Nt;

  dca::linalg::Vector<int, dca::linalg::CPU> b_ind;
  dca::linalg::Vector<int, dca::linalg::CPU> r_ind;
  dca::linalg::Vector<double, dca::linalg::CPU> tau;

  dca::linalg::Vector<int, dca::linalg::GPU> b_ind_GPU;
  dca::linalg::Vector<int, dca::linalg::GPU> r_ind_GPU;
  dca::linalg::Vector<double, dca::linalg::GPU> tau_GPU;

  using G0_INTERPOLATION_TEMPLATE<parameters_type>::beta;
};

template <typename parameters_type>
G0_INTERPOLATION<dca::linalg::GPU, parameters_type>::G0_INTERPOLATION(int id,
                                                                      parameters_type& parameters_ref)
    : G0_INTERPOLATION_TEMPLATE<parameters_type>(id, parameters_ref),

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
      tau_GPU("tau_GPU   G0_INTERPOLATION<dca::linalg::GPU>", 4096) {
  b_ind_GPU.setThreadAndStreamId(thread_id, stream_id);
  r_ind_GPU.setThreadAndStreamId(thread_id, stream_id);
  tau_GPU.setThreadAndStreamId(thread_id, stream_id);
}

/*!
 *  \brief  Set the functions 'G0_r_t_shifted' and 'grad_G0_r_t_shifted'
 */
template <typename parameters_type>
template <class MOMS_type>
void G0_INTERPOLATION<dca::linalg::GPU, parameters_type>::initialize(MOMS_type& MOMS) {
  G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize(MOMS);

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

  LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(G0_r_t_CPU, G0_r_t_GPU);
  LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(grad_G0_r_t_CPU, grad_G0_r_t_GPU);

  for (int t_ind = 0; t_ind < shifted_t::dmn_size(); t_ind++)
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
      for (int nu1_ind = 0; nu1_ind < b::dmn_size(); nu1_ind++)
        for (int nu0_ind = 0; nu0_ind < b::dmn_size(); nu0_ind++)
          for (int l_ind = 0; l_ind < 4; l_ind++)
            akima_coefficients_CPU(l_ind + 4 * t_ind, nu0_ind + Nb * (nu1_ind + Nb * r_ind)) =
                akima_coefficients(l_ind, nu0_ind, nu1_ind, r_ind, t_ind);

  LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(akima_coefficients_CPU,
                                                                  akima_coefficients_GPU);
}

template <typename parameters_type>
template <class configuration_type>
void G0_INTERPOLATION<dca::linalg::GPU, parameters_type>::build_G0_matrix(
    configuration_type& configuration, dca::linalg::Matrix<double, dca::linalg::GPU>& G0_e_spin,
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

  b_ind_GPU.resize(configuration_size);
  r_ind_GPU.resize(configuration_size);
  tau_GPU.resize(configuration_size);

  for (int l = 0; l < configuration_size; ++l) {
    b_ind[l] = configuration_e_spin[l].get_band();
    r_ind[l] = configuration_e_spin[l].get_r_site();
    tau[l] = configuration_e_spin[l].get_tau();
  }

  LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(
      b_ind.ptr(), b_ind_GPU.ptr(), configuration_size);  //, LIN_ALG::ASYNCHRONOUS);
  LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(
      r_ind.ptr(), r_ind_GPU.ptr(), configuration_size);  //, LIN_ALG::ASYNCHRONOUS);
  LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(
      tau.ptr(), tau_GPU.ptr(), configuration_size);  //, LIN_ALG:: SYNCHRONOUS);

  int first_shuffled_index = 0;  // configuration.get_first_shuffled_spin_index(e_spin);
  G0_INTERPOLATION_KERNELS::akima_interpolation_on_GPU(
      Nb, Nr, Nt, beta, first_shuffled_index, configuration_size, b_ind_GPU.ptr(), r_ind_GPU.ptr(),
      tau_GPU.ptr(), G0_e_spin.ptr(), G0_e_spin.size(), G0_e_spin.capacity(), r1_min_r0_GPU.ptr(),
      r1_min_r0_GPU.size(), r1_min_r0_GPU.capacity(), akima_coefficients_GPU.ptr(),
      akima_coefficients_GPU.size(), akima_coefficients_GPU.capacity());
}

template <typename parameters_type>
template <class configuration_type>
void G0_INTERPOLATION<dca::linalg::GPU, parameters_type>::update_G0_matrix(
    configuration_type& configuration, dca::linalg::Matrix<double, dca::linalg::GPU>& G0_e_spin,
    e_spin_states_type e_spin) {
  assert(G0_e_spin.get_thread_id() == thread_id);

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

  /*
    b_ind_GPU.resize(configuration_size);
    r_ind_GPU.resize(configuration_size);
    tau_GPU  .resize(configuration_size);
  */
  /*
    b_ind_GPU.reserve(configuration_size);
    r_ind_GPU.reserve(configuration_size);
    tau_GPU  .reserve(configuration_size);
  */

  for (int l = 0; l < configuration_size; ++l) {
    b_ind[l] = configuration_e_spin[l].get_band();
    r_ind[l] = configuration_e_spin[l].get_r_site();
    tau[l] = configuration_e_spin[l].get_tau();
  }

  /*
    LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(b_ind.ptr(), b_ind_GPU.ptr(),
    configuration_size, LIN_ALG::ASYNCHRONOUS);
    LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(r_ind.ptr(), r_ind_GPU.ptr(),
    configuration_size, LIN_ALG::ASYNCHRONOUS);
    LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(tau  .ptr(), tau_GPU  .ptr(),
    configuration_size, LIN_ALG:: SYNCHRONOUS);
  */

  /*
    b_ind_GPU.reserve(configuration_size);
    r_ind_GPU.reserve(configuration_size);
    tau_GPU  .reserve(configuration_size);

    LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(b_ind.ptr(), b_ind_GPU.ptr(),
    configuration_size, thread_id, stream_id);
    LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(r_ind.ptr(), r_ind_GPU.ptr(),
    configuration_size, thread_id, stream_id);
    LIN_ALG::COPY_FROM<dca::linalg::CPU, dca::linalg::GPU>::execute(tau  .ptr(), tau_GPU  .ptr(),
    configuration_size, thread_id, stream_id);
  */

  b_ind_GPU.set(b_ind, LIN_ALG::ASYNCHRONOUS);
  r_ind_GPU.set(r_ind, LIN_ALG::ASYNCHRONOUS);
  tau_GPU.set(tau, LIN_ALG::ASYNCHRONOUS);

  G0_INTERPOLATION_KERNELS::akima_interpolation_on_GPU(
      Nb, Nr, Nt, beta, first_shuffled_index, configuration_size, b_ind_GPU.ptr(), r_ind_GPU.ptr(),
      tau_GPU.ptr(), G0_e_spin.ptr(), G0_e_spin.size(), G0_e_spin.capacity(), r1_min_r0_GPU.ptr(),
      r1_min_r0_GPU.size(), r1_min_r0_GPU.capacity(), akima_coefficients_GPU.ptr(),
      akima_coefficients_GPU.size(), akima_coefficients_GPU.capacity(), thread_id, stream_id);
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G0_MATRIX_ROUTINES_CTAUX_G0_MATRIX_ROUTINES_GPU_H
