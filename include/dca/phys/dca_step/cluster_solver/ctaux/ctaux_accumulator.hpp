// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This class organizes the measurements in the CT-AUX QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_ACCUMULATOR_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/sp/sp_accumulator_nfft.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/accumulator_nonlocal_chi.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/accumulator_nonlocal_g.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/tp_equal_time_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/feynman_expansion_order_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ct_aux_hs_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_pair.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/mc_accumulator_data.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/quantum/numerical_error_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
class CtauxAccumulator : public MC_accumulator_data {
public:
  using this_type = CtauxAccumulator<device_t, parameters_type, Data>;
  using DataType = Data;

  typedef parameters_type my_parameters_type;

  typedef vertex_pair<parameters_type> vertex_pair_type;
  typedef vertex_singleton vertex_singleton_type;

  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;
  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  typedef RClusterDmn r_dmn_t;
  typedef KClusterDmn k_dmn_t;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef CT_AUX_HS_configuration<parameters_type> configuration_type;

  CtauxAccumulator(parameters_type& parameters_ref, Data& data_ref, int id);

  template <typename Writer>
  void write(Writer& writer);

  void initialize(int dca_iteration);

  template <typename walker_type>
  void update_from(walker_type& walker);

  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sum_to(this_type& other);

  void finalize();

  std::vector<vertex_singleton_type>& get_configuration(e_spin_states_type e_spin = e_UP);

  func::function<double, func::dmn_0<domains::numerical_error_domain>>& get_error_distribution() {
    return error;
  }

  func::function<double, func::dmn_0<Feynman_expansion_order_domain>>& get_visited_expansion_order_k() {
    return visited_expansion_order_k;
  }

  // equal time-measurements
  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& get_G_r_t() {
    return G_r_t;
  }
  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& get_G_r_t_stddev() {
    return G_r_t_stddev;
  }

  func::function<double, func::dmn_variadic<b, r_dmn_t>>& get_charge_cluster_moment() {
    return charge_cluster_moment;
  }
  func::function<double, func::dmn_variadic<b, r_dmn_t>>& get_magnetic_cluster_moment() {
    return magnetic_cluster_moment;
  }
  func::function<double, func::dmn_variadic<b, r_dmn_t>>& get_dwave_pp_correlator() {
    return dwave_pp_correlator;
  }

  // sp-measurements
  const auto& get_M_r_w() const {
    return M_r_w;
  }

  // TODO: Remove getter methods that return a non-const reference.
  auto& get_M_r_w() {
    return M_r_w;
  }
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& get_M_r_w_squared() {
    return M_r_w_squared;
  }

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& get_M_k_w() {
    return M_k_w;
  }

  // tp-measurements
  auto& get_G4() {
    return G4;
  }

#ifdef MEASURE_ERROR_BARS
  void store_standard_deviation(int nr_measurements, std::ofstream& points_file,
                                std::ofstream& norm_file);
  void update_sum_squares();
#endif

private:
  void compute_M_v_v(std::vector<vertex_singleton_type>& configuration_e_spin,
                     dca::linalg::Matrix<double, dca::linalg::CPU>& N,
                     dca::linalg::Matrix<double, dca::linalg::CPU>& M, int thread_id, int stream_id);

  void compute_M_v_v(std::vector<vertex_singleton_type>& configuration_e_spin,
                     dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                     dca::linalg::Matrix<double, dca::linalg::CPU>& M, int thread_id, int stream_id);

  void accumulate_single_particle_quantities();

  void accumulate_equal_time_quantities();

  void accumulate_two_particle_quantities();

protected:
  parameters_type& parameters;
  Data& data_;
  concurrency_type& concurrency;

  int thread_id;

  using MC_accumulator_data::GFLOP;

  using MC_accumulator_data::DCA_iteration;
  using MC_accumulator_data::number_of_measurements;

  using MC_accumulator_data::current_sign;
  using MC_accumulator_data::accumulated_sign;

  CV<parameters_type> CV_obj;

  dca::linalg::Vector<double, dca::linalg::CPU> exp_V_minus_one;

  std::vector<vertex_singleton_type> HS_configuration_e_UP;
  std::vector<vertex_singleton_type> HS_configuration_e_DN;

  dca::linalg::Matrix<double, dca::linalg::CPU> M_e_UP;
  dca::linalg::Matrix<double, dca::linalg::CPU> M_e_DN;

  func::function<double, func::dmn_0<domains::numerical_error_domain>> error;
  func::function<double, func::dmn_0<Feynman_expansion_order_domain>> visited_expansion_order_k;

  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>> G_r_t;
  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>> G_r_t_stddev;

  func::function<double, func::dmn_variadic<b, r_dmn_t>> charge_cluster_moment;
  func::function<double, func::dmn_variadic<b, r_dmn_t>> magnetic_cluster_moment;
  func::function<double, func::dmn_variadic<b, r_dmn_t>> dwave_pp_correlator;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>> M_r_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>> M_r_w_squared;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>> M_r_w_stddev;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>> M_k_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>> M_k_w_stddev;

  ctaux::SpAccumulatorNfft<parameters_type, Data> single_particle_accumulator_obj;

  ctaux::TpEqualTimeAccumulator<parameters_type, Data> MC_two_particle_equal_time_accumulator_obj;

  typename Data::TpGreensFunction G4;

  accumulator_nonlocal_G<parameters_type, Data> accumulator_nonlocal_G_obj;
  accumulator_nonlocal_chi<parameters_type, Data> accumulator_nonlocal_chi_obj;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
CtauxAccumulator<device_t, parameters_type, Data>::CtauxAccumulator(parameters_type& parameters_ref,
                                                                    Data& data_ref, int id)
    : MC_accumulator_data(),

      parameters(parameters_ref),
      data_(data_ref),
      concurrency(parameters.get_concurrency()),

      thread_id(id),

      CV_obj(parameters),

      exp_V_minus_one(64),

      //     HS_configuration(0),

      HS_configuration_e_UP(0),
      HS_configuration_e_DN(0),

      M_e_UP(0, 64),
      M_e_DN(0, 64),

      error("numerical-error-distribution-of-N-matrices"),
      visited_expansion_order_k("<k>"),

      G_r_t("G_r_t_measured"),
      G_r_t_stddev("G_r_t_stddev"),

      charge_cluster_moment("charge-cluster-moment"),
      magnetic_cluster_moment("magnetic-cluster-moment"),
      dwave_pp_correlator("dwave-pp-correlator"),

      M_r_w("M_r_w"),
      M_r_w_squared("M_r_w_squared"),
      M_r_w_stddev("M_r_w_stddev"),

      M_k_w("M_k_w"),
      M_k_w_stddev("M_k_w_stddev"),

      single_particle_accumulator_obj(parameters),

      MC_two_particle_equal_time_accumulator_obj(parameters, data_, id),

      G4("two_particle_function"),

      accumulator_nonlocal_G_obj(parameters, data_, id),
      accumulator_nonlocal_chi_obj(parameters, data_, id, G4) {}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::initialize(int dca_iteration) {
  profiler_type profiler(__FUNCTION__, "CT-AUX accumulator", __LINE__, thread_id);

  MC_accumulator_data::initialize(dca_iteration);

  CV_obj.initialize(data_);

  for (int i = 0; i < visited_expansion_order_k.size(); i++)
    visited_expansion_order_k(i) = 0;

  single_particle_accumulator_obj.initialize();

  M_r_w = 0.;
  M_r_w_squared = 0.;

  for (int i = 0; i < M_k_w.size(); i++)
    M_k_w(i) = 0;

  if (parameters.additional_time_measurements()) {
    G_r_t = 0.;
    G_r_t_stddev = 0.;

    charge_cluster_moment = 0;
    magnetic_cluster_moment = 0;
    dwave_pp_correlator = 0;

    MC_two_particle_equal_time_accumulator_obj.initialize();
  }

  if (parameters.get_four_point_type() != NONE)
    accumulator_nonlocal_chi_obj.initialize();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::finalize() {
  profiler_type profiler(__FUNCTION__, "CT-AUX accumulator", __LINE__, thread_id);

  single_particle_accumulator_obj.finalize(M_r_w, M_r_w_squared);

  {
    for (int l = 0; l < M_r_w_stddev.size(); l++)
      M_r_w_stddev(l) = std::sqrt(abs(M_r_w_squared(l)) - std::pow(abs(M_r_w(l)), 2));

    double factor = 1. / std::sqrt(parameters.get_measurements() - 1);

    M_r_w_stddev *= factor;
  }

  if (parameters.additional_time_measurements()) {
    MC_two_particle_equal_time_accumulator_obj.finalize();  // G_r_t, G_r_t_stddev);

    G_r_t = MC_two_particle_equal_time_accumulator_obj.get_G_r_t();
    G_r_t_stddev = MC_two_particle_equal_time_accumulator_obj.get_G_r_t_stddev();

    charge_cluster_moment = MC_two_particle_equal_time_accumulator_obj.get_charge_cluster_moment();
    magnetic_cluster_moment =
        MC_two_particle_equal_time_accumulator_obj.get_magnetic_cluster_moment();
    dwave_pp_correlator = MC_two_particle_equal_time_accumulator_obj.get_dwave_pp_correlator();
  }

  //       single_particle_accumulator_obj.compute_M_r_w(M_r_w);
  //       single_particle_accumulator_obj.compute_M_r_w(M_r_w);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
std::vector<vertex_singleton>& CtauxAccumulator<device_t, parameters_type, Data>::get_configuration(
    e_spin_states_type e_spin) {
  if (e_spin == e_UP)
    return HS_configuration_e_UP;
  else
    return HS_configuration_e_DN;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
template <typename Writer>
void CtauxAccumulator<device_t, parameters_type, Data>::write(Writer& writer) {
//       writer.open_group("CT-AUX-SOLVER-functions");

#ifdef DCA_WITH_QMC_BIT
  writer.execute(error);
#endif  // DCA_WITH_QMC_BIT

  writer.execute(visited_expansion_order_k);

  //       writer.execute(M_r_w);
  //       writer.execute(M_r_w_stddev);

  if (parameters.additional_time_measurements()) {
    writer.execute(charge_cluster_moment);
    writer.execute(magnetic_cluster_moment);
    writer.execute(dwave_pp_correlator);

    writer.execute(G_r_t);
    writer.execute(G_r_t_stddev);
  }

  //       writer.close_group();
}

/*!
 *  \brief Get all the information from the walker in order to start a measurement.
 *
 *   \f{eqnarray}{
 *    M_{i,j} &=& (e^{V_i}-1) N_{i,j}
 *   \f}
 */
template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
template <typename walker_type>
void CtauxAccumulator<device_t, parameters_type, Data>::update_from(walker_type& walker) {
  profiler_type profiler("update from", "CT-AUX accumulator", __LINE__, thread_id);

  GFLOP += walker.get_Gflop();

  current_sign = walker.get_sign();

  configuration_type& full_configuration = walker.get_configuration();

  {
    int k = full_configuration.get_number_of_interacting_HS_spins();

    if (k < visited_expansion_order_k.size())
      visited_expansion_order_k(k) += 1.;
  }

#ifdef DCA_WITH_QMC_BIT
  error += walker.get_error_distribution();
  walker.get_error_distribution() = 0;
#endif  // DCA_WITH_QMC_BIT

  HS_configuration_e_DN.resize(full_configuration.get(e_DN).size());
  copy(full_configuration.get(e_DN).begin(), full_configuration.get(e_DN).end(),
       HS_configuration_e_DN.begin());

  HS_configuration_e_UP.resize(full_configuration.get(e_UP).size());
  copy(full_configuration.get(e_UP).begin(), full_configuration.get(e_UP).end(),
       HS_configuration_e_UP.begin());

  compute_M_v_v(HS_configuration_e_DN, walker.get_N(e_DN), M_e_DN, walker.get_thread_id(), 0);
  compute_M_v_v(HS_configuration_e_UP, walker.get_N(e_UP), M_e_UP, walker.get_thread_id(), 0);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::measure() {
  number_of_measurements += 1;
  accumulated_sign += current_sign;

  accumulate_single_particle_quantities();

  if (DCA_iteration == parameters.get_dca_iterations() - 1 &&
      parameters.additional_time_measurements())
    accumulate_equal_time_quantities();

  if (DCA_iteration == parameters.get_dca_iterations() - 1 && parameters.get_four_point_type() != NONE)
    accumulate_two_particle_quantities();
}

#ifdef MEASURE_ERROR_BARS

/*!
 *  \brief Output and store standard deviation and error.
 *
 *  It computes and write to the given files the standard deviation of the measurements of the one
 * particle accumulator.
 *  It outputs the L1-Norm, i.e. \f$\sum_{i=1}^N \left|x_i\right|/N\f$, the L2-Norm, i.e.
 * \f$\sqrt{\sum_{i=1}^N \left|x_i\right|^2/N}\f$,
 *  and the Linf-Norm, i.e. \f$\max_{i=1}^N \left|x_i\right|\f$ of the standard deviation and of the
 * error.
 */
template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::store_standard_deviation(
    int nr_measurements, std::ofstream& points_file, std::ofstream& norm_file) {
  single_particle_accumulator_obj.store_standard_deviation(nr_measurements, points_file, norm_file);
}

/*!
 *  \brief Update the sum of the squares of the measurements of the single particle accumulator.
 *         It has to be called after each measurement.
 */
template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::update_sum_squares() {
  single_particle_accumulator_obj.update_sum_squares();
}
#endif

/*!
 *   \f{eqnarray}{
 *    M_{i,j} &=& (e^{V_i}-1) N_{i,j} \nonumber
 *   \f}
 */
template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::compute_M_v_v(
    std::vector<vertex_singleton_type>& configuration_e_spin,
    dca::linalg::Matrix<double, dca::linalg::CPU>& N,
    dca::linalg::Matrix<double, dca::linalg::CPU>& M, int /*walker_thread_id*/,
    int /*walker_stream_id*/) {
  assert(int(configuration_e_spin.size()) == N.nrRows() && N.is_square());

  // What happens if configuration_size = 0?
  int configuration_size = configuration_e_spin.size();

  M.resizeNoCopy(N.size());

  exp_V_minus_one.resize(configuration_size);

  for (int i = 0; i < configuration_size; ++i)
    exp_V_minus_one[i] = CV_obj.exp_V(configuration_e_spin[i]) - 1.;

  dca::linalg::matrixop::multiplyDiagonalLeft(exp_V_minus_one, N, M);
}

/*!
 *   \f{eqnarray}{
 *    M_{i,j} &=& (e^{V_i}-1) N_{i,j} \nonumber
 *   \f}
 */
#ifdef DCA_HAVE_CUDA
template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::compute_M_v_v(
    std::vector<vertex_singleton_type>& configuration_e_spin,
    dca::linalg::Matrix<double, dca::linalg::GPU>& N,
    dca::linalg::Matrix<double, dca::linalg::CPU>& M, int walker_thread_id, int walker_stream_id) {
  assert(int(configuration_e_spin.size()) == N.nrRows() && N.is_square());

  M.set(N, walker_thread_id, walker_stream_id);

  // What happens if configuration_size = 0?
  int configuration_size = configuration_e_spin.size();
  exp_V_minus_one.resize(configuration_size);

  for (int i = 0; i < configuration_size; ++i)
    exp_V_minus_one[i] = CV_obj.exp_V(configuration_e_spin[i]) - 1.;

  dca::linalg::matrixop::multiplyDiagonalLeft(exp_V_minus_one, M, M);
}
#endif  // DCA_HAVE_CUDA

/*************************************************************
 **                                                         **
 **                    G2 - MEASUREMENTS                    **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::accumulate_single_particle_quantities() {
  profiler_type profiler("sp-accumulation", "CT-AUX accumulator", __LINE__, thread_id);

  single_particle_accumulator_obj.accumulate_M_r_w(HS_configuration_e_DN, M_e_DN, current_sign, e_DN);
  single_particle_accumulator_obj.accumulate_M_r_w(HS_configuration_e_UP, M_e_UP, current_sign, e_UP);

  GFLOP += 2. * 8. * M_e_DN.size().first * M_e_DN.size().first * (1.e-9);
  GFLOP += 2. * 8. * M_e_UP.size().first * M_e_UP.size().first * (1.e-9);
}

/*************************************************************
 **                                                         **
 **                 equal-time - MEASUREMENTS               **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::accumulate_equal_time_quantities() {
  profiler_type profiler("equal-time-measurements", "CT-AUX accumulator", __LINE__, thread_id);

  MC_two_particle_equal_time_accumulator_obj.compute_G_r_t(HS_configuration_e_DN, M_e_DN,
                                                           HS_configuration_e_UP, M_e_UP);

  MC_two_particle_equal_time_accumulator_obj.accumulate_G_r_t(current_sign);

  MC_two_particle_equal_time_accumulator_obj.accumulate_moments(current_sign);

  MC_two_particle_equal_time_accumulator_obj.accumulate_dwave_pp_correlator(current_sign);

  GFLOP += MC_two_particle_equal_time_accumulator_obj.get_GFLOP();
}

/*************************************************************
 **                                                         **
 **                 nonlocal \chi - MEASUREMENTS            **
 **                                                         **
 *************************************************************/

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::accumulate_two_particle_quantities() {
  {
    profiler_type profiler("tp-accumulation nonlocal G", "CT-AUX accumulator", __LINE__, thread_id);

    accumulator_nonlocal_G_obj.execute(HS_configuration_e_UP, M_e_UP, HS_configuration_e_DN, M_e_DN);
  }

  GFLOP += accumulator_nonlocal_G_obj.get_GFLOP();

  {
    profiler_type profiler("tp-accumulation nonlocal chi", "CT-AUX accumulator", __LINE__, thread_id);
    accumulator_nonlocal_chi_obj.execute(current_sign, accumulator_nonlocal_G_obj);
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class Data>
void CtauxAccumulator<device_t, parameters_type, Data>::sum_to(this_type& other) {
  finalize();

  other.get_Gflop() += get_Gflop();

  other.get_sign() += get_sign();
  other.get_number_of_measurements() += get_number_of_measurements();

  {
    for (int i = 0; i < visited_expansion_order_k.size(); i++)
      other.get_visited_expansion_order_k()(i) += visited_expansion_order_k(i);

    for (int i = 0; i < error.size(); i++)
      other.get_error_distribution()(i) += error(i);
  }

  {  // equal time measurements
    for (int i = 0; i < G_r_t.size(); i++)
      other.get_G_r_t()(i) += G_r_t(i);

    for (int i = 0; i < G_r_t_stddev.size(); i++)
      other.get_G_r_t_stddev()(i) += G_r_t_stddev(i);

    for (int i = 0; i < charge_cluster_moment.size(); i++)
      other.get_charge_cluster_moment()(i) += charge_cluster_moment(i);

    for (int i = 0; i < magnetic_cluster_moment.size(); i++)
      other.get_magnetic_cluster_moment()(i) += magnetic_cluster_moment(i);

    for (int i = 0; i < dwave_pp_correlator.size(); i++)
      other.get_dwave_pp_correlator()(i) += dwave_pp_correlator(i);
  }

  {  // sp-measurements
    for (int i = 0; i < M_r_w.size(); i++)
      other.get_M_r_w()(i) += M_r_w(i);

    for (int i = 0; i < M_r_w_squared.size(); i++)
      other.get_M_r_w_squared()(i) += M_r_w_squared(i);
  }

  {  // tp-measurements
    for (int i = 0; i < G4.size(); i++)
      other.get_G4()(i) += G4(i);
  }
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_ACCUMULATOR_HPP
