// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
// Modified by Andrei Plamada (plamada@itp.phys.ethz.ch).
//
// This class organizes the MC walker in the SS CT-HYB QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_WALKER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_WALKER_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_typedefs.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_hybridization_solver_routines.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/walker_tools/anti_segment_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/walker_tools/full_line_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/walker_tools/segment_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/walker_tools/shift_segment_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/walker_tools/ss_hybridization_walker_routines.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/walker_tools/swap_segment_tools.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/util/plot.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class SsCtHybWalker {
public:
  typedef typename parameters_type::random_number_generator rng_type;

  typedef typename SsCtHybTypedefs<parameters_type, MOMS_type>::profiler_type profiler_type;
  typedef typename SsCtHybTypedefs<parameters_type, MOMS_type>::concurrency_type concurrency_type;

  typedef typename SsCtHybTypedefs<parameters_type, MOMS_type>::vertex_vertex_matrix_type
      vertex_vertex_matrix_type;
  typedef typename SsCtHybTypedefs<parameters_type, MOMS_type>::configuration_type configuration_type;

  using t = func::dmn_0<domains::time_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using nu_nu_r_DCA_t = func::dmn_variadic<nu, nu, RClusterDmn, t>;

  typedef func::function<vertex_vertex_matrix_type, nu> M_matrix_type;

  typedef ss_hybridization_solver_routines<parameters_type, MOMS_type> ss_hybridization_solver_routines_type;
  typedef ss_hybridization_walker_routines<parameters_type, MOMS_type, configuration_type, rng_type>
      ss_hybridization_walker_routines_type;

  typedef full_line_tools<ss_hybridization_walker_routines_type> full_line_tools_t;
  typedef anti_segment_tools<ss_hybridization_walker_routines_type> anti_segment_tools_t;
  typedef segment_tools<ss_hybridization_walker_routines_type> segment_tools_t;
  typedef shift_segment_tools<ss_hybridization_walker_routines_type> shift_segment_tools_t;
  typedef swap_segment_tools<ss_hybridization_walker_routines_type> swap_segment_tools_t;

public:
  SsCtHybWalker(parameters_type& parameters_ref, MOMS_type& MOMS_ref, rng_type& rng_ref, int id = 0);

  /*!
   *  \brief Initializes the configuration and sets \f$\mu_i = \frac12 \sum_j
   * \frac{U_{ij}+U_{ji}}{2}\f$.
   */
  void initialize();  // func::function<double, nu> mu_DC);

  /*!
   *  \brief Returns if the configuration has gone through a warm-up sweep.
   */
  bool& is_thermalized() {
    return thermalized;
  }

  /*!
   *  \brief Goes through an integration sweep. Do N insertion, removal, shift or swap steps.
   *  The insertion and removal steps include the insertion or removal of a full line, a segment or
   * an anti-segment.
   */
  void doSweep();

  /*!
   *  \brief Goes through an integration step. Do N insertion, removal, shift or swap steps.
   *  The insertion and removal steps include the insertion or removal of a full line, a segment or
   * an anti-segment.
   */
  void do_step();

  /*!
   *  \brief Returns the QMC sign.
   */
  double get_sign() {
    return sign;
  }

  /*!
   *  \brief Returns the current configuration.
   */
  configuration_type& get_configuration() {
    return configuration;
  }

  /*!
   *  \brief Returns the current inverse hybridization matrix \f$M = F_r(t)^{-1}\f$.
   */
  M_matrix_type& get_M_matrices() {
    return M;
  }

  /*!
   *  \brief Returns the hybridization_tools object
   */
  ss_hybridization_walker_routines_type& get_ss_hybridization_walker_routines() {
    return ss_hybridization_walker_routines_obj;
  }

  // Writes the current progress and the configuration size to stdout.
  // TODO: Before this method can be made const, SS_CT_HYB_configuration needs to be made const
  //       correct.
  void updateShell(const int done, const int total) /*const*/;

  // Writes a summary of the walker's Markov chain updates to stdout.
  void printSummary() const;

private:
  void test_interpolation();

  int get_random_interacting_flavor();

  void do_insert_remove(int j);
  void insert_or_remove_full_line(int j);
  void insert_or_remove_anti_segment(int j);
  void insert_or_remove_segment(int j);
  void shift_segment(int j);
  void swap_random_orbitals();

private:
  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

  rng_type& rng;

  int thread_id;

  configuration_type configuration;

  ss_hybridization_solver_routines_type ss_hybridization_solver_routines_obj;
  ss_hybridization_walker_routines_type ss_hybridization_walker_routines_obj;

  func::function<double, nu>& mu;
  func::function<double, nu_nu_r_DCA_t>& F_r_t;

  full_line_tools_t full_line_tools_obj;
  anti_segment_tools_t anti_segment_tools_obj;
  segment_tools_t segment_tools_obj;
  shift_segment_tools_t shift_segment_tools_obj;
  swap_segment_tools_t swap_segment_tools_obj;

  func::function<vertex_vertex_matrix_type, nu> M;

  bool thermalized;

  double sign;

  size_t nb_updates;
  size_t nb_successfull_updates;

  double percentage_steps;
  double percentage_shifts;
  double percentage_swaps;

  double total;

  double p_0;
  double p_1;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
SsCtHybWalker<device_t, parameters_type, MOMS_type>::SsCtHybWalker(parameters_type& parameters_ref,
                                                                   MOMS_type& MOMS_ref,
                                                                   rng_type& rng_ref, int id)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      rng(rng_ref),

      thread_id(id),

      configuration(),

      ss_hybridization_solver_routines_obj(parameters, MOMS),
      ss_hybridization_walker_routines_obj(parameters, MOMS, configuration, rng),

      mu(ss_hybridization_solver_routines_obj.get_mu()),
      F_r_t(ss_hybridization_solver_routines_obj.get_F_r_t()),

      full_line_tools_obj(ss_hybridization_walker_routines_obj),
      anti_segment_tools_obj(ss_hybridization_walker_routines_obj),
      segment_tools_obj(ss_hybridization_walker_routines_obj),
      shift_segment_tools_obj(ss_hybridization_walker_routines_obj),
      swap_segment_tools_obj(ss_hybridization_walker_routines_obj),

      M("M-matrices"),

      thermalized(false),

      sign(1) {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::printSummary() const {
  std::cout.unsetf(std::ios_base::floatfield);
  std::cout << "\n"
            << "Walker: process ID = " << concurrency.id() << ", thread ID = " << thread_id << "\n"
            << "-------------------------------------------\n"
            << "nb_successfull_updates/nb_updates : "
            << static_cast<double>(nb_successfull_updates) / static_cast<double>(nb_updates) << "\n"
            << std::endl;

  std::cout << std::scientific;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::initialize() {
  ss_hybridization_solver_routines_obj.initialize_functions();

  ss_hybridization_walker_routines_obj.initialize_akima_coefficients(F_r_t);

  // test_interpolation();

  {
    sign = 1;

    nb_updates = 0;
    nb_successfull_updates = 0;
  }

  {
    percentage_steps = parameters.get_steps_per_sweep();
    percentage_shifts = parameters.get_shifts_per_sweep();

    total = (percentage_steps + percentage_shifts);

    p_0 = (percentage_steps + 0) / total;
    p_1 = (percentage_steps + percentage_shifts) / total;
  }

  {
    configuration.initialize();

    is_thermalized() = false;

    for (int i = 0; i < M.size(); i++) {
      M(i).resize(0);
    }
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::test_interpolation() {
  std::cout << __FUNCTION__ << std::endl;

  util::Plot::plotBandsLinesPoints(F_r_t);

  {
    double beta = parameters.get_beta();

    int coor[2];

    for (int b_i = 0; b_i < b::dmn_size(); b_i++) {
      for (int s_i = 0; s_i < s::dmn_size(); s_i++) {
        coor[0] = b_i;
        coor[1] = s_i;

        std::vector<double> x(0);
        std::vector<double> y(0);

        for (int t_i = 0; t_i < 10 * t::dmn_size(); t_i++) {
          double t_val = -beta + 2 * beta / (10 * t::dmn_size() - 1) * t_i;

          if (t_i == 0)
            t_val += 1.e-6;

          if (t_i == 10 * t::dmn_size() - 1)
            t_val -= 1.e-6;

          double F_val = ss_hybridization_walker_routines_obj.interpolate_F(coor, t_val, F_r_t);

          x.push_back(t_val);
          y.push_back(F_val);
        }

        util::Plot::plotPoints(x, y);
      }
    }
  }

  throw std::logic_error(__FUNCTION__);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::doSweep() {
  double factor = 1.;
  if (thermalized)
    factor = parameters.get_sweeps_per_measurement();

  int nr_of_segments = std::max(16, configuration.size());

  double ratio = double(nb_successfull_updates + 1) / double(nb_updates + 1);
  int factor2 = std::max(1, int(std::ceil(1.0 / ratio)));

  int nb_steps = nr_of_segments * factor * factor2;

  for (int l = 0; l < nb_steps; l++)
    do_step();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
int SsCtHybWalker<device_t, parameters_type, MOMS_type>::get_random_interacting_flavor() {
  int spin = s::dmn_size() * rng();
  int int_band = parameters.get_interacting_orbitals().size() * rng();

  return parameters.get_interacting_orbitals()[int_band] + spin * b::dmn_size();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::do_step() {
  double p = rng();

  int so_ind = get_random_interacting_flavor();

  if (p < p_0) {
    do_insert_remove(so_ind);
  }

  if (p_0 < p and p < p_1) {
    shift_segment(so_ind);
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::do_insert_remove(int so_ind) {
  double rn = rng();

  if (configuration.get_vertices(so_ind).size() == 0) {
    bool succes;
    nb_updates += 1;

    if (configuration.get_full_line(so_ind)) {
      if (rn < 0.75)
        succes = full_line_tools_obj.insert_or_remove(so_ind, mu(so_ind));
      else
        succes = anti_segment_tools_obj.insert_anti_segment(so_ind, mu(so_ind), sign, M, F_r_t);
    }
    else {
      if (rn < 0.75)
        succes = full_line_tools_obj.insert_or_remove(so_ind, mu(so_ind));
      else
        succes = segment_tools_obj.insert_segment(so_ind, mu(so_ind), sign, M, F_r_t);
    }
    nb_successfull_updates += succes ? 1 : 0;
  }
  else {
    if (rn < 0.5) {
      insert_or_remove_anti_segment(so_ind);
    }
    else {
      insert_or_remove_segment(so_ind);
    }
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::insert_or_remove_full_line(int j) {
  nb_updates += 1;

  bool succes = full_line_tools_obj.insert_or_remove(j, mu(j));

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::insert_or_remove_anti_segment(int j) {
  nb_updates += 1;

  bool succes;
  if (rng() < 0.5)
    succes = anti_segment_tools_obj.insert_anti_segment(j, mu(j), sign, M, F_r_t);
  else
    succes = anti_segment_tools_obj.remove_anti_segment(j, mu(j), sign, M, F_r_t);

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::insert_or_remove_segment(int j) {
  nb_updates += 1;

  bool succes;
  if (rng() < 0.5)
    succes = segment_tools_obj.insert_segment(j, mu(j), sign, M, F_r_t);
  else
    succes = segment_tools_obj.remove_segment(j, mu(j), sign, M, F_r_t);

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::shift_segment(int j) {
  nb_updates += 1;

  bool succes;
  if (rng() < 0.5)
    succes = shift_segment_tools_obj.shift_segment_start_vertex(j, mu(j), sign, M, F_r_t);
  else
    succes = shift_segment_tools_obj.shift_segment_end_vertex(j, mu(j), sign, M, F_r_t);

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::swap_random_orbitals() {
  nb_updates += 1;

  int i = get_random_interacting_flavor();
  int j = get_random_interacting_flavor();

  bool succes = swap_segment_tools_obj.swap_orbitals(i, j, mu, sign, M, F_r_t);
  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void SsCtHybWalker<device_t, parameters_type, MOMS_type>::updateShell(const int done,
                                                                       const int total) {
  if (concurrency.id() == concurrency.first() && total > 10 && (done % (total / 10)) == 0) {
    std::cout.unsetf(std::ios_base::floatfield);

    std::cout << "\t\t\t" << std::setw(14)
              << static_cast<double>(done) / static_cast<double>(total) * 100. << " % completed"
              << "\t" << std::setw(11) << "<k> = " << configuration.size() << "\t"
              << dca::util::print_time() << std::endl;

    std::cout << std::scientific;
  }
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_WALKER_HPP
