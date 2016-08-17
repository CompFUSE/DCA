// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Bart Ydens
//         Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// This class organizes the MC walker in the SS CT-HYB QMC.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_WALKER_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_WALKER_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/mc_walker.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <utility>
#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg_device_types.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_solver_routines.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_type_definitions.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_walker_tools/ANTI_SEGMENT_TOOLS.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_walker_tools/FULL_LINE_TOOLS.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_walker_tools/SEGMENT_TOOLS.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_walker_tools/SHIFT_SEGMENT_TOOLS.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_walker_tools/SWAP_SEGMENT_TOOLS.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_ss_hybridization/ss_hybridization_walker_tools/ss_hybridization_walker_routines.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type> {
public:
  typedef typename parameters_type::random_number_generator rng_type;

  typedef
      typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::profiler_type profiler_type;
  typedef
      typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::concurrency_type concurrency_type;

  typedef
      typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::vertex_vertex_matrix_type
          vertex_vertex_matrix_type;
  typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::configuration_type
      configuration_type;

  using t = dmn_0<time_domain>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using r_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     REAL_SPACE, BRILLOUIN_ZONE>>;
  using nu_nu_r_DCA_t = dmn_variadic<nu, nu, r_DCA, t>;

  typedef FUNC_LIB::function<vertex_vertex_matrix_type, nu> M_matrix_type;

  typedef ss_hybridization_solver_routines<parameters_type, MOMS_type> ss_hybridization_solver_routines_type;
  typedef ss_hybridization_walker_routines<parameters_type, MOMS_type, configuration_type, rng_type>
      ss_hybridization_walker_routines_type;

  typedef full_line_tools<ss_hybridization_walker_routines_type> full_line_tools_t;
  typedef anti_segment_tools<ss_hybridization_walker_routines_type> anti_segment_tools_t;
  typedef segment_tools<ss_hybridization_walker_routines_type> segment_tools_t;
  typedef shift_segment_tools<ss_hybridization_walker_routines_type> shift_segment_tools_t;
  typedef swap_segment_tools<ss_hybridization_walker_routines_type> swap_segment_tools_t;

public:
  MC_walker(parameters_type& parameters_ref, MOMS_type& MOMS_ref, rng_type& rng_ref, int id = 0);

  ~MC_walker();

  /*!
   *  \brief Initializes the configuration and sets \f$\mu_i = \frac12 \sum_j
   * \frac{U_{ij}+U_{ji}}{2}\f$.
   */
  void initialize();  // FUNC_LIB::function<double, nu> mu_DC);

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
  void do_sweep();

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

  /*!
   *  \brief Print the hybridization functions \f$F_k(w)\f$, \f$F_k(t)\f$ and \f$F_r(t)\f$.
   */
  template <class stream_type>
  void to_JSON(stream_type& ss);

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

  FUNC_LIB::function<double, nu>& mu;
  FUNC_LIB::function<double, nu_nu_r_DCA_t>& F_r_t;

  full_line_tools_t full_line_tools_obj;
  anti_segment_tools_t anti_segment_tools_obj;
  segment_tools_t segment_tools_obj;
  shift_segment_tools_t shift_segment_tools_obj;
  swap_segment_tools_t swap_segment_tools_obj;

  FUNC_LIB::function<vertex_vertex_matrix_type, nu> M;

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
MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::MC_walker(parameters_type& parameters_ref,
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
MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::~MC_walker() {
  if (concurrency.id() == 0 and thread_id == 0) {
    std::stringstream ss;
    ss << "\n\n\t\t walker died --> nb_successfull_updates/nb_updates : "
       << double(nb_successfull_updates) / double(nb_updates) << "\n\n";
    std::cout << ss.str();
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::initialize() {
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
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::test_interpolation() {
  std::cout << __FUNCTION__ << std::endl;

  SHOW::execute_on_bands(F_r_t);

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

        SHOW::plot_points(x, y);
      }
    }
  }

  throw std::logic_error(__FUNCTION__);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::do_sweep() {
  double factor = 1.;
  if (thermalized)
    factor = parameters.get_number_of_sweeps_per_measurement();

  int nr_of_segments = std::max(16, configuration.size());

  double ratio = double(nb_successfull_updates + 1) / double(nb_updates + 1);
  int factor2 = std::max(1, int(std::ceil(1.0 / ratio)));

  int nb_steps = nr_of_segments * factor * factor2;

  for (int l = 0; l < nb_steps; l++)
    do_step();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
int MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::get_random_interacting_flavor() {
  int spin = s::dmn_size() * rng();
  int int_band = parameters.get_interacting_bands().size() * rng();

  return parameters.get_interacting_bands()[int_band] + spin * b::dmn_size();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::do_step() {
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
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::do_insert_remove(int so_ind) {
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
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::insert_or_remove_full_line(int j) {
  nb_updates += 1;

  bool succes = full_line_tools_obj.insert_or_remove(j, mu(j));

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::insert_or_remove_anti_segment(int j) {
  nb_updates += 1;

  bool succes;
  if (rng() < 0.5)
    succes = anti_segment_tools_obj.insert_anti_segment(j, mu(j), sign, M, F_r_t);
  else
    succes = anti_segment_tools_obj.remove_anti_segment(j, mu(j), sign, M, F_r_t);

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::insert_or_remove_segment(int j) {
  nb_updates += 1;

  bool succes;
  if (rng() < 0.5)
    succes = segment_tools_obj.insert_segment(j, mu(j), sign, M, F_r_t);
  else
    succes = segment_tools_obj.remove_segment(j, mu(j), sign, M, F_r_t);

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::shift_segment(int j) {
  nb_updates += 1;

  bool succes;
  if (rng() < 0.5)
    succes = shift_segment_tools_obj.shift_segment_start_vertex(j, mu(j), sign, M, F_r_t);
  else
    succes = shift_segment_tools_obj.shift_segment_end_vertex(j, mu(j), sign, M, F_r_t);

  nb_successfull_updates += succes ? 1 : 0;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::swap_random_orbitals() {
  nb_updates += 1;

  int i = get_random_interacting_flavor();
  int j = get_random_interacting_flavor();

  bool succes = swap_segment_tools_obj.swap_orbitals(i, j, mu, sign, M, F_r_t);
  nb_successfull_updates += succes ? 1 : 0;
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_WALKER_H
