// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// This class implements the insertion and removal of antisegments \f$(c-c^{\dagger} pair)\f$.
//
// \f{eqnarray}{
//   W_{acc}(c_k \rightarrow c_{k+1}) = min \left(1, \frac{\beta l_{max}}{k+1} \frac{det F(k+1)}{det
//   F(k)}\frac{W_{loc}(k+1)}{W_{loc}(k)} \right)
// \f}

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_ANTI_SEGMENT_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_ANTI_SEGMENT_TOOLS_HPP

#include <cmath>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/structures/hybridization_vertex.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <typename hybridization_routines_type>
class anti_segment_tools {
public:
  typedef typename hybridization_routines_type::parameters_type parameters_type;
  typedef typename hybridization_routines_type::MOMS_type MOMS_type;
  typedef typename hybridization_routines_type::configuration_type configuration_type;
  typedef typename hybridization_routines_type::rng_type rng_type;

  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;

public:
  anti_segment_tools(hybridization_routines_type& hybridization_routines_ref);

  template <typename function_type_1, typename function_type_2>
  bool insert_anti_segment(int this_flavor, double mu, double& sign, function_type_1& M,
                           function_type_2& F);

  template <typename function_type_1, typename function_type_2>
  bool remove_anti_segment(int this_flavor, double mu, double& sign, function_type_1& M,
                           function_type_2& F);

private:
  double get_other_length_u(int j, Hybridization_vertex& vertex);

  /*******************
   ***  INSERTION  ***
   *******************/

  template <typename function_type_1, typename function_type_2>
  bool insert_anti_segment_into_full_line(int this_flavor, double mu, double& sign,
                                          function_type_1& M, function_type_2& F);

  template <typename function_type_1, typename function_type_2>
  bool insert_anti_segment_into_segmented_line(int this_flavor, double mu, double& sign,
                                               function_type_1& M, function_type_2& F);

  /*******************
   ***  REMOVAL    ***
   *******************/

  template <typename function_type_1, typename function_type_2>
  bool remove_anti_segment_for_one_segment(int this_flavor, double mu, double& sign,
                                           function_type_1& M, function_type_2& F);

  template <typename function_type_1, typename function_type_2>
  bool remove_anti_segment_for_multiple_segments(int this_flavor, double mu, double& sign,
                                                 function_type_1& M, function_type_2& F);

private:
  hybridization_routines_type& hybridization_routines;

  const parameters_type& parameters;
  MOMS_type& MOMS;

  configuration_type& configuration;
  rng_type& rng;

  int FLAVORS;
  double BETA;
};

template <typename hybridization_routines_type>
anti_segment_tools<hybridization_routines_type>::anti_segment_tools(
    hybridization_routines_type& hybridization_routines_ref)
    :

      hybridization_routines(hybridization_routines_ref),

      parameters(hybridization_routines.get_parameters()),
      MOMS(hybridization_routines.get_MOMS()),

      configuration(hybridization_routines.get_configuration()),
      rng(hybridization_routines.get_rng()) {
  FLAVORS = b::dmn_size() * s::dmn_size();
  BETA = parameters.get_beta();
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool anti_segment_tools<hybridization_routines_type>::insert_anti_segment(int this_flavor,
                                                                          double mu, double& sign,
                                                                          function_type_1& M,
                                                                          function_type_2& F) {
  if (configuration.get_full_line(this_flavor) == true) {
    return insert_anti_segment_into_full_line(this_flavor, mu, sign, M, F);
  }
  else {
    if (configuration.get_vertices(this_flavor).size() > 0) {
      return insert_anti_segment_into_segmented_line(this_flavor, mu, sign, M, F);
    }
    else {
      return false;
    }
  }
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool anti_segment_tools<hybridization_routines_type>::insert_anti_segment_into_full_line(
    int this_flavor, double mu, double& sign, function_type_1& M, function_type_2& F) {
  double t = BETA * rng();
  double max_length = BETA;

  double length = hybridization_routines.compute_length(rng(), max_length, 0);

  double t_end = (t + length < BETA ? t + length : t + length - BETA);

  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  Hybridization_vertex segment_insert(t_end, t);
  Hybridization_vertex segment_remove(t, t_end);

  double log_prob, overlap, det_rat, det_rat_sign;
  std::vector<double> R(vertices.size()), Q(vertices.size());

  double otherlength_u = get_other_length_u(this_flavor, segment_remove);
  det_rat = hybridization_routines.det_rat_up(this_flavor, segment_insert, M(this_flavor), vertices,
                                              F, R, Q, det_rat_sign, overlap);

  log_prob = std::log(BETA * max_length * det_rat) - length * mu + otherlength_u;

  if (std::log(rng()) < log_prob) {
    hybridization_routines.compute_M_up(0, 0, M(this_flavor), vertices, F, R, Q, det_rat * overlap);
    sign *= det_rat_sign;
    vertices.push_back(segment_insert);
    configuration.get_full_line(this_flavor) = false;

    // return exp(-length*mu+otherlength_u);
    return true;
  }

  // return -1;
  return false;
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool anti_segment_tools<hybridization_routines_type>::insert_anti_segment_into_segmented_line(
    int this_flavor, double mu, double& /*sign*/, function_type_1& M, function_type_2& F) {
  double t = BETA * rng();

  double t_up;    // distance to next segment up (t_start)
  double t_down;  // distance to next segment down (t_end)

  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  typename orbital_configuration_type::iterator s_up;    // iterator of the segment up
  typename orbital_configuration_type::iterator s_down;  // iterator of the segment down

  hybridization_routines.compute_intervals(t, BETA, t_up, t_down, vertices, s_up, s_down);

  if (t_down <
      0) {  // t does lie on a segment -> it's possible to insert an anti-segment starting from t
    // double length = BETA*rng();
    double length = hybridization_routines.compute_length(rng(), -t_down, 0);
    if (length < -t_down) {
      Hybridization_vertex segment_shrink(s_down->t_start(), t);

      double t_start = t + length;
      if (t_start > BETA)
        t_start -= BETA;

      Hybridization_vertex segment_insert(t_start, s_down->t_end());
      Hybridization_vertex anti_segment(t, t_start);
      Hybridization_vertex remove_segment(t_start, t);

      double otherlength_u = get_other_length_u(this_flavor, anti_segment);

      double log_prob, overlap, det_rat, det_rat_sign;
      std::vector<double> Q(vertices.size()), R(vertices.size());

      det_rat = hybridization_routines.det_rat_up(this_flavor, remove_segment, M(this_flavor),
                                                  vertices, F, R, Q, det_rat_sign, overlap);

      log_prob =
          std::log(BETA * (-t_down) / (vertices.size() + 1) * det_rat) - length * mu + otherlength_u;
      // log_prob = std::log(BETA*(BETA)/(vertices.size()+1)*det_rat)-length*mu+otherlength_u;
      if (std::log(rng()) < log_prob) {
        int s, r;  // s is the segment which is shifted, r the segment which is inserted
        s = 0;
        for (typename orbital_configuration_type::iterator it = vertices.begin(); it != s_down; it++)
          s++;
        if (anti_segment.t_end() > segment_shrink.t_start())
          r = s + 1;
        else {
          r = 0;
          s++;
        }

        hybridization_routines.compute_M_up(r, s, M(this_flavor), vertices, F, R, Q,
                                            det_rat * overlap);

        s_down->set_t_end(t);
        if (segment_insert.t_start() > vertices.begin()->t_start()) {
          s_down++;
          vertices.insert(s_down, segment_insert);
        }
        else {
          vertices.insert(vertices.begin(), segment_insert);
        }

        // return exp(-length*mu+otherlength_u);
        return true;
      }
    }
  }

  // return -1;
  return false;
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool anti_segment_tools<hybridization_routines_type>::remove_anti_segment(int this_flavor,
                                                                          double mu, double& sign,
                                                                          function_type_1& M,
                                                                          function_type_2& F) {
  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  switch (vertices.size()) {
    case 0:
      return false;
      break;

    case 1:
      return remove_anti_segment_for_one_segment(this_flavor, mu, sign, M, F);
      break;

    default:
      return remove_anti_segment_for_multiple_segments(this_flavor, mu, sign, M, F);
  };
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool anti_segment_tools<hybridization_routines_type>::remove_anti_segment_for_one_segment(
    int this_flavor, double mu, double& /*sign*/, function_type_1& M, function_type_2& /*F*/) {
  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);
  assert(vertices.size() == 1);
  typename orbital_configuration_type::iterator s_down = vertices.begin();

  double det_rat = std::fabs(M(this_flavor)(0, 0));

  double length = s_down->t_start() - s_down->t_end();
  if (length < 0)
    length += BETA;

  Hybridization_vertex anti_segment(s_down->t_end(), s_down->t_start());

  double otherlength_u = get_other_length_u(this_flavor, anti_segment);

  double log_prob = std::log((BETA * BETA / det_rat)) - length * mu + otherlength_u;

  if (std::log(rng()) < -log_prob) {
    configuration.get_full_line(this_flavor) = true;
    vertices.erase(s_down);
    hybridization_routines.compute_M_down(0, 0, M(this_flavor));

    // return exp(length*mu-otherlength_u);
    return true;
  }

  // return -1;
  return false;
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool anti_segment_tools<hybridization_routines_type>::remove_anti_segment_for_multiple_segments(
    int this_flavor, double mu, double& /*sign*/, function_type_1& M, function_type_2& /*F*/) {
  typename orbital_configuration_type::iterator s_up;    // iterator of the segment up
  typename orbital_configuration_type::iterator s_down;  // iterator of the segment down

  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  assert(vertices.size() > 1);

  int r = rng() * vertices.size();
  s_up = vertices.begin();
  for (int i = 0; i < r; i++)
    s_up++;

  int s = r - 1;
  if (s < 0) {
    s = vertices.size() - 1;
    s_down = vertices.end();
    s_down--;
  }
  else {
    s_down = s_up;
    s_down--;
  }

  double length = s_up->t_start() - s_down->t_end();
  if (length < 0)
    length += BETA;

  double t_total = s_up->t_end() - s_down->t_end();
  if (t_total < 0)
    t_total += BETA;

  Hybridization_vertex anti_segment(s_down->t_end(), s_up->t_start());

  double otherlength_u = get_other_length_u(this_flavor, anti_segment);

  double log_prob, det_rat, det_rat_sign;

  det_rat = hybridization_routines.det_rat_down(r, s, M(this_flavor), vertices, det_rat_sign);

  log_prob = std::log(BETA * t_total / vertices.size() / det_rat) - length * mu + otherlength_u;
  // log_prob = std::log(BETA*BETA/vertices.size()/det_rat)-length*mu+otherlength_u;

  if (std::log(rng()) < -log_prob) {
    hybridization_routines.compute_M_down(r, s, M(this_flavor));

    double t_end = s_up->t_end();
    vertices.erase(s_up);

    if (r > 0) {
      s_up = vertices.begin();
      for (int k = 0; k < s; k++)
        s_up++;
    }
    else {
      s = vertices.size() - 1;
      s_up = vertices.end();
      s_up--;
    }

    s_up->set_t_end(t_end);

    // return exp(length*mu-otherlength_u);
    return true;
  }

  // return -1;
  return false;
}

template <typename hybridization_routines_type>
double anti_segment_tools<hybridization_routines_type>::get_other_length_u(
    int j, Hybridization_vertex& vertex) {
  double otherlength_u = 0;

  for (int i = 0; i < FLAVORS; i++) {
    if (i == j)
      continue;

    double other_length = hybridization_routines.compute_overlap(
        vertex, configuration.get_vertices(i), configuration.get_full_line(i), BETA);
    otherlength_u += other_length * MOMS.H_interactions(i, j, 0);
  }

  return otherlength_u;
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_ANTI_SEGMENT_TOOLS_HPP
