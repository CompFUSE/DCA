// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// This class implements the insertion and removal of segments \f$(c^{\dagger}-c\ pair)\f$.
//
// \f{eqnarray}{
//   W_{acc}(c_k \rightarrow c_{k+1}) = min \left(1, \frac{\beta l_{max}}{k+1} \frac{det F(k+1)}{det
//   F(k)}\frac{W_{loc}(k+1)}{W_{loc}(k)} \right)
// \f}

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SEGMENT_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SEGMENT_TOOLS_HPP

#include <cassert>
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
class segment_tools {
public:
  typedef typename hybridization_routines_type::parameters_type parameters_type;
  typedef typename hybridization_routines_type::MOMS_type MOMS_type;
  typedef typename hybridization_routines_type::configuration_type configuration_type;
  typedef typename hybridization_routines_type::rng_type rng_type;

  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;

public:
  segment_tools(hybridization_routines_type& hybridization_routines_ref);

  template <typename function_type_1, typename function_type_2>
  bool insert_segment(int this_flavor, double mu, double& sign, function_type_1& M,
                      function_type_2& F);

  template <typename function_type_1, typename function_type_2>
  bool remove_segment(int this_flavor, double mu, double& sign, function_type_1& M,
                      function_type_2& F);

private:
  hybridization_routines_type& hybridization_routines;

  const parameters_type& parameters;
  MOMS_type& MOMS;

  configuration_type& configuration;
  rng_type& rng;

  int spin_orbitals;
  double beta;
};

template <typename hybridization_routines_type>
segment_tools<hybridization_routines_type>::segment_tools(
    hybridization_routines_type& hybridization_routines_ref)
    :

      hybridization_routines(hybridization_routines_ref),

      parameters(hybridization_routines.get_parameters()),
      MOMS(hybridization_routines.get_MOMS()),

      configuration(hybridization_routines.get_configuration()),
      rng(hybridization_routines.get_rng()) {
  spin_orbitals = b::dmn_size() * s::dmn_size();
  beta = parameters.get_beta();
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool segment_tools<hybridization_routines_type>::insert_segment(int this_flavor, double mu,
                                                                double& sign, function_type_1& M,
                                                                function_type_2& F) {
  //     cout << __FUNCTION__ << endl;
  assert(this_flavor >= 0 && this_flavor < spin_orbitals);

  if (configuration.get_full_line(this_flavor) == true)
    return -1;

  double t = beta * rng();
  // cout << "t : " << t << endl;

  double t_up;    // distance to next segment up
  double t_down;  // distance to next segment down

  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  typename orbital_configuration_type::iterator s_up;    // iterator of the segment up
  typename orbital_configuration_type::iterator s_down;  // iterator of the segment down

  hybridization_routines.compute_intervals(t, beta, t_up, t_down, vertices, s_up, s_down);

  // double trace = -1;
  bool succes = false;

  if (t_down >
      0) {  // t does not lie on a segment -> it's possible to insert a new one starting from t
    // double length = beta*rng();
    double length = hybridization_routines.compute_length(rng(), t_up, 0);
    // cout << "length : " << length << endl;
    if (length < t_up) {
      Hybridization_vertex segment_insert;
      segment_insert.set_t_start(t);
      double t_final = t + length;
      if (t_final > beta)
        segment_insert.set_t_end(t_final - beta);
      else
        segment_insert.set_t_end(t_final);

      double otherlength_u = 0;

      for (int i = 0; i < spin_orbitals; i++) {
        if (i == this_flavor)
          continue;
        double other_length = hybridization_routines.compute_overlap(
            segment_insert, configuration.get_vertices(i), configuration.get_full_line(i),
            beta);  // compute_overlap(segment_insert, other_segments[i], other_full_line[i], beta);
        otherlength_u +=
            other_length * MOMS.H_interactions(i, this_flavor, 0);  // u(i, this_flavor);
      }
      // cout << "otherlength_u : " << otherlength_u << endl;

      double log_prob, overlap, det_rat, det_rat_sign;
      std::vector<double> R(vertices.size()), Q(vertices.size());

      det_rat = hybridization_routines.det_rat_up(this_flavor, segment_insert, M(this_flavor),
                                                  vertices, F, R, Q, det_rat_sign, overlap);
      // cout << "det_rat : " << det_rat << endl;
      // log_prob = std::log(beta*beta/(vertices.size()+1)*det_rat)+mu*length-otherlength_u;
      log_prob =
          std::log(beta * t_up / (vertices.size() + 1) * det_rat) + mu * length - otherlength_u;
      // cout << "log_prob : " << log_prob << endl;
      if (std::log(rng()) < log_prob) {
        int position = 0;
        for (typename orbital_configuration_type::iterator it = vertices.begin(); it != s_up; it++)
          position++;

        hybridization_routines.compute_M_up(position, position, M(this_flavor), vertices, F, R, Q,
                                            det_rat * overlap);
        sign *= det_rat_sign;
        vertices.insert(s_up, segment_insert);

        // trace = std::exp(mu*length-otherlength_u);
        succes = true;
      }
    }
  }

  // return trace;
  return succes;
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool segment_tools<hybridization_routines_type>::remove_segment(int this_flavor, double mu,
                                                                double& sign, function_type_1& M,
                                                                function_type_2& /*F*/) {
  assert(this_flavor >= 0 && this_flavor < spin_orbitals);

  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  if (vertices.size() == 0)
    return -1;

  typename orbital_configuration_type::iterator s_up;    // iterator of the segment up
  typename orbital_configuration_type::iterator s_down;  // iterator of the segment down

  int position = rng() * vertices.size();

  s_down = vertices.begin();
  for (int i = 0; i < position; i++)
    s_down++;
  s_up = s_down;
  s_up++;
  if (s_up == vertices.end())
    s_up = vertices.begin();

  double length = s_down->t_end() - s_down->t_start();
  if (length < 0)
    length += beta;

  double t_total = s_up->t_start() - s_down->t_start();
  if (t_total <= 0)
    t_total += beta;

  Hybridization_vertex segment_remove = *s_down;

  double otherlength_u = 0;
  for (int i = 0; i < spin_orbitals; i++) {
    if (i == this_flavor)
      continue;
    double other_length = hybridization_routines.compute_overlap(
        segment_remove, configuration.get_vertices(i), configuration.get_full_line(i),
        beta);  // compute_overlap(segment_remove, other_segments[i], other_full_line[i], beta);
    otherlength_u += other_length * MOMS.H_interactions(i, this_flavor, 0);
  }

  double log_prob, det_rat, det_rat_sign;

  det_rat = hybridization_routines.det_rat_down(position, position, M(this_flavor), vertices,
                                                det_rat_sign);

  // log_prob = std::log(beta*beta/vertices.size()/det_rat)+length*mu-otherlength_u;
  log_prob = std::log(beta * t_total / vertices.size() / det_rat) + length * mu - otherlength_u;

  if (std::log(rng()) < -log_prob) {
    hybridization_routines.compute_M_down(position, position, M(this_flavor));
    sign *= det_rat_sign;
    vertices.erase(s_down);
    // return std::exp(-mu*length+otherlength_u);

    return true;
  }

  // return -1;
  return false;
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SEGMENT_TOOLS_HPP
