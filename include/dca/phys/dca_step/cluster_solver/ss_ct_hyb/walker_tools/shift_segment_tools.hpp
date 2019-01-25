// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
// Modified by Andrei Plamada (plamada@itp.phys.ethz.ch).
//
// This class implements the shifting of a \f$c^{\dagger}\f$.
//
// \f{eqnarray}{
//   W_{acc}(c_k \rightarrow c'_k) = min \left(1, \frac{det F'}{det F}\frac{W'_{loc}}{W_{loc}}
//   \right)
// \f}

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SHIFT_SEGMENT_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SHIFT_SEGMENT_TOOLS_HPP

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
class shift_segment_tools {
public:
  typedef typename hybridization_routines_type::parameters_type parameters_type;
  typedef typename hybridization_routines_type::MOMS_type MOMS_type;
  typedef typename hybridization_routines_type::configuration_type configuration_type;
  typedef typename hybridization_routines_type::rng_type rng_type;

  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;

public:
  shift_segment_tools(hybridization_routines_type& hybridization_routines_ref);

  template <typename function_type_1, typename function_type_2>
  bool shift_segment_end_vertex(int this_flavor, double mu, double& sign, function_type_1& M,
                                function_type_2& F);

  template <typename function_type_1, typename function_type_2>
  bool shift_segment_start_vertex(int this_flavor, double mu, double& sign, function_type_1& M,
                                  function_type_2& F);

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
shift_segment_tools<hybridization_routines_type>::shift_segment_tools(
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
bool shift_segment_tools<hybridization_routines_type>::shift_segment_end_vertex(
    int this_flavor, double mu, double& sign, function_type_1& M, function_type_2& F) {
  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  int size = vertices.size();

  if (size < 1)
    return -1;

  int n = size * rng();

  typename orbital_configuration_type::iterator s, s_up;
  s = vertices.begin() + n;
  s_up = vertices.begin() + (n + 1) % size;

  double interval = s_up->t_start() - s->t_start();
  if (interval <= 0)
    interval += BETA;

  double length = hybridization_routines.compute_length(rng(), interval, 0);
  double length_old = s->t_end() - s->t_start();
  if (length_old < 0)
    length_old += BETA;

  double new_t_end = s->t_start() + length;
  if (new_t_end > BETA)
    new_t_end -= BETA;

  Hybridization_vertex segment_insert(s->t_start(), new_t_end);
  Hybridization_vertex segment_remove(s->t_start(), s->t_end());

  double otherlength_u = 0;
  for (int i = 0; i < FLAVORS; i++) {
    if (i == this_flavor)
      continue;

    double other_length =
        hybridization_routines.compute_overlap(segment_insert, configuration.get_vertices(i),
                                               configuration.get_full_line(i), BETA) -
        hybridization_routines.compute_overlap(segment_remove, configuration.get_vertices(i),
                                               configuration.get_full_line(i), BETA);

    otherlength_u += other_length * MOMS.H_interactions(i, this_flavor, 0);
  }

  double det_rat, det_rat_sign, overlap;
  std::vector<double> R(vertices.size()), Q(vertices.size());

  det_rat = hybridization_routines.det_rat_shift_end(this_flavor, segment_insert, n, M(this_flavor),
                                                     vertices, F, R, Q, det_rat_sign, overlap);

  if (std::log(rng()) < std::log(det_rat) + (length - length_old) * mu - otherlength_u) {
    if (det_rat_sign < 0) {
      M(this_flavor).print();
    }
    hybridization_routines.compute_M_shift_end(n, M(this_flavor), R, Q, det_rat * overlap);
    sign *= det_rat_sign;
    s->set_t_end(new_t_end);

    // return exp((length-length_old)*mu-otherlength_u);
    return true;
  }

  // return -1;
  return false;
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
bool shift_segment_tools<hybridization_routines_type>::shift_segment_start_vertex(
    int this_flavor, double mu, double& sign, function_type_1& M, function_type_2& F) {
  orbital_configuration_type& vertices = configuration.get_vertices(this_flavor);

  int size = vertices.size();

  if (size < 1)
    return -1;

  int n = size * rng();

  typename orbital_configuration_type::iterator s, s_up;
  s = vertices.begin() + n;
  s_up = vertices.begin() + (n + 1) % size;

  n = (n + 1) % size;

  double interval = s_up->t_end() - s->t_end();
  if (interval <= 0)
    interval += BETA;

  double length = hybridization_routines.compute_length(rng(), interval, 0);
  double length_old = s_up->t_end() - s_up->t_start();
  if (length_old < 0)
    length_old += BETA;

  double new_t_start = s_up->t_end() - length;
  if (new_t_start < 0)
    new_t_start += BETA;

  Hybridization_vertex segment_insert(new_t_start, s_up->t_end());
  Hybridization_vertex segment_remove(s_up->t_start(), s_up->t_end());

  double otherlength_u = 0;
  for (int i = 0; i < FLAVORS; i++) {
    if (i == this_flavor)
      continue;

    double other_length =
        hybridization_routines.compute_overlap(segment_insert, configuration.get_vertices(i),
                                               configuration.get_full_line(i), BETA) -
        hybridization_routines.compute_overlap(segment_remove, configuration.get_vertices(i),
                                               configuration.get_full_line(i), BETA);

    otherlength_u += other_length * MOMS.H_interactions(i, this_flavor, 0);
  }

  double det_rat, det_rat_sign, overlap;
  std::vector<double> R(vertices.size()), Q(vertices.size());

  det_rat = hybridization_routines.det_rat_shift_start(
      this_flavor, segment_insert, n, M(this_flavor), vertices, F, R, Q, det_rat_sign, overlap);

  if (std::log(rng()) < std::log(det_rat) + (length - length_old) * mu - otherlength_u) {
    if (det_rat_sign < 0) {
      M(this_flavor).print();
    }
    hybridization_routines.compute_M_shift_start(n, M(this_flavor), R, Q, det_rat * overlap);
    sign *= det_rat_sign;
    s_up->set_t_start(new_t_start);

    if (n == 0 and
        (segment_insert.t_end() - segment_insert.t_start()) *
                (segment_remove.t_end() - segment_remove.t_start()) <
            0) {
      hybridization_routines.cycle_column_backward(M(this_flavor));
      hybridization_routines.cycle_row_backward(M(this_flavor));

      Hybridization_vertex aux;
      aux = vertices[0];

      for (int i = 0; i < size - 1; i++)
        vertices[i] = vertices[i + 1];

      vertices[size - 1] = aux;
    }
    if (n == size - 1 and
        (segment_insert.t_end() - segment_insert.t_start()) *
                (segment_remove.t_end() - segment_remove.t_start()) <
            0) {
      hybridization_routines.cycle_column_forward(M(this_flavor));
      hybridization_routines.cycle_row_forward(M(this_flavor));

      Hybridization_vertex aux;
      aux = vertices[size - 1];

      for (int i = size - 1; i > 0; i--)
        vertices[i] = vertices[i - 1];

      vertices[0] = aux;
    }

    // return exp((length-length_old)*mu-otherlength_u);
    return true;
  }

  // return -1;
  return false;
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SHIFT_SEGMENT_TOOLS_HPP
