// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
// Modified by Andrei Plamada (plamada@itp.phys.ethz.ch).
//
// This class implements the swap of segemnts between two hybridization lines.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SWAP_SEGMENT_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SWAP_SEGMENT_TOOLS_HPP

#include <cmath>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/structures/hybridization_vertex.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <typename hybridization_routines_type>
class swap_segment_tools {
public:
  typedef dca::linalg::Matrix<double, dca::linalg::CPU> vertex_vertex_matrix_type;

  typedef typename hybridization_routines_type::parameters_type parameters_type;
  typedef typename hybridization_routines_type::MOMS_type MOMS_type;
  typedef typename hybridization_routines_type::configuration_type configuration_type;
  typedef typename hybridization_routines_type::rng_type rng_type;

  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

public:
  swap_segment_tools(hybridization_routines_type& hybridization_routines_ref);

  template <typename function_type_0, typename function_type_1, typename function_type_2>
  bool swap_orbitals(int i, int j, function_type_0& mu, double& sign, function_type_1& M,
                     function_type_2& F);

  template <typename function_type_1, typename function_type_2>
  double construct_inverse(function_type_1& M, double beta, function_type_2& F, int flavor_1,
                           int flavor_2);

  template <typename matrix_type_1, typename function_type_2>
  void construct_matrix(matrix_type_1& M, double beta, function_type_2& F, int flavor_1,
                        int flavor_2);

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
swap_segment_tools<hybridization_routines_type>::swap_segment_tools(
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
template <typename function_type_0, typename function_type_1, typename function_type_2>
bool swap_segment_tools<hybridization_routines_type>::swap_orbitals(
    int i, int j, function_type_0& mu, double& /*sign*/, function_type_1& M, function_type_2& F) {
  orbital_configuration_type& vertices = configuration.get_vertices(i);
  orbital_configuration_type& swap_vertices = configuration.get_vertices(j);

  bool nonsymmetric_orbitals = false;
  double det_rat;

  if (nonsymmetric_orbitals) {
    vertex_vertex_matrix_type M_new_this, M_new_swap;

    double det_old_this = construct_inverse(M_new_this, beta, F, i, i);  // before swap
    double det_old_swap = construct_inverse(M_new_swap, beta, F, j, j);
    double det_new_this = construct_inverse(M_new_this, beta, F, i, j);  // after swap
    double det_new_swap = construct_inverse(M_new_swap, beta, F, j, i);
    det_rat = (det_new_this / det_old_this) * (det_new_swap / det_old_swap);
  }
  else {
    det_rat = 1.0;
  }
  Hybridization_vertex full_segment(0, beta);

  double length_this = hybridization_routines.compute_overlap(full_segment, vertices,
                                                              configuration.get_full_line(i), beta);
  double length_swap = hybridization_routines.compute_overlap(full_segment, swap_vertices,
                                                              configuration.get_full_line(j), beta);
  double e_site =
      length_this * (mu(i) - mu(j)) + length_swap * (mu(j) - mu(i));  // site energy contribution

  double overlap_u = 0;
  for (int k = 0; k < spin_orbitals; k++) {
    if (k != i && k != j) {
      double overlap = 0;
      if (configuration.get_full_line(i)) {
        overlap += hybridization_routines.compute_overlap(
            full_segment, configuration.get_vertices(k), configuration.get_full_line(k), beta);
      }
      else
        for (typename orbital_configuration_type::iterator it = vertices.begin();
             it != vertices.end(); it++) {
          overlap += hybridization_routines.compute_overlap(*it, configuration.get_vertices(k),
                                                            configuration.get_full_line(k), beta);
        }
      overlap_u += (MOMS.H_interactions(j, k, 0) - MOMS.H_interactions(i, k, 0)) *
                   overlap;  //(u(j,k)-u(i,k))*overlap;
      overlap = 0;
      if (configuration.get_full_line(j)) {
        overlap += hybridization_routines.compute_overlap(
            full_segment, configuration.get_vertices(k), configuration.get_full_line(k), beta);
      }
      else
        for (typename orbital_configuration_type::iterator it = swap_vertices.begin();
             it != swap_vertices.end(); it++) {
          overlap += hybridization_routines.compute_overlap(*it, configuration.get_vertices(k),
                                                            configuration.get_full_line(k), beta);
        }
      overlap_u += (MOMS.H_interactions(i, k, 0) - MOMS.H_interactions(j, k, 0)) *
                   overlap;  //(u(i,k)-u(j,k))*overlap;
    }
  }

  double log_prob = std::log(std::fabs(det_rat)) + (-overlap_u - e_site);
  if (std::log(rng()) < log_prob) {
    if (nonsymmetric_orbitals) {
      construct_inverse(M(i), beta, F, i, j);
      construct_inverse(M(j), beta, F, j, i);
    }
    else {
      vertex_vertex_matrix_type M_aux;
      M_aux = M(i);
      M(i).copy_from(M(j));
      M(j).copy_from(M_aux);
    }

    swap(vertices, swap_vertices);
    int dummy1 = configuration.get_full_line(i);
    configuration.get_full_line(i) = configuration.get_full_line(j);
    configuration.get_full_line(j) = dummy1;

    // return exp(-overlap_u-e_site);
    return true;
  }
  return false;
  //    return -1;
}

template <typename hybridization_routines_type>
template <typename function_type_1, typename function_type_2>
double swap_segment_tools<hybridization_routines_type>::construct_inverse(
    function_type_1& M, double beta, function_type_2& F, int flavor_1, int flavor_2) {
  double det = 1;
  orbital_configuration_type segments = configuration.get_vertices(flavor_2);

  if (segments.size() > 0) {
    construct_matrix(M, beta, F, flavor_1, flavor_2);

    dca::linalg::matrixop::inverse(M);

    for (int i = 0; i < M.size().first; i++)
      det *= M(i, i);
  }
  else
    M.resize(0);

  return det;
}

template <typename hybridization_routines_type>
template <typename matrix_type_1, typename function_type_2>
void swap_segment_tools<hybridization_routines_type>::construct_matrix(matrix_type_1& M, double beta,
                                                                       function_type_2& F,
                                                                       int flavor_1, int flavor_2) {
  orbital_configuration_type segments = configuration.get_vertices(flavor_2);

  int N = segments.size();

  matrix_type_1 M_new(N, N);

  int row = -1;
  int column = -1;

  int coor[2];
  nu nu_obj;

  nu_obj.linind_2_subind(flavor_1, coor);

  for (typename orbital_configuration_type::iterator it1 = segments.begin(); it1 != segments.end();
       it1++) {
    row++;
    for (typename orbital_configuration_type::iterator it2 = segments.begin();
         it2 != segments.end(); it2++) {
      column++;

      double argument = it1->t_end() - it2->t_start();
      double sign = 1;
      if (argument < 0) {
        argument += beta;
        sign = -1;
      }
      M_new(row, column) = hybridization_routines.interpolate_F(coor, argument, F) * sign;
    }
    column = -1;
  }
  M.copy_from(M_new);
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SWAP_SEGMENT_TOOLS_HPP
