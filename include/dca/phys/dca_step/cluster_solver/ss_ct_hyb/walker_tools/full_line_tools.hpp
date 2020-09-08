// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// This class implements the insertion and removal of a full line. The insertion of a full line is
// only possible if no segement is present on that line.
//
// \f{eqnarray}{
//   W_{acc}(empty \rightarrow full) = min \left(1, \frac{W_{loc}(full)}{W_{loc}(empty)} \right)
// \f}

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_FULL_LINE_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_FULL_LINE_TOOLS_HPP

#include <cmath>
#include <cstdlib>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <typename hybridization_routines_type>
class full_line_tools {
public:
  typedef typename hybridization_routines_type::parameters_type parameters_type;
  typedef typename hybridization_routines_type::MOMS_type MOMS_type;
  typedef typename hybridization_routines_type::configuration_type configuration_type;
  typedef typename hybridization_routines_type::rng_type rng_type;

  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;

public:
  full_line_tools(hybridization_routines_type& hybridization_routines_ref);

  bool insert_or_remove(int j, double mu);

private:
  bool insert_full_line(int j, double mu);
  bool remove_full_line(int j, double mu);

  double get_other_length_u(int j);

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
full_line_tools<hybridization_routines_type>::full_line_tools(
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
bool full_line_tools<hybridization_routines_type>::insert_or_remove(int j, double mu) {
  if (configuration.get_vertices(j).size())
    return false;

  if (configuration.get_full_line(j) == true)
    return remove_full_line(j, mu);
  else
    return insert_full_line(j, mu);
}

template <typename hybridization_routines_type>
bool full_line_tools<hybridization_routines_type>::insert_full_line(int j, double mu) {
  //     cout << __FUNCTION__ << endl;

  if (configuration.get_full_line(j) == true)
    return false;

  double otherlength_u = get_other_length_u(j);

  if (std::log(rng()) < beta * mu - otherlength_u) {
    configuration.get_full_line(j) = true;
    return true;
  }
  return false;
}

template <typename hybridization_routines_type>
bool full_line_tools<hybridization_routines_type>::remove_full_line(int j, double mu) {
  //     cout << __FUNCTION__ << endl;

  if (configuration.get_full_line(j) == false)
    return false;

  double otherlength_u = get_other_length_u(j);

  if (std::log(rng()) < -beta * mu + otherlength_u) {
    configuration.get_full_line(j) = false;
    return true;
  }
  return false;
}

template <typename hybridization_routines_type>
double full_line_tools<hybridization_routines_type>::get_other_length_u(int j) {
  double otherlength_u = 0;

  for (int i = 0; i < spin_orbitals; i++) {
    if (i == j)
      continue;

    double other_length = 0;

    if (configuration.get_full_line(i) == true)
      other_length = beta;
    else {
      orbital_configuration_type& vertices = configuration.get_vertices(i);

      for (size_t l = 0; l < vertices.size(); l++)
        other_length += (vertices[l].t_end() - vertices[l].t_start() > 0
                             ? vertices[l].t_end() - vertices[l].t_start()
                             : vertices[l].t_end() - vertices[l].t_start() + beta);
    }

    otherlength_u += other_length * MOMS.H_interactions(i, j, 0);
  }

  return otherlength_u;
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_FULL_LINE_TOOLS_HPP
