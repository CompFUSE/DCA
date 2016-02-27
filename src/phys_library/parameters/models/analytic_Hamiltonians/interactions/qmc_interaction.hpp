// Copyright 2016 ETH Zurich.
//
// Author: Urs Haehner (haehneru@itp.phys.ethz.ch), ETH Zurich
//
// First draft.

#ifndef QMC_INTERACTION_HPP
#define QMC_INTERACTION_HPP

#include <array>
#include <vector>
#include "ctaux_vertex_pair.h"

template <typename parameters_type, typename H_interaction_type>
class QMC_interaction {
public:
  using rng_type = typename parameters_type::rng_type;
  using e_band_dmn = typename parameters_type::b;
  using e_spin_dmn = typename parameters_type::s;
  using r_dmn = typename parameters_type::r_DCA;
  using nu = typename parameters_type::nu;

  QMC_interaction(parameters_type& parameters_ref, rng_type& rng_ref, H_interaction_type& H_interaction);
  
  void set_vertex();
  
private:
  void initialize_correlated_orbitals();
  
  parameters_type& parameters;
  rng_type& rng;

  std::vector<int> correlated_orbitals;
};

#endif  // QMC_INTERACTION_HPP
