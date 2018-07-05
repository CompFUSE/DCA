// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides cluster specifications.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SPECIFICATIONS_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SPECIFICATIONS_HPP

#include "dca/math/function_transform/domain_specifications.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename scalar_type, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
struct cluster_specifications {};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, CLUSTER, MOMENTUM_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::DISCRETE,
      math::transform::KRONECKER_DELTA, math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, CLUSTER, REAL_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::EXPANSION, math::transform::HARMONICS,
      math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, LATTICE_SP, MOMENTUM_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::DISCRETE,
      math::transform::KRONECKER_DELTA, math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, LATTICE_SP, REAL_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::EXPANSION, math::transform::HARMONICS,
      math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, LATTICE_TP, MOMENTUM_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::DISCRETE,
      math::transform::KRONECKER_DELTA, math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, LATTICE_TP, REAL_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::EXPANSION, math::transform::HARMONICS,
      math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, VASP_LATTICE, MOMENTUM_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::DISCRETE,
      math::transform::KRONECKER_DELTA, math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, VASP_LATTICE, REAL_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::EXPANSION, math::transform::HARMONICS,
      math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, TMP_CLUSTER, MOMENTUM_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::DISCRETE,
      math::transform::KRONECKER_DELTA, math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

template <typename scalar_type, CLUSTER_SHAPE S>
struct cluster_specifications<scalar_type, TMP_CLUSTER, REAL_SPACE, S> {
  typedef math::transform::domain_specifications<
      scalar_type, std::vector<scalar_type>, math::transform::EXPANSION, math::transform::HARMONICS,
      math::transform::PERIODIC, math::transform::EQUIDISTANT>
      dmn_specifications_type;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SPECIFICATIONS_HPP
