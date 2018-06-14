// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides cluster operations.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_OPERATIONS_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_OPERATIONS_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class cluster_operations {
public:
  template <typename scalar_type>
  static int index(const std::vector<scalar_type>& element,
                   const std::vector<std::vector<scalar_type>>& elements, const CLUSTER_SHAPE shape);

  template <typename scalar_type>
  static int origin_index(const std::vector<std::vector<scalar_type>>& elements,
                          const CLUSTER_SHAPE shape);

  template <typename scalar_type>
  static std::vector<scalar_type> translate_inside_cluster(
      const std::vector<scalar_type>& r, const std::vector<std::vector<scalar_type>>& basis);

  template <typename scalar_type>
  static bool test_translate_inside_cluster(const std::vector<std::vector<scalar_type>>& elements,
                                            const std::vector<std::vector<scalar_type>>& basis);

  template <typename scalar_type>
  static scalar_type minimal_distance(std::vector<scalar_type> vec_0, std::vector<scalar_type> vec_1,
                                      const std::vector<std::vector<scalar_type>>& basis);
  template <typename scalar_type>
  static bool is_minimal(const std::vector<scalar_type> R_vec,
                         const std::vector<std::vector<scalar_type>>& basis);

  template <typename scalar_type>
  static std::vector<std::vector<scalar_type>> equivalent_vectors(
      std::vector<scalar_type> R_vec, const std::vector<std::vector<scalar_type>>& basis);

  // Finds and returns the cluster vector whose distance (L2 norm) squared to input_vec is minimal.
  template <typename scalar_type>
  static std::vector<scalar_type> find_closest_cluster_vector(
      const std::vector<scalar_type>& input_vec,
      const std::vector<std::vector<scalar_type>>& cluster_vectors,
      const std::vector<std::vector<scalar_type>>& super_basis);

  // Finds and returns the cluster vector whose distance (L2 norm) squared to input_vec is minimal.
  // If the distance is larger than 'tolerance' a logic_error is thrown.
  template <typename scalar_type>
  static std::vector<scalar_type> find_closest_cluster_vector(
      const std::vector<scalar_type>& input_vec,
      const std::vector<std::vector<scalar_type>>& cluster_vectors,
      const std::vector<std::vector<scalar_type>>& super_basis, const scalar_type tolerance);
};

template <typename scalar_type>
int cluster_operations::index(const std::vector<scalar_type>& element,
                              const std::vector<std::vector<scalar_type>>& elements,
                              const CLUSTER_SHAPE shape) {
  int index = -1;

  switch (shape) {
    case BRILLOUIN_ZONE:  // the k-vectors in the brillouin zone are sorted according to
                          // math::util::isLessVector !
      index = lower_bound(elements.begin(), elements.end(), element,
                          math::util::isLessVector<scalar_type>) -
              elements.begin();
      break;

    case PARALLELLEPIPEDUM:
      index = find(elements.begin(), elements.end(), element) - elements.begin();
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  assert(index > -1 and index < elements.size());

  if (math::util::distance2(element, elements[index]) > 1.e-6) {
    std::cout << "\n\t " << __FUNCTION__ << "\t" << index << "\n";
    math::util::print(element);
    std::cout << "\n";
    math::util::print(elements[index]);
    std::cout << "\n";
    std::cout << "\n\n";
  }

  assert(math::util::distance2(element, elements[index]) < 1.e-6);

  return index;
}

// INTERNAL: This function doesn't work with scalar_type != double since origin is of type
//       std::vector<double>.
template <typename scalar_type>
int cluster_operations::origin_index(const std::vector<std::vector<scalar_type>>& elements,
                                     const CLUSTER_SHAPE shape) {
  std::vector<double> origin(elements[0].size(), 0.);

  return index(origin, elements, shape);
}

template <typename scalar_type>
std::vector<scalar_type> cluster_operations::translate_inside_cluster(
    const std::vector<scalar_type>& r, const std::vector<std::vector<scalar_type>>& basis) {
  int DIMENSION = r.size();

  std::vector<scalar_type> r_affine = math::util::coordinates(r, basis);

  for (size_t d = 0; d < r.size(); d++) {
    while (r_affine[d] < -1.e-6) {
      r_affine[d] += 1.;
    }

    while (r_affine[d] > 1 - 1.e-6) {
      r_affine[d] -= 1.;
    }
  }

  std::vector<scalar_type> r_vec(r.size(), 0.);

  for (int d1 = 0; d1 < DIMENSION; ++d1) {
    for (int d0 = 0; d0 < DIMENSION; ++d0) {
      r_vec[d0] += basis[d1][d0] * r_affine[d1];
    }
  }

  return r_vec;
}

template <typename scalar_type>
bool cluster_operations::test_translate_inside_cluster(
    const std::vector<std::vector<scalar_type>>& elements,
    const std::vector<std::vector<scalar_type>>& basis) {
  static bool passed_test = false;

  if (!passed_test) {
    std::vector<scalar_type> k1, k2;

    for (size_t l = 0; l < elements.size(); l++) {
      k1 = elements[l];
      k2 = translate_inside_cluster(k1, basis);

      if (math::util::distance2(k1, k2) > 1.e-6) {
        throw std::logic_error(__FUNCTION__);
      }
    }

    passed_test = true;
  }

  return passed_test;
}

template <typename scalar_type>
scalar_type cluster_operations::minimal_distance(std::vector<scalar_type> vec_0,
                                                 std::vector<scalar_type> vec_1,
                                                 const std::vector<std::vector<scalar_type>>& basis) {
  int DIMENSION = vec_0.size();

  vec_0 = translate_inside_cluster(vec_0, basis);
  vec_1 = translate_inside_cluster(vec_1, basis);

  scalar_type MIN_DISTANCE = math::util::distance2(vec_0, vec_1);

  switch (DIMENSION) {
    case 1:
      for (int l0 = -1; l0 <= 1; l0++) {
        std::vector<scalar_type> vec = vec_0;

        for (int d = 0; d < DIMENSION; d++) {
          vec[d] += (l0 * basis[0][d]);
        }

        scalar_type distance = math::util::distance2(vec, vec_1);

        if (distance < MIN_DISTANCE) {
          MIN_DISTANCE = distance;
        }
      }
      break;

    case 2:
      for (int l0 = -1; l0 <= 1; l0++) {
        for (int l1 = -1; l1 <= 1; l1++) {
          std::vector<scalar_type> vec = vec_0;

          for (int d = 0; d < DIMENSION; d++) {
            vec[d] += (l0 * basis[0][d] + l1 * basis[1][d]);
          }

          scalar_type distance = math::util::distance2(vec, vec_1);

          if (distance < MIN_DISTANCE) {
            MIN_DISTANCE = distance;
          }
        }
      }
      break;

    case 3:
      for (int l0 = -1; l0 <= 1; l0++) {
        for (int l1 = -1; l1 <= 1; l1++) {
          for (int l2 = -1; l2 <= 1; l2++) {
            std::vector<scalar_type> vec = vec_0;

            for (int d = 0; d < DIMENSION; d++) {
              vec[d] += (l0 * basis[0][d] + l1 * basis[1][d] + l2 * basis[2][d]);
            }

            scalar_type distance = math::util::distance2(vec, vec_1);

            if (distance < MIN_DISTANCE) {
              MIN_DISTANCE = distance;
            }
          }
        }
      }
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  return MIN_DISTANCE;
}

template <typename scalar_type>
bool cluster_operations::is_minimal(const std::vector<scalar_type> R_vec,
                                    const std::vector<std::vector<scalar_type>>& basis) {
  int DIMENSION = R_vec.size();

  std::vector<scalar_type> origin(DIMENSION, 0.);

  scalar_type MIN_DISTANCE = minimal_distance(origin, R_vec, basis);

  bool minimal = math::util::l2Norm2(R_vec) > (MIN_DISTANCE + 1.e-6) ? false : true;

  return minimal;
}

template <typename scalar_type>
std::vector<std::vector<scalar_type>> cluster_operations::equivalent_vectors(
    std::vector<scalar_type> R_vec, const std::vector<std::vector<scalar_type>>& basis) {
  const static scalar_type EPS = 1.e-3;
  const static scalar_type ONE_PLUS_EPS = 1. + EPS;
  const static scalar_type ONE_MIN_EPS = 1. - EPS;

  int DIMENSION = R_vec.size();

  R_vec = translate_inside_cluster(R_vec, basis);

  std::vector<scalar_type> origin(DIMENSION, 0.);

  scalar_type MIN_DISTANCE = minimal_distance(origin, R_vec, basis);

  // bool IS_MINIMAL = L2_norm(R_vec)>(MIN_DISTANCE+1.e-6)? false : true;

  std::vector<std::vector<scalar_type>> r_min;

  switch (DIMENSION) {
    case 1: {
      for (int l0 = -2; l0 <= 2; l0++) {
        std::vector<scalar_type> vec = R_vec;

        for (int d = 0; d < DIMENSION; d++)
          vec[d] += l0 * basis[0][d];

        scalar_type distance = std::sqrt(math::util::l2Norm2(vec));

        if (distance > sqrt(MIN_DISTANCE) * ONE_MIN_EPS - EPS &&
            distance < sqrt(MIN_DISTANCE) * ONE_PLUS_EPS + EPS)
          r_min.push_back(vec);
      }
    } break;

    case 2: {
      for (int l0 = -2; l0 <= 2; l0++) {
        for (int l1 = -2; l1 <= 2; l1++) {
          std::vector<scalar_type> vec = R_vec;

          for (int d = 0; d < DIMENSION; d++)
            vec[d] += (l0 * basis[0][d] + l1 * basis[1][d]);

          scalar_type distance = std::sqrt(math::util::l2Norm2(vec));

          if (distance > sqrt(MIN_DISTANCE) * ONE_MIN_EPS - EPS &&
              distance < sqrt(MIN_DISTANCE) * ONE_PLUS_EPS + EPS)
            r_min.push_back(vec);
        }
      }
    } break;

    case 3: {
      for (int l0 = -2; l0 <= 2; l0++) {
        for (int l1 = -2; l1 <= 2; l1++) {
          for (int l2 = -2; l2 <= 2; l2++) {
            std::vector<scalar_type> vec = R_vec;

            for (int d = 0; d < DIMENSION; d++)
              vec[d] += (l0 * basis[0][d] + l1 * basis[1][d] + l2 * basis[2][d]);

            scalar_type distance = std::sqrt(math::util::l2Norm2(vec));

            if (distance > sqrt(MIN_DISTANCE) * ONE_MIN_EPS - EPS &&
                distance < sqrt(MIN_DISTANCE) * ONE_PLUS_EPS + EPS)
              r_min.push_back(vec);
          }
        }
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  sort(r_min.begin(), r_min.end(), math::util::isSameVector<scalar_type>);
  int vec_size =
      unique(r_min.begin(), r_min.end(), math::util::isSameVector<scalar_type>) - r_min.begin();

  r_min.resize(vec_size);

  return r_min;
}

template <typename scalar_type>
std::vector<scalar_type> cluster_operations::find_closest_cluster_vector(
    const std::vector<scalar_type>& input_vec,
    const std::vector<std::vector<scalar_type>>& cluster_vectors,
    const std::vector<std::vector<scalar_type>>& super_basis) {
  if (cluster_vectors.size() == 0) {
    throw std::logic_error(__FUNCTION__);
  }

  double min_distance = minimal_distance(input_vec, cluster_vectors[0], super_basis);
  int min_index = 0;

  for (int i = 1; i < cluster_vectors.size(); ++i) {
    double distance = minimal_distance(input_vec, cluster_vectors[i], super_basis);
    if (distance < min_distance) {
      min_distance = distance;
      min_index = i;
    }
  }

  return cluster_vectors[min_index];
}

template <typename scalar_type>
std::vector<scalar_type> cluster_operations::find_closest_cluster_vector(
    const std::vector<scalar_type>& input_vec,
    const std::vector<std::vector<scalar_type>>& cluster_vectors,
    const std::vector<std::vector<scalar_type>>& super_basis, const scalar_type tolerance) {
  std::vector<scalar_type> result_vec =
      find_closest_cluster_vector(input_vec, cluster_vectors, super_basis);

  if (std::sqrt(minimal_distance(input_vec, result_vec, super_basis)) > tolerance) {
    throw std::logic_error(__FUNCTION__);
  }
  return result_vec;
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_OPERATIONS_HPP
