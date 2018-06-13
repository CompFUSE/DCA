// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class initializes the cluster domain.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_INITIALIZER_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_INITIALIZER_HPP

#include <cassert>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/util/coordinate_transformation.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename cluster_type>
class cluster_domain_initializer {};

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
class cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>> {
public:
  typedef std::vector<scalar_type> element_type;

  typedef cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE> r_dmn;
  typedef cluster_domain<scalar_type, DIMENSION, NAME, MOMENTUM_SPACE, SHAPE> k_dmn;

  static void execute(scalar_type* r_basis, std::vector<int> R_basis,
                      bool init_add_and_subtract_matrices = false);

  static void execute(scalar_type* r_basis, std::vector<std::vector<int>> R_basis,
                      bool init_add_and_subtract_matrices = true);

private:
  static void allocate_data(int N);

  static void initialize_basis(scalar_type* r_basis, std::vector<int> R_basis);

  static void initialize_basis(scalar_type* r_basis, std::vector<std::vector<int>> R_basis);

  static void compute_basis();

  static void initialize_basis_vectors();

  static void initialize_elements();

  static void initialize_elements_1D(scalar_type shift);
  static void initialize_elements_2D(scalar_type shift);
  static void initialize_elements_3D(scalar_type shift);

  static void initialize_elements(std::vector<int> R_basis);

  static void initialize_elements_1D(std::vector<int> R_basis);
  static void initialize_elements_2D(std::vector<int> R_basis);
  static void initialize_elements_3D(std::vector<int> R_basis);

  static void initialize_add(scalar_type* basis, std::vector<std::vector<scalar_type>>& elements,
                             dca::linalg::Matrix<int, dca::linalg::CPU>& A);

  static void initialize_subtract(scalar_type* basis, std::vector<std::vector<scalar_type>>& elements,
                                  dca::linalg::Matrix<int, dca::linalg::CPU>& A);

  static void initialize_volume();
};

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<cluster_domain<
    scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::execute(scalar_type* r_basis,
                                                                std::vector<int> grid_size,
                                                                bool init_add_and_subtract_matrices) {
  assert(SHAPE == PARALLELLEPIPEDUM);

  {
    r_dmn::get_dimensions() = grid_size;
    k_dmn::get_dimensions() = grid_size;
  }

  {
    allocate_data(grid_size.size());

    initialize_basis(r_basis, grid_size);

    compute_basis();

    initialize_basis_vectors();
  }

  initialize_elements(grid_size);

  r_dmn::is_initialized() = true;
  k_dmn::is_initialized() = true;

  if (init_add_and_subtract_matrices) {
    //       initialize_add(grid_size, r_dmn::get_elements(), r_dmn::get_add_matrix());
    //       initialize_add(grid_size, k_dmn::get_elements(), k_dmn::get_add_matrix());

    //       initialize_subtract(grid_size, r_dmn::get_elements(), r_dmn::get_subtract_matrix());
    //       initialize_subtract(grid_size, k_dmn::get_elements(), k_dmn::get_subtract_matrix());
  }

  initialize_volume();
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<cluster_domain<
    scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::execute(scalar_type* r_basis,
                                                                std::vector<std::vector<int>> R_basis,
                                                                bool init_add_and_subtract_matrices) {
  assert(SHAPE == BRILLOUIN_ZONE);

  {
    r_dmn::get_dimensions().resize(0);
    k_dmn::get_dimensions().resize(0);
  }

  {
    allocate_data(R_basis.size());

    initialize_basis(r_basis, R_basis);

    compute_basis();

    initialize_basis_vectors();
  }

  initialize_elements();

  r_dmn::is_initialized() = true;
  k_dmn::is_initialized() = true;

  if (init_add_and_subtract_matrices) {
    initialize_add(r_dmn::get_super_basis(), r_dmn::get_elements(), r_dmn::get_add_matrix());
    initialize_add(k_dmn::get_super_basis(), k_dmn::get_elements(), k_dmn::get_add_matrix());

    initialize_subtract(r_dmn::get_super_basis(), r_dmn::get_elements(),
                        r_dmn::get_subtract_matrix());
    initialize_subtract(k_dmn::get_super_basis(), k_dmn::get_elements(),
                        k_dmn::get_subtract_matrix());
  }

  initialize_volume();
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<
    cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::allocate_data(int /*N*/) {
  // assert(N == DIMENSION);

  {
    r_dmn::get_basis() = new scalar_type[DIMENSION * DIMENSION];
    r_dmn::get_super_basis() = new scalar_type[DIMENSION * DIMENSION];

    r_dmn::get_inverse_basis() = new scalar_type[DIMENSION * DIMENSION];
    r_dmn::get_inverse_super_basis() = new scalar_type[DIMENSION * DIMENSION];
  }

  {
    k_dmn::get_basis() = new scalar_type[DIMENSION * DIMENSION];
    k_dmn::get_super_basis() = new scalar_type[DIMENSION * DIMENSION];

    k_dmn::get_inverse_basis() = new scalar_type[DIMENSION * DIMENSION];
    k_dmn::get_inverse_super_basis() = new scalar_type[DIMENSION * DIMENSION];
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<cluster_domain<
    scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_basis(scalar_type* r_basis,
                                                                         std::vector<int> R_basis) {
  {
    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        r_dmn::get_basis()[d0 + d1 * DIMENSION] = r_basis[d0 + d1 * DIMENSION];
  }

  {
    for (int d1 = 0; d1 < DIMENSION; d1++)
      for (int d0 = 0; d0 < DIMENSION; d0++)
        r_dmn::get_super_basis()[d0 + d1 * DIMENSION] = r_basis[d0 + d1 * DIMENSION] * R_basis[d1];
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::
    initialize_basis(scalar_type* r_basis, std::vector<std::vector<int>> R_basis) {
  {
    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        r_dmn::get_basis()[d0 + d1 * DIMENSION] = r_basis[d0 + d1 * DIMENSION];
  }

  {
    for (int d0 = 0; d0 < DIMENSION; d0++) {
      for (int d1 = 0; d1 < DIMENSION; d1++) {
        r_dmn::get_super_basis()[d0 + d1 * DIMENSION] = 0;

        for (int d2 = 0; d2 < DIMENSION; d2++)
          r_dmn::get_super_basis()[d0 + d1 * DIMENSION] +=
              r_dmn::get_basis()[d0 + d2 * DIMENSION] * R_basis[d1][d2];
      }
    }
  }
}

/*!
 *   The convention is that the basis and superbasis vectors are stored in the columns!
 *
 *   The convention is that the inverse basis (and inverse superbasis) is defined as the inverse of
 * the basis (superbasis) matrix!
 */
template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::compute_basis() {
  {
    dca::linalg::Matrix<scalar_type, dca::linalg::CPU> A("A", DIMENSION);

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        A(d0, d1) = r_dmn::get_basis()[d0 + d1 * DIMENSION];

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        k_dmn::get_inverse_super_basis()[d0 + d1 * DIMENSION] = A(d1, d0) / (2 * M_PI);

    dca::linalg::matrixop::inverse(A);

    //     for(int d0=0; d0<DIMENSION; d0++)
    //       for(int d1=0; d1<DIMENSION; d1++)
    //  r_dmn::get_inverse_basis()[d0+d1*DIMENSION] = A(d0,d1);

    //     for(int d0=0; d0<DIMENSION; d0++)
    //       for(int d1=0; d1<DIMENSION; d1++)
    //  k_dmn::get_super_basis()[d0+d1*DIMENSION] = A(d0,d1)*(2*M_PI);

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        r_dmn::get_inverse_basis()[d0 + d1 * DIMENSION] = A(d0, d1);

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        k_dmn::get_super_basis()[d0 + d1 * DIMENSION] = A(d1, d0) * (2 * M_PI);

    if (true)  // test
    {
      for (int d0 = 0; d0 < DIMENSION; d0++) {
        for (int d1 = 0; d1 < DIMENSION; d1++) {
          scalar_type result = 0;
          for (int d2 = 0; d2 < DIMENSION; d2++)
            result += r_dmn::get_basis()[d0 + d2 * DIMENSION] *
                      r_dmn::get_inverse_basis()[d2 + d1 * DIMENSION];

          if (d0 == d1)
            result -= 1;

          if (std::abs(result) > 1.e-6)
            throw std::logic_error(__FUNCTION__);
        }
      }
    }
  }

  {
    dca::linalg::Matrix<scalar_type, dca::linalg::CPU> A("A", DIMENSION);

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        A(d0, d1) = r_dmn::get_super_basis()[d0 + d1 * DIMENSION];

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        k_dmn::get_inverse_basis()[d0 + d1 * DIMENSION] = A(d1, d0) / (2 * M_PI);

    dca::linalg::matrixop::inverse(A);

    //     for(int d0=0; d0<DIMENSION; d0++)
    //       for(int d1=0; d1<DIMENSION; d1++)
    //  r_dmn::get_inverse_super_basis()[d0+d1*DIMENSION] = A(d0,d1);

    //     for(int d0=0; d0<DIMENSION; d0++)
    //       for(int d1=0; d1<DIMENSION; d1++)
    //  k_dmn::get_basis()[d0+d1*DIMENSION] = A(d0,d1)*(2*M_PI);

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        r_dmn::get_inverse_super_basis()[d0 + d1 * DIMENSION] = A(d0, d1);

    for (int d0 = 0; d0 < DIMENSION; d0++)
      for (int d1 = 0; d1 < DIMENSION; d1++)
        k_dmn::get_basis()[d0 + d1 * DIMENSION] = A(d1, d0) * (2 * M_PI);

    if (true)  // test
    {
      for (int d0 = 0; d0 < DIMENSION; d0++) {
        for (int d1 = 0; d1 < DIMENSION; d1++) {
          scalar_type result = 0;
          for (int d2 = 0; d2 < DIMENSION; d2++)
            result += r_dmn::get_super_basis()[d0 + d2 * DIMENSION] *
                      r_dmn::get_inverse_super_basis()[d2 + d1 * DIMENSION];

          if (d0 == d1)
            result -= 1;

          if (std::abs(result) > 1.e-6)
            throw std::logic_error(__FUNCTION__);
        }
      }
    }
  }
}

/*!
 *   The convention is that the basis and superbasis vectors are stored in the columns!
 *
 *   The convention is that the inverse basis (and inverse superbasis) is defined as the inverse of
 * the basis (superbasis) matrix!
 */
template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<
    cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_basis_vectors() {
  {
    r_dmn::get_basis_vectors().resize(0);

    for (int d0 = 0; d0 < DIMENSION; d0++) {
      element_type tmp(0);
      for (int d1 = 0; d1 < DIMENSION; d1++)
        tmp.push_back(r_dmn::get_basis()[d1 + d0 * DIMENSION]);

      r_dmn::get_basis_vectors().push_back(tmp);
    }
  }

  {
    k_dmn::get_basis_vectors().resize(0);

    for (int d0 = 0; d0 < DIMENSION; d0++) {
      element_type tmp(0);
      for (int d1 = 0; d1 < DIMENSION; d1++)
        tmp.push_back(k_dmn::get_basis()[d1 + d0 * DIMENSION]);

      k_dmn::get_basis_vectors().push_back(tmp);
    }
  }

  {
    r_dmn::get_super_basis_vectors().resize(0);

    for (int d0 = 0; d0 < DIMENSION; d0++) {
      element_type tmp(0);
      for (int d1 = 0; d1 < DIMENSION; d1++)
        tmp.push_back(r_dmn::get_super_basis()[d1 + d0 * DIMENSION]);

      r_dmn::get_super_basis_vectors().push_back(tmp);
    }
  }

  {
    k_dmn::get_super_basis_vectors().resize(0);

    for (int d0 = 0; d0 < DIMENSION; d0++) {
      element_type tmp(0);
      for (int d1 = 0; d1 < DIMENSION; d1++)
        tmp.push_back(k_dmn::get_super_basis()[d1 + d0 * DIMENSION]);

      k_dmn::get_super_basis_vectors().push_back(tmp);
    }
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_elements() {
  const static scalar_type shift = 0.0;

  switch (DIMENSION) {
    case 1:
      initialize_elements_1D(shift);
      break;

    case 2:
      initialize_elements_2D(shift);
      break;

    case 3:
      initialize_elements_3D(shift);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  {
    r_dmn::get_size() = r_dmn::get_elements().size();
    k_dmn::get_size() = k_dmn::get_elements().size();

    sort(r_dmn::get_elements().begin(), r_dmn::get_elements().end(),
         math::util::isLessVector<scalar_type>);
    sort(k_dmn::get_elements().begin(), k_dmn::get_elements().end(),
         math::util::isLessVector<scalar_type>);
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<cluster_domain<
    scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_elements_1D(scalar_type shift) {
  {
    scalar_type t[1];
    scalar_type c[1];

    k_dmn::get_elements().resize(0);
    for (int d0 = -100; d0 < 100; d0++) {
      t[0] = d0 * k_dmn::get_basis()[0];

      c[0] = t[0] / k_dmn::get_super_basis()[0];

      if (c[0] > -1.e-6 and c[0] < 1 - 1.e-6) {
        element_type tmp(1, 0);

        tmp[0] = t[0];

        k_dmn::get_elements().push_back(tmp);
      }
    }
  }

  {
    scalar_type t[1];
    scalar_type c[1];

    r_dmn::get_elements().resize(0);
    for (int d0 = -100; d0 < 100; d0++) {
      t[0] = d0 * r_dmn::get_basis()[0];

      c[0] = t[0] / r_dmn::get_super_basis()[0];

      if (c[0] > shift - 1.e-6 and c[0] < shift + 1. - 1.e-6) {
        element_type tmp(1, 0);

        tmp[0] = t[0];

        r_dmn::get_elements().push_back(tmp);
      }
    }
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<cluster_domain<
    scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_elements_2D(scalar_type shift) {
  math::util::coordinate_transformation<scalar_type> coordinate_trafo(2);

  {
    coordinate_trafo.set_basis(k_dmn::get_super_basis());

    scalar_type t[2];
    scalar_type c[2];

    k_dmn::get_elements().resize(0);
    for (int d0 = -100; d0 < 100; d0++) {
      for (int d1 = -100; d1 < 100; d1++) {
        t[0] = d0 * k_dmn::get_basis()[0] + d1 * k_dmn::get_basis()[2];
        t[1] = d0 * k_dmn::get_basis()[1] + d1 * k_dmn::get_basis()[3];

        coordinate_trafo.execute(t, c);

        if (c[0] > shift - 1.e-6 and c[0] < shift + 1. - 1.e-6 and c[1] > shift - 1.e-6 and
            c[1] < shift + 1. - 1.e-6) {
          element_type tmp(2, 0);

          tmp[0] = t[0];
          tmp[1] = t[1];

          k_dmn::get_elements().push_back(tmp);
        }
      }
    }
  }

  {
    coordinate_trafo.set_basis(r_dmn::get_super_basis());

    scalar_type t[2];
    scalar_type c[2];

    r_dmn::get_elements().resize(0);
    for (int d0 = -100; d0 < 100; d0++) {
      for (int d1 = -100; d1 < 100; d1++) {
        t[0] = d0 * r_dmn::get_basis()[0] + d1 * r_dmn::get_basis()[2];
        t[1] = d0 * r_dmn::get_basis()[1] + d1 * r_dmn::get_basis()[3];

        coordinate_trafo.execute(t, c);

        if (c[0] > shift - 1.e-6 and c[0] < shift + 1. - 1.e-6 and c[1] > shift - 1.e-6 and
            c[1] < shift + 1. - 1.e-6) {
          element_type tmp(2, 0);

          tmp[0] = t[0];
          tmp[1] = t[1];

          r_dmn::get_elements().push_back(tmp);
        }
      }
    }
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<cluster_domain<
    scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_elements_3D(scalar_type shift) {
  math::util::coordinate_transformation<scalar_type> coordinate_trafo(3);

  {
    coordinate_trafo.set_basis(k_dmn::get_super_basis());

    scalar_type t[3];
    scalar_type c[3];

    k_dmn::get_elements().resize(0);
    for (int d0 = -100; d0 < 100; d0++) {
      for (int d1 = -100; d1 < 100; d1++) {
        for (int d2 = -100; d2 < 100; d2++) {
          t[0] = d0 * k_dmn::get_basis()[0] + d1 * k_dmn::get_basis()[3] + d2 * k_dmn::get_basis()[6];
          t[1] = d0 * k_dmn::get_basis()[1] + d1 * k_dmn::get_basis()[4] + d2 * k_dmn::get_basis()[7];
          t[2] = d0 * k_dmn::get_basis()[2] + d1 * k_dmn::get_basis()[5] + d2 * k_dmn::get_basis()[8];

          coordinate_trafo.execute(t, c);

          if (c[0] > shift - 1.e-6 and c[0] < shift + 1. - 1.e-6 and c[1] > shift - 1.e-6 and
              c[1] < shift + 1. - 1.e-6 and c[2] > shift - 1.e-6 and c[2] < shift + 1. - 1.e-6) {
            element_type tmp(3, 0);

            tmp[0] = t[0];
            tmp[1] = t[1];
            tmp[2] = t[2];

            k_dmn::get_elements().push_back(tmp);
          }
        }
      }
    }
  }

  {
    coordinate_trafo.set_basis(r_dmn::get_super_basis());

    scalar_type t[3];
    scalar_type c[3];

    r_dmn::get_elements().resize(0);
    for (int d0 = -100; d0 < 100; d0++) {
      for (int d1 = -100; d1 < 100; d1++) {
        for (int d2 = -100; d2 < 100; d2++) {
          t[0] = d0 * r_dmn::get_basis()[0] + d1 * r_dmn::get_basis()[3] + d2 * r_dmn::get_basis()[6];
          t[1] = d0 * r_dmn::get_basis()[1] + d1 * r_dmn::get_basis()[4] + d2 * r_dmn::get_basis()[7];
          t[2] = d0 * r_dmn::get_basis()[2] + d1 * r_dmn::get_basis()[5] + d2 * r_dmn::get_basis()[8];

          coordinate_trafo.execute(t, c);

          if (c[0] > shift - 1.e-6 and c[0] < shift + 1. - 1.e-6 and c[1] > shift - 1.e-6 and
              c[1] < shift + 1. - 1.e-6 and c[2] > shift - 1.e-6 and c[2] < shift + 1. - 1.e-6) {
            element_type tmp(3, 0);

            tmp[0] = t[0];
            tmp[1] = t[1];
            tmp[2] = t[2];

            r_dmn::get_elements().push_back(tmp);
          }
        }
      }
    }
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<func::dmn_0<cluster_domain<
    scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_elements(std::vector<int> grid_size) {
  switch (DIMENSION) {
    case 1:
      initialize_elements_1D(grid_size);
      break;

    case 2:
      initialize_elements_2D(grid_size);
      break;

    case 3:
      initialize_elements_3D(grid_size);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  {
    r_dmn::get_size() = r_dmn::get_elements().size();
    k_dmn::get_size() = k_dmn::get_elements().size();
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::
    initialize_elements_1D(std::vector<int> grid_size) {
  std::vector<std::vector<double>>& r_basis = r_dmn::get_basis_vectors();
  std::vector<std::vector<double>>& k_basis = k_dmn::get_basis_vectors();

  for (int j = 0; j < grid_size[0]; j++) {
    std::vector<double> r_vec(1, 0);
    r_vec[0] = j * r_basis[0][0];

    r_dmn::get_elements().push_back(r_vec);

    std::vector<double> k_vec(1, 0);
    k_vec[0] = j * k_basis[0][0];

    k_dmn::get_elements().push_back(k_vec);
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::
    initialize_elements_2D(std::vector<int> grid_size) {
  std::vector<std::vector<double>>& r_basis = r_dmn::get_basis_vectors();
  std::vector<std::vector<double>>& k_basis = k_dmn::get_basis_vectors();

  for (int i = 0; i < grid_size[1]; i++) {
    for (int j = 0; j < grid_size[0]; j++) {
      std::vector<double> r_vec(2, 0);
      r_vec[0] = j * r_basis[0][0] + i * r_basis[1][0];
      r_vec[1] = j * r_basis[0][1] + i * r_basis[1][1];

      r_dmn::get_elements().push_back(r_vec);

      std::vector<double> k_vec(2, 0);
      k_vec[0] = j * k_basis[0][0] + i * k_basis[1][0];
      k_vec[1] = j * k_basis[0][1] + i * k_basis[1][1];

      k_dmn::get_elements().push_back(k_vec);
    }
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::
    initialize_elements_3D(std::vector<int> grid_size) {
  std::vector<std::vector<double>>& r_basis = r_dmn::get_basis_vectors();
  std::vector<std::vector<double>>& k_basis = k_dmn::get_basis_vectors();

  for (int l = 0; l < grid_size[2]; l++) {
    for (int i = 0; i < grid_size[1]; i++) {
      for (int j = 0; j < grid_size[0]; j++) {
        std::vector<double> r_vec(3, 0);
        r_vec[0] = j * r_basis[0][0] + i * r_basis[1][0] + l * r_basis[2][0];
        r_vec[1] = j * r_basis[0][1] + i * r_basis[1][1] + l * r_basis[2][1];
        r_vec[2] = j * r_basis[0][2] + i * r_basis[1][2] + l * r_basis[2][2];

        r_dmn::get_elements().push_back(r_vec);

        std::vector<double> k_vec(3, 0);
        k_vec[0] = j * k_basis[0][0] + i * k_basis[1][0] + l * k_basis[2][0];
        k_vec[1] = j * k_basis[0][1] + i * k_basis[1][1] + l * k_basis[2][1];
        k_vec[2] = j * k_basis[0][2] + i * k_basis[1][2] + l * k_basis[2][2];

        k_dmn::get_elements().push_back(k_vec);
      }
    }
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::
    initialize_add(scalar_type* basis, std::vector<std::vector<scalar_type>>& elements,
                   dca::linalg::Matrix<int, dca::linalg::CPU>& A) {
  assert(SHAPE == BRILLOUIN_ZONE);

  std::vector<std::vector<scalar_type>> basis_vecs(0);

  for (int i = 0; i < DIMENSION; i++) {
    std::vector<scalar_type> b_vec;
    for (int j = 0; j < DIMENSION; j++)
      b_vec.push_back(basis[j + i * DIMENSION]);

    basis_vecs.push_back(b_vec);
  }

  A.resizeNoCopy(elements.size());

  for (int i = 0; i < elements.size(); i++) {
    for (int j = 0; j < elements.size(); j++) {
      A(i, j) = -1;

      std::vector<scalar_type>& x_i = elements[i];
      std::vector<scalar_type>& x_j = elements[j];

      std::vector<scalar_type> x_i_plus_x_j = x_i;
      for (int d = 0; d < x_i_plus_x_j.size(); d++)
        x_i_plus_x_j[d] += x_j[d];

      x_i_plus_x_j = cluster_operations::translate_inside_cluster(x_i_plus_x_j, basis_vecs);

      A(i, j) = cluster_operations::index(x_i_plus_x_j, elements, SHAPE);

      if (A(i, j) == -1)
        throw std::logic_error(__FUNCTION__);
    }
  }

  // A.print();
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::
    initialize_subtract(scalar_type* basis, std::vector<std::vector<scalar_type>>& elements,
                        dca::linalg::Matrix<int, dca::linalg::CPU>& A) {
  assert(SHAPE == BRILLOUIN_ZONE);

  std::vector<std::vector<scalar_type>> basis_vecs(0);

  for (int i = 0; i < DIMENSION; i++) {
    std::vector<scalar_type> b_vec;
    for (int j = 0; j < DIMENSION; j++)
      b_vec.push_back(basis[j + i * DIMENSION]);

    basis_vecs.push_back(b_vec);
  }

  A.resizeNoCopy(elements.size());

  for (int i = 0; i < elements.size(); i++) {
    for (int j = 0; j < elements.size(); j++) {
      A(i, j) = -1;

      std::vector<scalar_type>& x_i = elements[i];
      std::vector<scalar_type>& x_j = elements[j];

      std::vector<scalar_type> x_j_min_x_i = x_j;
      for (int d = 0; d < DIMENSION; d++)
        x_j_min_x_i[d] -= x_i[d];

      x_j_min_x_i = cluster_operations::translate_inside_cluster(x_j_min_x_i, basis_vecs);

      A(i, j) = cluster_operations::index(x_j_min_x_i, elements, SHAPE);

      if (A(i, j) == -1)
        throw std::logic_error(__FUNCTION__);
    }
  }
}

template <typename scalar_type, int DIMENSION, CLUSTER_NAMES NAME, CLUSTER_SHAPE SHAPE>
void cluster_domain_initializer<
    func::dmn_0<cluster_domain<scalar_type, DIMENSION, NAME, REAL_SPACE, SHAPE>>>::initialize_volume() {
  switch (DIMENSION) {
    case 1:
      r_dmn::get_volume() = r_dmn::get_super_basis_vectors()[0][0];
      k_dmn::get_volume() = k_dmn::get_super_basis_vectors()[0][0];
      break;

    case 2:
      r_dmn::get_volume() =
          math::util::area(r_dmn::get_super_basis_vectors()[0], r_dmn::get_super_basis_vectors()[1]);

      k_dmn::get_volume() =
          math::util::area(k_dmn::get_super_basis_vectors()[0], k_dmn::get_super_basis_vectors()[1]);
      break;

    case 3:
      r_dmn::get_volume() = math::util::volume(r_dmn::get_super_basis_vectors()[0],
                                               r_dmn::get_super_basis_vectors()[1],
                                               r_dmn::get_super_basis_vectors()[2]);

      k_dmn::get_volume() = math::util::volume(k_dmn::get_super_basis_vectors()[0],
                                               k_dmn::get_super_basis_vectors()[1],
                                               k_dmn::get_super_basis_vectors()[2]);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_INITIALIZER_HPP
