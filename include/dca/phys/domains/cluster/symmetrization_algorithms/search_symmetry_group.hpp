// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// The operation is lattice compatible if there exists a integer matrix I, such that
//     O*V = V*I .
// The operation is lattice compatible if there exists a integer matrix I and a permutation matrix
// P, such that
//     O*a = a*P .

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_SEARCH_SYMMETRY_GROUP_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_SEARCH_SYMMETRY_GROUP_HPP

#include <stdexcept>
#include <vector>

#include "dca/linalg/linalg.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/point_group_symmetry_domain.hpp"
#include "dca/phys/domains/quantum/symmetry_group_level.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <class base_cluster_type, class point_group,
          domains::symmetry_group_level_type symmetry_group_level, int INDEX>
class search_symmetry_group {
private:
  typedef typename dca::util::TypeAt<INDEX, point_group>::type symmetry_type;
  typedef typename symmetry_type::base_type group_action_type;

  const static int DIMENSION = base_cluster_type::DIMENSION;

  typedef domains::electron_band_domain b_dmn_t;

  //   typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
  //   typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

  typedef typename base_cluster_type::r_cluster_type r_cluster_type;
  typedef typename base_cluster_type::k_cluster_type k_cluster_type;

  typedef domains::point_group_symmetry_domain<symmetry_group_level, base_cluster_type> sym_dmn_t;

  typedef typename sym_dmn_t::element_type symmetry_element_type;

private:
  template <domains::symmetry_group_level_type my_symmetry_group_level>
  struct intialize_T_matrices {
    static void execute(double* T, double* T_inv) {
      switch (my_symmetry_group_level) {
        case domains::UNIT_CELL:
          for (int i = 0; i < DIMENSION; ++i)
            for (int j = 0; j < DIMENSION; ++j)
              T[j + DIMENSION * i] = r_cluster_type::get_basis_vectors()[i][j];
          break;

        case domains::SUPER_CELL:
          for (int i = 0; i < DIMENSION; ++i)
            for (int j = 0; j < DIMENSION; ++j)
              T[j + DIMENSION * i] = r_cluster_type::get_super_basis_vectors()[i][j];
          break;

        default:
          throw std::logic_error(__FUNCTION__);
      }

      dca::linalg::lapack::lacpy("A", DIMENSION, DIMENSION, T, DIMENSION, T_inv, DIMENSION);
      dca::linalg::lapack::inverse(DIMENSION, T_inv, DIMENSION);
    }
  };

  static void back_inside_cluster(double* T_inv, double* v, double* w) {
    for (int i = 0; i < DIMENSION; ++i)
      w[i] = 0.;

    for (int i = 0; i < DIMENSION; i++)
      for (int j = 0; j < DIMENSION; j++)
        w[i] += T_inv[i + DIMENSION * j] * v[j];

    for (int i = 0; i < DIMENSION; ++i) {
      while (w[i] < -1.e-6)
        w[i] += 1.;

      while (w[i] > 1 - 1.e-6)
        w[i] -= 1.;
    }
  }

  static bool is_duplicate() {
    if (sym_dmn_t::get_size() == 0)
      return false;
    else {
      bool is_a_duplicate = false;

      for (int l = 0; l < sym_dmn_t::get_size(); l++) {
        bool is_this_a_duplicate = true;

        for (int i = 0; i < DIMENSION; i++)
          for (int j = 0; j < DIMENSION; j++)
            if (std::fabs(sym_dmn_t::get_elements()[l].O[i + j * DIMENSION] -
                          symmetry_type::matrix()[i + j * DIMENSION]) > 1.e-6)
              is_this_a_duplicate = false;

        if (is_this_a_duplicate)
          is_a_duplicate = true;
      }

      return is_a_duplicate;
    }
  }

  static bool is_lattice_compatible(double* T, double* T_inv) {
    double* permutation = new double[DIMENSION * DIMENSION];

    for (int i = 0; i < DIMENSION; i++)
      for (int j = 0; j < DIMENSION; j++)
        permutation[i + DIMENSION * j] = 0.;

    for (int i = 0; i < DIMENSION; i++)
      for (int j = 0; j < DIMENSION; j++)
        for (int l = 0; l < DIMENSION; l++)
          permutation[i + DIMENSION * j] +=
              symmetry_type::matrix()[i + DIMENSION * l] * T[l + DIMENSION * j];

    bool is_compatible = true;

    std::vector<double> v(DIMENSION, 0);
    std::vector<double> v_c(DIMENSION, 0);

    for (int i = 0; i < DIMENSION; i++) {
      for (int j = 0; j < DIMENSION; j++)
        v[j] = permutation[j + DIMENSION * i];

      back_inside_cluster(&T_inv[0], &v[0], &v_c[0]);

      if (math::util::l2Norm2(v_c) > 1.e-6)
        is_compatible = false;
    }

    delete[] permutation;

    return is_compatible;
  }

  static bool is_unit_cell_compatible(double* /*T*/, double* T_inv) {
      // NOTE: This function does not work as intended. Presumbaly, it's purpose is to check if the unit call maps onto itself under
      // a certain symmetry operation. However, it will always return 0 if there is jsut one site in unit cell at the origin,
      // which will always map onto itself. In this case, it will return true, even though the other sites may not map onto themselves.
      // Also, if the unit cell does not contain a site at the origin, even for the identity operation, it will return false,
      // presumbaly becaise there is something wrong with the back_inside_cluster function.

      int Na = b_dmn_t::get_size();

    // std::cout << "avec0 BEFORE!!!!!: \n";
    // for (int j = 0; j < DIMENSION; j++)        std::cout << b_dmn_t::get_elements()[0].a_vec[j] << "\n";



    double* permutation = new double[Na * DIMENSION];

    for (int j = 0; j < Na; j++)
      for (int i = 0; i < DIMENSION; i++)
        permutation[i + DIMENSION * j] = 0;  // b_dmn_t::get_elements()[i].a_vec[j];

    for (int j = 0; j < Na; j++)
      for (int i = 0; i < DIMENSION; i++)
        for (int l = 0; l < DIMENSION; l++) {
          permutation[i + DIMENSION * j] += symmetry_type::matrix()[i + DIMENSION * l] *
                                            b_dmn_t::get_elements()[j].a_vec[l];  // T[l+DIMENSION*j];
         // if (j==0) std::cout << "avec INSIDE: " << b_dmn_t::get_elements()[j].a_vec[l] << "\n";
         // if (j==0) std::cout << "permutation INSIDE: " << permutation[i + DIMENSION * j] << "\n";
        }

    bool is_compatible = false;

    std::vector<double> v(DIMENSION, 0);
    std::vector<double> v_c(DIMENSION, 0);

    for (int l = 0; l < Na; l++) {
      for (int j = 0; j < DIMENSION; j++)
        v[j] = permutation[j + DIMENSION * l];

      // if (l==5) {
      // std::cout << "v before back_inside_cluster: \n";
      // for (int j = 0; j < DIMENSION; j++)        std::cout << "\t" << v[j] << "\n";
      //   }

      // v_c is affine coordinates
      back_inside_cluster(&T_inv[0], &v[0], &v_c[0]);


      for (int j = 0; j < DIMENSION; j++)
        v[j] = 0;  // permutation[j+DIMENSION*i];

      for (int i = 0; i < DIMENSION; i++)
        for (int j = 0; j < DIMENSION; j++)
          v[i] += symmetry_type::matrix()[i + DIMENSION * j] * v_c[j];

      // if (l==5) {
      // std::cout << "vc after back_inside_cluster: \n";
      // for (int j = 0; j < DIMENSION; j++)        std::cout << "\t" << v_c[j] << "\n";
      // std::cout << "v after back_inside_cluster: \n";
      // for (int j = 0; j < DIMENSION; j++)        std::cout << "\t" << v[j] << "\n";
      //   }

    // if (l==0) {
    //     std::cout << "v : \n";
    //     for (int j = 0; j < DIMENSION; j++)        std::cout << "\t" << v[j] << "\n";
    //     std::cout << "vc : \n";
    //     for (int j = 0; j < DIMENSION; j++)        std::cout << "\t" << v_c[j] << "\n";
    //     std::cout << "avec0: \n";
    //     for (int j = 0; j < DIMENSION; j++)        std::cout << b_dmn_t::get_elements()[0].a_vec[j] << "\n";
    // }

      for (int j = 0; j < Na; j++) {
        // std::cout << "l, j, distance "<< l << "," << j << "," << math::util::distance2(v, b_dmn_t::get_elements()[j].a_vec) << "\n";
        if (math::util::distance2(v, b_dmn_t::get_elements()[j].a_vec) < 1.e-6 and
            b_dmn_t::get_elements()[l].flavor == b_dmn_t::get_elements()[j].flavor) {
          is_compatible = true;
          // std::cout << "is_compatible for l, j, distance "<< l << "," << j << "," << math::util::distance2(v, b_dmn_t::get_elements()[j].a_vec) << "\n";
            }
      }
    }

    delete[] permutation;

    return is_compatible;
  }

public:
  static void execute() {
    double* T = new double[DIMENSION * DIMENSION];
    double* T_inv = new double[DIMENSION * DIMENSION];

    double* T_cell = new double[DIMENSION * DIMENSION];
    double* T_cell_inv = new double[DIMENSION * DIMENSION];

    intialize_T_matrices<domains::UNIT_CELL>::execute(T_cell, T_cell_inv);
    intialize_T_matrices<symmetry_group_level>::execute(T, T_inv);

    bool is_a_duplicate = is_duplicate();

    // bool unit_cell_compatible = is_unit_cell_compatible(T_cell, T_cell_inv);
    // NOTE: this function does not work as intended, see the comment in the function definition. For now, we set it to true, which is a necessary condition for the symmetry operation to be a symmetry of the system, but not a sufficient one.
    bool unit_cell_compatible = true;
    bool lattice_compatible = is_lattice_compatible(T, T_inv);

    // std::cout << "is_unit_cell_compatible " << unit_cell_compatible << "\n";
    // std::cout << "is_lattice_compatible " << lattice_compatible << "\n";
    // std::cout << "is_a_duplicate " << is_a_duplicate << "\n";

    if (lattice_compatible && unit_cell_compatible && not is_a_duplicate) {
      symmetry_element_type symmetry_element(DIMENSION);

      for (int i = 0; i < DIMENSION; i++)
        for (int j = 0; j < DIMENSION; j++)
          symmetry_element.O[i + j * DIMENSION] = symmetry_type::matrix()[i + j * DIMENSION];

      for (int i = 0; i < DIMENSION; i++)
        symmetry_element.t[i] = 0.;

      sym_dmn_t::get_elements().push_back(symmetry_element);

      sym_dmn_t::get_size() += 1;
    }

    delete[] T;
    delete[] T_inv;

    delete[] T_cell;
    delete[] T_cell_inv;

    search_symmetry_group<base_cluster_type, point_group, symmetry_group_level, INDEX - 1>::execute();
  }
};

template <class cluster_type, class point_group, domains::symmetry_group_level_type symmetry_group_level>
class search_symmetry_group<cluster_type, point_group, symmetry_group_level, -1> {
public:
  static void execute() {}
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_SEARCH_SYMMETRY_GROUP_HPP
