// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the tetrahedron mesh initializer object that is templated on the dimension of
// the tetrahedron mesh and the momentum space cluster type.

#ifndef DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_HPP
#define DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_HPP

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "dca/math/geometry/tetrahedron_mesh/facet.hpp"
#include "dca/math/geometry/tetrahedron_mesh/simplex.hpp"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron.hpp"
#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

//
// Empty class template.
//
template <int dimension, class k_cluster_type>
class tetrahedron_mesh_initializer {};

//
// Template specialization for the initialization of a 2D tetrahedron mesh.
//
template <class k_cluster_type>
class tetrahedron_mesh_initializer<2, k_cluster_type> {
public:
  const static int DIMENSION = 2;

  typedef tetrahedron<2> tetrahedron_t;
  typedef simplex<2> simplex_t;
  typedef facet<2> facet_t;
  typedef std::vector<double> vector_t;

  tetrahedron_mesh_initializer(std::vector<simplex_t>& simplices_ref,
                               std::vector<facet_t>& facets_ref, std::vector<vector_t>& mesh_ref,
                               std::vector<tetrahedron_t>& tetrahedra_ref, int N_recursion_ref);

  void execute();

private:
  void make_convex_hull();
  void find_facets();
  void make_mesh_points();

  std::vector<simplex_t>& simplices;
  std::vector<facet_t>& facets;

  std::vector<vector_t>& mesh;
  std::vector<tetrahedron_t>& tetrahedra;

  int N_recursion;
};

//
// Template specialization for the initialization of a 3D tetrahedron mesh.
//
template <class k_cluster_type>
class tetrahedron_mesh_initializer<3, k_cluster_type> {
public:
  const static int DIMENSION = 3;

  typedef tetrahedron<3> tetrahedron_t;
  typedef simplex<3> simplex_t;
  typedef facet<3> facet_t;
  typedef std::vector<double> vector_t;

  tetrahedron_mesh_initializer(std::vector<simplex_t>& simplices_ref,
                               std::vector<facet_t>& facets_ref, std::vector<vector_t>& mesh_ref,
                               std::vector<tetrahedron_t>& tetrahedra_ref, int N_recursion_ref);

  void execute();

private:
  void make_convex_hull();
  void find_facets();
  void make_mesh_points();

  std::vector<simplex_t>& simplices;
  std::vector<facet_t>& facets;

  std::vector<vector_t>& mesh;
  std::vector<tetrahedron_t>& tetrahedra;

  int N_recursion;
};

//
// Definition of member functions of the 2D tetrahedron mesh initializer.
//
template <class k_cluster_type>
tetrahedron_mesh_initializer<2, k_cluster_type>::tetrahedron_mesh_initializer(
    std::vector<simplex_t>& simplices_ref, std::vector<facet_t>& facets_ref,
    std::vector<vector_t>& mesh_ref, std::vector<tetrahedron_t>& tetrahedra_ref, int N_recursion_ref)
    : simplices(simplices_ref),
      facets(facets_ref),
      mesh(mesh_ref),
      tetrahedra(tetrahedra_ref),
      N_recursion(N_recursion_ref) {}

template <class k_cluster_type>
void tetrahedron_mesh_initializer<2, k_cluster_type>::execute() {
  make_convex_hull();
  find_facets();
  make_mesh_points();
}

template <class k_cluster_type>
void tetrahedron_mesh_initializer<2, k_cluster_type>::make_convex_hull() {
  std::vector<double>& b0 = k_cluster_type::get_basis_vectors()[0];
  std::vector<double>& b1 = k_cluster_type::get_basis_vectors()[1];

  std::vector<std::vector<double>> B_collection(0);
  {
    std::vector<double> B(2, 0);

    for (int t0 = -1; t0 <= 1; t0++) {
      for (int t1 = -1; t1 <= 1; t1++) {
        if (t0 != 0 || t1 != 0) {
          B[0] = t0 * b0[0] + t1 * b1[0];
          B[1] = t0 * b0[1] + t1 * b1[1];

          B_collection.push_back(B);
        }
      }
    }
  }

  double* A = new double[2 * 2];
  double* B = new double[2];

  for (std::size_t l0 = 0; l0 < B_collection.size(); l0++) {
    for (std::size_t l1 = 0; l1 < B_collection.size(); l1++) {
      A[0 + 2 * 0] = B_collection[l0][0];
      A[0 + 2 * 1] = B_collection[l0][1];
      A[1 + 2 * 0] = B_collection[l1][0];
      A[1 + 2 * 1] = B_collection[l1][1];

      B[0] = B_collection[l0][0] * B_collection[l0][0] / 2. +
             B_collection[l0][1] * B_collection[l0][1] / 2.;
      B[1] = B_collection[l1][0] * B_collection[l1][0] / 2. +
             B_collection[l1][1] * B_collection[l1][1] / 2.;

      try {
        linalg::lapack::solve(2, A, 2, B);  // overwrites A with its LU factorization.
      }
      catch (linalg::lapack::util::LapackException& e) {
        if (e.info() < 0)
          // Argument error: re-throw.
          throw;
        // A is singular.
        continue;
      }

      double det_A = A[0] * A[3];

      if (std::abs(det_A) > 1.e-6) {
        simplex<DIMENSION> s;
        s.k_vec.resize(2, 0);

        s.k_vec[0] = B[0];
        s.k_vec[1] = B[1];

        simplices.push_back(s);
      }
    }
  }

  delete[] A;
  delete[] B;

  {
    std::vector<double> K(2, 0);

    for (std::size_t B_ind = 0; B_ind < B_collection.size(); B_ind++) {
      for (std::size_t s_ind = 0; s_ind < simplices.size(); s_ind++) {
        double diff_k_K = math::util::distance2(simplices[s_ind].k_vec, K);
        double diff_k_B = math::util::distance2(simplices[s_ind].k_vec, B_collection[B_ind]);

        if (diff_k_K > diff_k_B + 1.e-6) {
          simplices.erase(simplices.begin() + s_ind);
          s_ind--;
        }
      }
    }

    {
      for (std::size_t s_ind0 = 0; s_ind0 < simplices.size(); s_ind0++) {
        for (std::size_t s_ind1 = s_ind0 + 1; s_ind1 < simplices.size(); s_ind1++) {
          if (math::util::distance2(simplices[s_ind0].k_vec, simplices[s_ind1].k_vec) < 1.e-6) {
            simplices.erase(simplices.begin() + s_ind1);
            s_ind1--;
          }
        }
      }
    }
  }
}

template <class k_cluster_type>
void tetrahedron_mesh_initializer<2, k_cluster_type>::find_facets() {
  int coor[2];

  for (std::size_t l0 = 0; l0 < simplices.size(); l0++) {
    for (std::size_t l1 = 0; l1 < simplices.size(); l1++) {
      coor[0] = l0;
      coor[1] = l1;

      if (facet<DIMENSION>::is_facet(coor, simplices)) {
        facet<DIMENSION> f0;
        if (l0 < l1) {
          f0.index.push_back(l0);
          f0.index.push_back(l1);
        }
        else {
          f0.index.push_back(l1);
          f0.index.push_back(l0);
        }

        facets.push_back(f0);
      }
    }
  }

  for (std::size_t l0 = 0; l0 < facets.size(); l0++) {
    for (std::size_t l1 = l0 + 1; l1 < facets.size(); l1++) {
      if (facet<DIMENSION>::equal(facets[l0], facets[l1], simplices)) {
        facets.erase(facets.begin() + l1);
        l1--;
      }
    }
  }
}

template <class k_cluster_type>
void tetrahedron_mesh_initializer<2, k_cluster_type>::make_mesh_points() {
  // Initial mesh = [origin, simplex_1, simplex_2, ...]
  // Note that the indices of the simplices are shifted by +1 w.r.t. to the simplices vector.
  // Add origin.
  mesh.resize(1, std::vector<double>(DIMENSION, 0.));

  // Add simplices.
  for (std::size_t l = 0; l < simplices.size(); l++)
    mesh.push_back(simplices[l].k_vec);

  // Create tetrahedra formed by connecting the facets to the origin.
  for (std::size_t l = 0; l < facets.size(); l++) {
    tetrahedron<DIMENSION> tet;

    {
      // Compute the positions of the three vectors defining the tetrahedron in the mesh vector.
      // Note that tet.index corresponds to the mesh vector, while facets[l].index corresponds to
      // simplices vector.
      tet.index[0] = 0;                       // Position of the origin in the mesh vector.
      tet.index[1] = facets[l].index[0] + 1;  // Position of the first simplex defining the facet.
      tet.index[2] = facets[l].index[1] + 1;  // Position of second simplex defining the facet.

      // Set the three vectors.
      tet.vec_0 = mesh[tet.index[0]];
      tet.vec_1 = mesh[tet.index[1]];
      tet.vec_2 = mesh[tet.index[2]];

      // Compute the volume of tetrahedron.
      tet.volume = tet.compute_volume(&tet.vec_0[0], &tet.vec_1[0], &tet.vec_2[0]);
    }

    // Compute the normal.
    {
      std::vector<double> normal(DIMENSION, 0.);

      for (std::size_t i = 0; i < facets[l].index.size(); i++) {
        normal[0] += simplices[facets[l].index[i]].k_vec[0] / double(facets[l].index.size());
        normal[1] += simplices[facets[l].index[i]].k_vec[1] / double(facets[l].index.size());
      }

      tet.normal = normal;
    }

    tetrahedra.push_back(tet);
  }

  tetrahedra.reserve(int(tetrahedra.size()) * int(std::pow(8., N_recursion)));
  tetrahedra.reserve(int(4 * tetrahedra.size()) * int(std::pow(2., N_recursion)));

  for (int i = 0; i < N_recursion; i++) {
    int n_tet = tetrahedra.size();
    for (int l = 0; l < n_tet; l++)
      tetrahedra[l].do_recursion(tetrahedra, mesh);

    tetrahedra.erase(tetrahedra.begin(), tetrahedra.begin() + n_tet);
    //       for(int l=0; l<n_tet; l++)
    // 	tetrahedra.erase(tetrahedra.begin());
  }

  {  // get rid of mesh-redundancy
    std::vector<std::vector<double>>::iterator it;
    std::vector<std::vector<double>> mesh_old = mesh;

    std::sort(mesh.begin(), mesh.end());
    it = std::unique(mesh.begin(), mesh.end());
    mesh.erase(it, mesh.end());

    std::vector<int> index(mesh_old.size(), -1);

    for (std::size_t i = 0; i < mesh_old.size(); i++) {
      it = lower_bound(mesh.begin(), mesh.end(), mesh_old[i]);  // --> complexity log(N) !
      index[i] = it - mesh.begin();
      assert(index[i] < int(mesh.size()));
    }

    for (std::size_t l = 0; l < tetrahedra.size(); l++)
      for (int z = 0; z < DIMENSION + 1; z++)
        tetrahedra[l].index[z] = index[tetrahedra[l].index[z]];
  }
}

//
// Definition of member functions of the 3D tetrahedron mesh initializer.
//
template <class k_cluster_type>
tetrahedron_mesh_initializer<3, k_cluster_type>::tetrahedron_mesh_initializer(
    std::vector<simplex_t>& simplices_ref, std::vector<facet_t>& facets_ref,
    std::vector<vector_t>& mesh_ref, std::vector<tetrahedron_t>& tetrahedra_ref, int N_recursion_ref)
    : simplices(simplices_ref),
      facets(facets_ref),
      mesh(mesh_ref),
      tetrahedra(tetrahedra_ref),
      N_recursion(N_recursion_ref) {}

template <class k_cluster_type>
void tetrahedron_mesh_initializer<3, k_cluster_type>::execute() {
  make_convex_hull();
  find_facets();
  make_mesh_points();
}

template <class k_cluster_type>
void tetrahedron_mesh_initializer<3, k_cluster_type>::make_convex_hull() {
  std::vector<double>& b0 = k_cluster_type::get_basis_vectors()[0];
  std::vector<double>& b1 = k_cluster_type::get_basis_vectors()[1];
  std::vector<double>& b2 = k_cluster_type::get_basis_vectors()[2];

  std::vector<std::vector<double>> B_collection(0);
  {
    std::vector<double> B(3, 0);

    for (int t0 = -1; t0 <= 1; t0++) {
      for (int t1 = -1; t1 <= 1; t1++) {
        for (int t2 = -1; t2 <= 1; t2++) {
          if (t0 != 0 || t1 != 0 || t2 != 0) {
            B[0] = t0 * b0[0] + t1 * b1[0] + t2 * b2[0];
            B[1] = t0 * b0[1] + t1 * b1[1] + t2 * b2[1];
            B[2] = t0 * b0[2] + t1 * b1[2] + t2 * b2[2];

            B_collection.push_back(B);
          }
        }
      }
    }
  }

  double* A = new double[3 * 3];
  double* B = new double[3];

  for (std::size_t l0 = 0; l0 < B_collection.size(); l0++) {
    for (std::size_t l1 = 0; l1 < B_collection.size(); l1++) {
      for (std::size_t l2 = 0; l2 < B_collection.size(); l2++) {
        A[0 + 3 * 0] = B_collection[l0][0];
        A[0 + 3 * 1] = B_collection[l0][1];
        A[0 + 3 * 2] = B_collection[l0][2];
        A[1 + 3 * 0] = B_collection[l1][0];
        A[1 + 3 * 1] = B_collection[l1][1];
        A[1 + 3 * 2] = B_collection[l1][2];
        A[2 + 3 * 0] = B_collection[l2][0];
        A[2 + 3 * 1] = B_collection[l2][1];
        A[2 + 3 * 2] = B_collection[l2][2];

        B[0] = B_collection[l0][0] * B_collection[l0][0] / 2. +
               B_collection[l0][1] * B_collection[l0][1] / 2. +
               B_collection[l0][2] * B_collection[l0][2] / 2.;
        B[1] = B_collection[l1][0] * B_collection[l1][0] / 2. +
               B_collection[l1][1] * B_collection[l1][1] / 2. +
               B_collection[l1][2] * B_collection[l1][2] / 2.;
        B[2] = B_collection[l2][0] * B_collection[l2][0] / 2. +
               B_collection[l2][1] * B_collection[l2][1] / 2. +
               B_collection[l2][2] * B_collection[l2][2] / 2.;

        try {
          linalg::lapack::solve(3, A, 3, B);  // overwrites A with its LU factorization.
        }
        catch (linalg::lapack::util::LapackException e) {
          if (e.info() < 0)
            // Argument error: re-throw.
            throw;
          // A is singular.
          continue;
        }

        double det_A = A[0] * A[4] * A[8];

        if (std::abs(det_A) > 1.e-6) {
          simplex<DIMENSION> s;
          s.k_vec.resize(3, 0);

          s.k_vec[0] = B[0];
          s.k_vec[1] = B[1];
          s.k_vec[2] = B[2];

          simplices.push_back(s);
        }
      }
    }
  }

  delete[] A;
  delete[] B;

  {
    std::vector<double> K(3, 0);

    for (std::size_t B_ind = 0; B_ind < B_collection.size(); B_ind++) {
      for (std::size_t s_ind = 0; s_ind < simplices.size(); s_ind++) {
        double diff_k_K = math::util::distance2(simplices[s_ind].k_vec, K);
        double diff_k_B = math::util::distance2(simplices[s_ind].k_vec, B_collection[B_ind]);

        if (diff_k_K > diff_k_B + 1.e-6) {
          simplices.erase(simplices.begin() + s_ind);
          s_ind--;
        }
      }
    }
  }

  {
    for (std::size_t s_ind0 = 0; s_ind0 < simplices.size(); s_ind0++) {
      for (std::size_t s_ind1 = s_ind0 + 1; s_ind1 < simplices.size(); s_ind1++) {
        if (math::util::distance2(simplices[s_ind0].k_vec, simplices[s_ind1].k_vec) < 1.e-6) {
          simplices.erase(simplices.begin() + s_ind1);
          s_ind1--;
        }
      }
    }
  }
}

template <class k_cluster_type>
void tetrahedron_mesh_initializer<3, k_cluster_type>::find_facets() {
  int coor[3];

  for (std::size_t l0 = 0; l0 < simplices.size(); l0++) {
    for (std::size_t l1 = 0; l1 < simplices.size(); l1++) {
      for (std::size_t l2 = 0; l2 < simplices.size(); l2++) {
        coor[0] = l0;
        coor[1] = l1;
        coor[2] = l2;

        if (facet<DIMENSION>::is_facet(coor, simplices)) {
          facet<DIMENSION> f0;
          f0.index.push_back(l0);
          f0.index.push_back(l1);
          f0.index.push_back(l2);

          facets.push_back(f0);
        }
      }
    }
  }

  for (std::size_t l0 = 0; l0 < facets.size(); l0++) {
    for (std::size_t l1 = l0 + 1; l1 < facets.size(); l1++) {
      bool equal = facet<DIMENSION>::equal(facets[l0], facets[l1], simplices);

      if (equal) {
        for (std::size_t l = 0; l < facets[l1].index.size(); l++)
          facets[l0].index.push_back(facets[l1].index[l]);

        facets.erase(facets.begin() + l1);
        l1--;
      }
    }
  }

  for (std::size_t l0 = 0; l0 < facets.size(); l0++) {
    std::sort(facets[l0].index.begin(), facets[l0].index.end());
    std::vector<int>::iterator it = std::unique(facets[l0].index.begin(), facets[l0].index.end());
    facets[l0].index.erase(it, facets[l0].index.end());
  }
}

// TODO: Test
template <class k_cluster_type>
void tetrahedron_mesh_initializer<3, k_cluster_type>::make_mesh_points() {
  assert(DIMENSION == 3);

  std::vector<int>::iterator result_i, result_j;

  // Initial mesh = [origin, centroid of facet 1, centroid of facet 2, ..., simplex 1, simplex 2,
  // ...].
  // Add origin.
  mesh.resize(1, std::vector<double>(3, 0.));

  // Add centroids of facets.
  for (std::size_t l = 0; l < facets.size(); l++) {
    std::vector<double> k(3, 0.);

    for (std::size_t i = 0; i < facets[l].index.size(); i++) {
      k[0] += simplices[facets[l].index[i]].k_vec[0] / double(facets[l].index.size());
      k[1] += simplices[facets[l].index[i]].k_vec[1] / double(facets[l].index.size());
      k[2] += simplices[facets[l].index[i]].k_vec[2] / double(facets[l].index.size());
    }

    mesh.push_back(k);
  }

  // Add simplices.
  for (std::size_t l = 0; l < simplices.size(); l++)
    mesh.push_back(simplices[l].k_vec);

  // Add tetrahedra formed by an edge of a facet, its centroid and the origin.
  for (std::size_t l = 0; l < facets.size(); l++) {
    for (std::size_t i = 0; i < facets[l].index.size(); i++) {
      for (std::size_t j = i + 1; j < facets[l].index.size(); j++) {
        bool form_edge_line = false;

        for (std::size_t k = 0; k < facets.size(); k++)
          if (k != l) {
            result_i = find(facets[k].index.begin(), facets[k].index.end(), facets[l].index[i]);
            result_j = find(facets[k].index.begin(), facets[k].index.end(), facets[l].index[j]);

            bool i_is_element_of_facet_k = result_i == facets[k].index.end() ? false : true;
            bool j_is_element_of_facet_k = result_j == facets[k].index.end() ? false : true;

            if (i_is_element_of_facet_k && j_is_element_of_facet_k)
              form_edge_line = true;
          }

        if (form_edge_line) {
          tetrahedron<DIMENSION> tet;

          {
            tet.index[0] = 0;      // Position of the origin in mesh vector.
            tet.index[1] = l + 1;  // Position of centroid of current facet.
            // Position of first simplex defining the current edge of the facet.
            tet.index[2] = facets[l].index[i] + 1 + facets.size();
            // Position of second simplex defining the current edge of the facet.
            tet.index[3] = facets[l].index[j] + 1 + facets.size();

            tet.vec_0 = mesh[tet.index[0]];
            tet.vec_1 = mesh[tet.index[1]];
            tet.vec_2 = mesh[tet.index[2]];
            tet.vec_3 = mesh[tet.index[3]];

            tet.volume =
                tet.compute_volume(&tet.vec_0[0], &tet.vec_1[0], &tet.vec_2[0], &tet.vec_3[0]);
          }

          // Compute the normal.
          {
            std::vector<double> normal(3, 0.);

            for (std::size_t i = 0; i < facets[l].index.size(); i++) {
              normal[0] += simplices[facets[l].index[i]].k_vec[0] / double(facets[l].index.size());
              normal[1] += simplices[facets[l].index[i]].k_vec[1] / double(facets[l].index.size());
              normal[2] += simplices[facets[l].index[i]].k_vec[2] / double(facets[l].index.size());
            }

            tet.normal = normal;
          }
          tetrahedra.push_back(tet);
        }
      }
    }
  }

  std::cout << "tetrahedra.size()" << tetrahedra.size() << std::endl;

  tetrahedra.reserve(int(tetrahedra.size()) * int(std::pow(8., N_recursion)));
  tetrahedra.reserve(int(4 * tetrahedra.size()) * int(std::pow(2., N_recursion)));

  for (int i = 0; i < N_recursion; i++) {
    int n_tet = tetrahedra.size();
    for (int l = 0; l < n_tet; l++)
      tetrahedra[l].do_recursion(tetrahedra, mesh);

    tetrahedra.erase(tetrahedra.begin(), tetrahedra.begin() + n_tet);
  }

  {  // get rid of mesh-redundancy
    std::vector<std::vector<double>>::iterator it;
    std::vector<std::vector<double>> mesh_old = mesh;

    std::sort(mesh.begin(), mesh.end());
    it = std::unique(mesh.begin(), mesh.end());
    mesh.erase(it, mesh.end());

    std::vector<int> index(mesh_old.size(), -1);

    for (std::size_t i = 0; i < mesh_old.size(); i++) {
      it = lower_bound(mesh.begin(), mesh.end(), mesh_old[i]);  // --> complexity log(N) !
      index[i] = it - mesh.begin();
      assert(index[i] < int(mesh.size()));
    }

    for (std::size_t l = 0; l < tetrahedra.size(); l++)
      for (int z = 0; z < 3 + 1; z++)
        tetrahedra[l].index[z] = index[tetrahedra[l].index[z]];
  }
}

}  // geometry
}  // math
}  // dca

#endif  // DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_HPP
