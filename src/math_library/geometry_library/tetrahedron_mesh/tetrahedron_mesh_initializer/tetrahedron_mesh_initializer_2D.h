// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization for the initialization of a 2D tetrahedron mesh.

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_TETRAHEDRON_MESH_INITIALIZER_2D_H
#define MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_TETRAHEDRON_MESH_INITIALIZER_2D_H

#include "math_library/geometry_library/tetrahedron_mesh/tetrahedron_mesh_initializer/tetrahedron_mesh_initializer_template.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "math_library/geometry_library/tetrahedron_mesh/facet.h"
#include "math_library/geometry_library/tetrahedron_mesh/tetrahedron/tetrahedron_2D.h"
#include "math_library/geometry_library/tetrahedron_mesh/simplex.h"
#include "math_library/geometry_library/vector_operations/vector_operations.hpp"

namespace math_algorithms {

template <class k_cluster_type>
class tetrahedron_mesh_initializer<2, k_cluster_type> {
  const static int DIMENSION = 2;

  typedef tetrahedron<2> tetrahedron_t;
  typedef simplex<2> simplex_t;
  typedef facet<2> facet_t;
  typedef std::vector<double> vector_t;

public:
  tetrahedron_mesh_initializer(std::vector<simplex_t>& simplices_ref,
                               std::vector<facet_t>& facets_ref, std::vector<vector_t>& mesh_ref,
                               std::vector<tetrahedron_t>& tetrahedra_ref, int N_recursion_ref);

  void execute();

private:
  void make_convex_hull();
  void find_facets();
  void make_mesh_points();

private:
  std::vector<simplex_t>& simplices;
  std::vector<facet_t>& facets;

  std::vector<vector_t>& mesh;
  std::vector<tetrahedron_t>& tetrahedra;

  int N_recursion;
};

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

  solve_plan<double> slv_pln(2, 1);

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

      {
        memcpy(slv_pln.matrix, A, sizeof(double) * 2 * 2);
        memcpy(slv_pln.solved_matrix, B, sizeof(double) * 2);

        slv_pln.execute_plan();

        double det_A = slv_pln.matrix[0] * slv_pln.matrix[3];

        if (std::fabs(det_A) > 1.e-6) {
          simplex<DIMENSION> s;
          s.k_vec.resize(2, 0);

          s.k_vec[0] = slv_pln.solved_matrix[0];
          s.k_vec[1] = slv_pln.solved_matrix[1];

          simplices.push_back(s);
        }
      }
    }
  }

  delete[] A;
  delete[] B;

  {
    std::vector<double> K(2, 0);

    for (std::size_t B_ind = 0; B_ind < B_collection.size(); B_ind++) {
      for (std::size_t s_ind = 0; s_ind < simplices.size(); s_ind++) {
        double diff_k_K = VECTOR_OPERATIONS::L2_NORM(simplices[s_ind].k_vec, K);
        double diff_k_B = VECTOR_OPERATIONS::L2_NORM(simplices[s_ind].k_vec, B_collection[B_ind]);

        if (diff_k_K > diff_k_B + 1.e-6) {
          simplices.erase(simplices.begin() + s_ind);
          s_ind--;
        }
      }
    }

    {
      for (std::size_t s_ind0 = 0; s_ind0 < simplices.size(); s_ind0++) {
        for (std::size_t s_ind1 = s_ind0 + 1; s_ind1 < simplices.size(); s_ind1++) {
          if (VECTOR_OPERATIONS::L2_NORM(simplices[s_ind0].k_vec, simplices[s_ind1].k_vec) < 1.e-6) {
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
  mesh.resize(1, std::vector<double>(DIMENSION, 0.));

  for (std::size_t l = 0; l < simplices.size(); l++)
    mesh.push_back(simplices[l].k_vec);

  for (std::size_t l = 0; l < facets.size(); l++) {
    tetrahedron<DIMENSION> tet;

    {
      tet.index[0] = 0;
      tet.index[1] = facets[l].index[0] + 1;
      tet.index[2] = facets[l].index[1] + 1;
    }

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
}
#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_TETRAHEDRON_MESH_INITIALIZER_2D_H
