// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Gaussian quadrature domain.

#ifndef DCA_MATH_GEOMETRY_GAUSSIAN_QUADRATURE_GAUSSIAN_QUADRATURE_DOMAIN_HPP
#define DCA_MATH_GEOMETRY_GAUSSIAN_QUADRATURE_GAUSSIAN_QUADRATURE_DOMAIN_HPP

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/math/function_transform/domain_specifications.hpp"
#include "dca/math/geometry/gaussian_quadrature/compute_weights_and_abscissas.hpp"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron_mesh.hpp"
#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

template <typename dmn_type>
class gaussian_quadrature_domain {};

template <typename cluster_type>
class gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>> {
public:
  const static int DIMENSION = cluster_type::DIMENSION;

  using dmn_specifications_type =
      transform::domain_specifications<double, std::vector<double>, transform::DISCRETE,
                                       transform::KRONECKER_DELTA, transform::INTERVAL,
                                       transform::NONEQUIDISTANT>;

  typedef typename dmn_specifications_type::scalar_type scalar_type;
  typedef typename dmn_specifications_type::element_type element_type;

public:
  static bool& is_initialized();

  static int& get_size();

  static std::string& get_name();

  static std::vector<scalar_type>& get_weights();
  static std::vector<element_type>& get_elements();

  static void reset();

  static void initialize_Brillouin_zone(int N_recursion, int rule);
  static void initialize_parallellepipedum(int N_recursion, int rule);

  static void initialize_Brillouin_zone(int N_recursion, int rule, double period);

  static void initialize_flat_mesh(tetrahedron_mesh<cluster_type>& mesh);

  static void initialize_gaussian_mesh(int rule, tetrahedron_mesh<cluster_type>& mesh);

  static void translate_according_to_period(double period, tetrahedron_mesh<cluster_type>& mesh);

  static void initialize_elements(tetrahedron_mesh<cluster_type>& mesh);

  static void plot_tetrahedra(tetrahedron_mesh<cluster_type>& mesh);

  static void plot_q_vecs(tetrahedron_mesh<cluster_type>& mesh);
};

template <typename cluster_type>
bool& gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::is_initialized() {
  static bool initialized = false;
  return initialized;
}

template <typename cluster_type>
int& gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::get_size() {
  static int size = 0;
  return size;
}

template <typename cluster_type>
std::string& gaussian_quadrature_domain<
    func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::get_name() {
  static std::string name = "gaussian_quadrature_domain";
  return name;
}

template <typename cluster_type>
std::vector<typename gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::scalar_type>& gaussian_quadrature_domain<
    func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::get_weights() {
  static std::vector<scalar_type> weights(0);
  return weights;
}

template <typename cluster_type>
std::vector<typename gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::element_type>& gaussian_quadrature_domain<
    func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::get_elements() {
  static std::vector<element_type> elements(0);
  return elements;
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::reset() {
  get_size() = 0;
  get_name() = "";

  get_elements().resize(0);

  is_initialized() = false;
}

template <typename cluster_type>
void gaussian_quadrature_domain<
    func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::initialize_Brillouin_zone(int N_recursion,
                                                                                         int rule) {
  tetrahedron_mesh<cluster_type> mesh(N_recursion);

  if (rule < 0)
    initialize_flat_mesh(mesh);
  else
    initialize_gaussian_mesh(rule, mesh);

  initialize_elements(mesh);
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::initialize_Brillouin_zone(
    int N_recursion, int rule, double period) {
  tetrahedron_mesh<cluster_type> mesh(N_recursion);

  if (rule < 0)
    initialize_flat_mesh(mesh);
  else
    initialize_gaussian_mesh(rule, mesh);

  translate_according_to_period(period, mesh);

  initialize_elements(mesh);

  if (false) {
    plot_tetrahedra(mesh);

    plot_q_vecs(mesh);
  }
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::initialize_flat_mesh(
    tetrahedron_mesh<cluster_type>& mesh) {
  std::cout << __FUNCTION__ << std::endl;

  std::vector<tetrahedron<DIMENSION>>& tets = mesh.get_tetrahedra();

  for (int l = 0; l < tets.size(); l++) {
    std::vector<double> cm = tets[l].compute_cm();

    tets[l].N_q = 1;

    tets[l].q_w = new double[1];
    tets[l].q_vecs = new double[DIMENSION];

    tets[l].q_w[0] = 1.;
    for (int d = 0; d < DIMENSION; d++)
      tets[l].q_vecs[d] = cm[d];
  }

  get_size() = get_elements().size();
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::initialize_gaussian_mesh(
    int rule, tetrahedron_mesh<cluster_type>& mesh) {
  std::vector<tetrahedron<DIMENSION>>& tets = mesh.get_tetrahedra();

  for (int l = 0; l < tets.size(); l++)
    computeWeightsAndAbscissas(rule, tets[l]);
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::
    translate_according_to_period(double period, tetrahedron_mesh<cluster_type>& mesh) {
  std::vector<double> cm(DIMENSION, 0.);
  std::vector<double> normal(DIMENSION, 0.);

  std::vector<tetrahedron<DIMENSION>>& tetrahedra = mesh.get_tetrahedra();

  for (size_t j = 0; j < tetrahedra.size(); j++) {
    cm = tetrahedra[j].compute_cm();
    normal = tetrahedra[j].normal;

    double x = math::util::innerProduct(cm, normal) / math::util::l2Norm2(normal);
    assert(x > -1.e-6 and x < 1. + 1.e-6);

    for (int d = 0; d < DIMENSION; d++)
      normal[d] *= -2.;

    if (std::sin(period * (2. * M_PI) * x) < 0.)
      tetrahedra[j].translate(normal);
  }
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::initialize_elements(
    tetrahedron_mesh<cluster_type>& mesh) {
  get_weights().resize(0);
  get_elements().resize(0);

  std::vector<tetrahedron<DIMENSION>>& tetrahedra = mesh.get_tetrahedra();

  std::vector<std::vector<double>> k_points(0);
  for (size_t l = 0; l < tetrahedra.size(); l++) {
    std::vector<double> k_point(DIMENSION + 1, 0);

    for (int i = 0; i < tetrahedra[l].N_q; i++) {
      for (int d = 0; d < DIMENSION; d++)
        k_point[d] = tetrahedra[l].q_vecs[d + i * DIMENSION];

      k_point[DIMENSION] = (tetrahedra[l].q_w[i]) * (tetrahedra[l].volume);

      k_points.push_back(k_point);
    }
  }

  std::sort(k_points.begin(), k_points.end());

  for (size_t l = 0; l < k_points.size(); l++) {
    std::vector<double> k_point(DIMENSION, 0);

    for (int d = 0; d < DIMENSION; d++)
      k_point[d] = k_points[l][d];

    get_elements().push_back(k_point);
    get_weights().push_back(k_points[l][DIMENSION]);
  }

  get_size() = get_elements().size();
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::plot_tetrahedra(
    tetrahedron_mesh<cluster_type>& mesh) {
  std::vector<tetrahedron<DIMENSION>>& tetrahedra = mesh.get_tetrahedra();

  dca::util::Plot plot("lines");

  for (size_t j = 0; j < tetrahedra.size(); j++)
    tetrahedra[j].plot(plot);
}

template <typename cluster_type>
void gaussian_quadrature_domain<func::dmn_0<tetrahedron_mesh<func::dmn_0<cluster_type>>>>::plot_q_vecs(
    tetrahedron_mesh<cluster_type>& mesh) {
  std::vector<tetrahedron<DIMENSION>>& tetrahedra = mesh.get_tetrahedra();

  dca::util::Plot plot("points");

  for (size_t j = 0; j < tetrahedra.size(); j++)
    tetrahedra[j].plot_q_vecs(plot);
}

}  // geometry
}  // math
}  // dca

#endif  // DCA_MATH_GEOMETRY_GAUSSIAN_QUADRATURE_GAUSSIAN_QUADRATURE_DOMAIN_HPP
