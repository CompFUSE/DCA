// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Tetrahedron integration class.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_INTEGRATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_INTEGRATION_HPP

#include <complex>
#include <functional>
#include <future>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/parallel/stdthread/thread_pool/thread_pool.hpp"
#include "dca/parallel/util/get_bounds.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_domain.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_integration_data.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_routines_inverse_matrix_function.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename parameters_type>
class tetrahedron_integration {
public:
  using ThisType = tetrahedron_integration<parameters_type>;
  using K_dmn = typename ClusterDomainAliases<parameters_type::lattice_dimension>::KClusterDmn;

  using k_cluster_type = typename K_dmn::parameter_type;

  using tet_dmn_type = func::dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_K>>;
  using tet_0_dmn_type = func::dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using Threading = typename parameters_type::ThreadingType;

  const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

public:
  tetrahedron_integration(parameters_type& parameters_ref);

  template <typename scalar_type>
  void execute(func::function<scalar_type, tet_dmn_type>& w_tet,
               func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet,
               func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>>& G_int) const;

protected:
  template <typename scalar_type>
  auto tetrahedron_integration_2D(
      const int id, const int nr_threads, func::function<scalar_type, tet_dmn_type>& w_tet,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet) const;

  template <typename scalar_type>
  auto tetrahedron_integration_3D(
      const int id, const int nr_threads, func::function<scalar_type, tet_dmn_type>& w_tet,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet) const;

  // TODO: implement 1D version.
  //  template <typename scalar_type>
  //  void tetrahedron_integration_1D(
  //      func::function<scalar_type, tet_dmn_type>& w_tet,
  //      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>&
  //      G_tet,
  //      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>>& G_int) const;

private:
  parameters_type& parameters;
};

template <typename parameters_type>
tetrahedron_integration<parameters_type>::tetrahedron_integration(parameters_type& parameters_ref)
    : parameters(parameters_ref) {}

template <typename parameters_type>
template <typename scalar_type>
void tetrahedron_integration<parameters_type>::execute(
    func::function<scalar_type, tet_dmn_type>& w_tet,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>>& G_int) const {
  const int nr_threads = parameters.get_coarsegraining_threads();

  std::function<void(int)> task;

  switch (DIMENSION) {
    //       case 1:
    //         parallelization_obj.execute(nr_threads, tetrahedron_integration_mt_1D<scalar_type>,
    //         (void*) &tetrahedron_integration_functions_obj);
    //         break;

    case 2:
      task =
          std::bind(&ThisType::tetrahedron_integration_2D<scalar_type>, this, std::placeholders::_1,
                    std::placeholders::_2, std::ref(w_tet), std::ref(G_tet));
      break;

    case 3:
      task =
          std::bind(&ThisType::tetrahedron_integration_3D<scalar_type>, this, std::placeholders::_1,
                    std::placeholders::_2, std::ref(w_tet), std::ref(G_tet));
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  Threading threads;
  G_int = threads.sumReduction(nr_threads, task);
}

// template <typename parameters_type>
// template <typename scalar_type>
// void tetrahedron_integration<parameters_type>::tetrahedron_integration_1D(
//    func::function<scalar_type, tet_dmn_type>& w_tet,
//    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet,
//    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>>& G_int) {
//  for (int j = 0; j < nu::dmn_size(); j++)
//    for (int i = 0; i < nu::dmn_size(); i++)
//      G_int(i, j) = 0;
//
//  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>> G_tmp("G_tmp");
//
//  tetrahedron_integration_data<scalar_type> data_obj(nu::dmn_size());
//
//  for (int tet_ind = 0; tet_ind < tet_dmn_type::dmn_size(); tet_ind += 4) {
//    scalar_type volume =
//        w_tet(tet_ind) + w_tet(tet_ind + 1) + w_tet(tet_ind + 2) + w_tet(tet_ind + 3);
//
//    std::complex<scalar_type>* G_0 = &G_tet(0, 0, tet_ind + 0);
//    std::complex<scalar_type>* G_1 = &G_tet(0, 0, tet_ind + 1);
//
//    tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1,
//                                                          &G_tmp(0, 0), data_obj);
//
//    for (int j = 0; j < nu::dmn_size(); j++)
//      for (int i = 0; i < nu::dmn_size(); i++)
//        G_int(i, j) += G_tmp(i, j);
//  }
//}

template <typename parameters_type>
template <typename scalar_type>
auto tetrahedron_integration<parameters_type>::tetrahedron_integration_2D(
    const int id, const int nr_threads, func::function<scalar_type, tet_dmn_type>& w_tet,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet) const {
  tet_dmn_type tet_dmn;
  std::pair<int, int> tet_bounds = dca::parallel::util::getBounds(id, nr_threads, tet_dmn);

  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>> G_int;

  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>> G_tmp("G_tmp");

  tetrahedron_integration_data<scalar_type> data_obj(nu::dmn_size());

  for (int tet_ind = tet_bounds.first; tet_ind < tet_bounds.second; tet_ind += 3) {
    scalar_type volume = w_tet(tet_ind) + w_tet(tet_ind + 1) + w_tet(tet_ind + 2);

    std::complex<scalar_type>* G_0 = &G_tet(0, 0, tet_ind + 0);
    std::complex<scalar_type>* G_1 = &G_tet(0, 0, tet_ind + 1);
    std::complex<scalar_type>* G_2 = &G_tet(0, 0, tet_ind + 2);

    tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, G_2,
                                                          &G_tmp(0, 0), data_obj);

    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        G_int(i, j) += G_tmp(i, j);
  }

  return G_int;
}

template <typename parameters_type>
template <typename scalar_type>
auto tetrahedron_integration<parameters_type>::tetrahedron_integration_3D(
    const int id, const int nr_threads, func::function<scalar_type, tet_dmn_type>& w_tet,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet) const {
  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>> G_int;
  tet_dmn_type tet_dmn;
  std::pair<int, int> tet_bounds = dca::parallel::util::getBounds(id, nr_threads, tet_dmn);

  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>> G_tmp("G_tmp");

  tetrahedron_integration_data<scalar_type> data_obj(nu::dmn_size());

  for (int tet_ind = tet_bounds.first; tet_ind < tet_bounds.second; tet_ind += 4) {
    scalar_type volume =
        w_tet(tet_ind) + w_tet(tet_ind + 1) + w_tet(tet_ind + 2) + w_tet(tet_ind + 3);

    std::complex<scalar_type>* G_0 = &G_tet(0, 0, tet_ind + 0);
    std::complex<scalar_type>* G_1 = &G_tet(0, 0, tet_ind + 1);
    std::complex<scalar_type>* G_2 = &G_tet(0, 0, tet_ind + 2);
    std::complex<scalar_type>* G_3 = &G_tet(0, 0, tet_ind + 3);

    tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, G_2,
                                                          G_3, &G_tmp(0, 0), data_obj);

    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        G_int(i, j) += G_tmp(i, j);
  }
  return G_int;
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_INTEGRATION_HPP
