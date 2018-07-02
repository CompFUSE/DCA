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
#include "dca/parallel/util/get_bounds.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_domain.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_integration_data.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_routines_inverse_matrix_function.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename parameters_type, typename K_dmn>
class tetrahedron_integration {
public:
  using k_cluster_type = typename K_dmn::parameter_type;

  using tet_dmn_type = func::dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_K>>;
  using tet_0_dmn_type = func::dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

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

template <typename parameters_type, typename K_dmn>
tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration(parameters_type& parameters_ref)
    : parameters(parameters_ref) {}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type>
void tetrahedron_integration<parameters_type, K_dmn>::execute(
    func::function<scalar_type, tet_dmn_type>& w_tet,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, tet_dmn_type>>& G_tet,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>>& G_int) const {
  const int nr_threads = parameters.get_coarsegraining_threads();

  const auto policy = nr_threads > 1 ? std::launch::async : std::launch::deferred;
  std::function<void(int)> task;

  switch (DIMENSION) {
    //       case 1:
    //         parallelization_obj.execute(nr_threads, tetrahedron_integration_mt_1D<scalar_type>,
    //         (void*) &tetrahedron_integration_functions_obj);
    //         break;

    case 2:
      task = std::bind(tetrahedron_integration_2D<scalar_type>, std::placeholders::_1, nr_threads,
                       std::ref(w_tet), std::ref(G_tet));
      break;

    case 3:
      task = std::bind(tetrahedron_integration_3D<scalar_type>, std::placeholders::_1, nr_threads,
                       std::ref(w_tet), std::ref(G_tet));
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  std::vector<std::future<func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu>>>> futures;
  for (int id = 0; id < nr_threads; id++) {
    futures.emplace_back(std::async(policy, task, id));
  }

  G_int = 0;
  for (int l = 0; l < nr_threads; l++) {
    G_int += futures[l].get();
  }
}

// template <typename parameters_type, typename K_dmn>
// template <typename scalar_type>
// void tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_1D(
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

template <typename parameters_type, typename K_dmn>
template <typename scalar_type>
auto tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_2D(
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

template <typename parameters_type, typename K_dmn>
template <typename scalar_type>
auto tetrahedron_integration<parameters_type, K_dmn>::tetrahedron_integration_3D(
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

    /*
    {
      std::vector<double> x(3);
      std::vector<double> y(3);
      std::vector<double> z(3);

      for(int d=0; d<3; d++){
        x[d] =
    tet_dmn_type::get_elements()[tet_ind+1][d]-tet_dmn_type::get_elements()[tet_ind+0][d];
        y[d] =
    tet_dmn_type::get_elements()[tet_ind+2][d]-tet_dmn_type::get_elements()[tet_ind+0][d];
        z[d] =
    tet_dmn_type::get_elements()[tet_ind+3][d]-tet_dmn_type::get_elements()[tet_ind+0][d];
      }

      if(abs(volume-math::util::volume(x,y,z)/6)>1.e-6)
        {
          std::cout << tet_ind << "\t" << volume << "\t" << math::util::volume(x,y,z) <<
    "\n";

          math::util::print(x); std::cout << "\n";
          math::util::print(y); std::cout << "\n";
          math::util::print(z); std::cout << "\n";
          std::cout << "\n\n";

          math::util::print(tet_dmn_type::get_elements()[tet_ind+0]); std::cout << "\t" <<
    G_tet(0,0,tet_ind+0) << "\n";
          math::util::print(tet_dmn_type::get_elements()[tet_ind+1]); std::cout << "\t" <<
    G_tet(0,0,tet_ind+1) << "\n";
          math::util::print(tet_dmn_type::get_elements()[tet_ind+2]); std::cout << "\t" <<
    G_tet(0,0,tet_ind+2) << "\n";
          math::util::print(tet_dmn_type::get_elements()[tet_ind+3]); std::cout << "\t" <<
    G_tet(0,0,tet_ind+3) << "\n";
          std::cout << "\n\n";

          assert(false);
        }
    }
    */

    std::complex<scalar_type>* G_0 = &G_tet(0, 0, tet_ind + 0);
    std::complex<scalar_type>* G_1 = &G_tet(0, 0, tet_ind + 1);
    std::complex<scalar_type>* G_2 = &G_tet(0, 0, tet_ind + 2);
    std::complex<scalar_type>* G_3 = &G_tet(0, 0, tet_ind + 3);

    tetrahedron_routines_inverse_matrix_function::execute(nu::dmn_size(), volume, G_0, G_1, G_2,
                                                          G_3, &G_tmp(0, 0), data_obj);

    //         for(int j=0; j<nu::dmn_size(); j++)
    //           for(int i=0; i<nu::dmn_size(); i++)
    // 	    G_tmp(i,j) = (G_tet(i,j,tet_ind+0)
    // 			  +G_tet(i,j,tet_ind+1)
    // 			  +G_tet(i,j,tet_ind+2)
    // 			  +G_tet(i,j,tet_ind+3))*volume/4.;

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
