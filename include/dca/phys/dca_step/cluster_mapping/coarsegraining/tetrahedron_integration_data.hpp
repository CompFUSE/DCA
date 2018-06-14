// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Helper class for the tetrahedron integration.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_INTEGRATION_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_INTEGRATION_DATA_HPP

#include <algorithm>
#include <complex>

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename scalar_type>
class tetrahedron_integration_data {
public:
  tetrahedron_integration_data(int N);
  ~tetrahedron_integration_data();

public:
  int INFO;

  std::complex<scalar_type>* G_inv_0;
  std::complex<scalar_type>* G_inv_1;
  std::complex<scalar_type>* G_inv_2;
  std::complex<scalar_type>* G_inv_3;

  std::complex<scalar_type>* VR_0;
  std::complex<scalar_type>* VR_1;
  std::complex<scalar_type>* VR_2;
  std::complex<scalar_type>* VR_3;

  std::complex<scalar_type>* VR_inv_0;
  std::complex<scalar_type>* VR_inv_1;
  std::complex<scalar_type>* VR_inv_2;
  std::complex<scalar_type>* VR_inv_3;

  std::complex<scalar_type>* W_0;
  std::complex<scalar_type>* W_1;
  std::complex<scalar_type>* W_2;
  std::complex<scalar_type>* W_3;

  // inverse calculation workspaces
  int inv_lwork;
  int* inv_ipiv;
  std::complex<scalar_type>* inv_work;

  // eigensolver calculation
  int eig_lwork;
  std::complex<scalar_type>* eig_work;
  scalar_type* eig_rwork;
};

template <typename scalar_type>
tetrahedron_integration_data<scalar_type>::tetrahedron_integration_data(int N)
    : G_inv_0(NULL),
      G_inv_1(NULL),
      G_inv_2(NULL),
      G_inv_3(NULL),

      VR_0(NULL),
      VR_1(NULL),
      VR_2(NULL),
      VR_3(NULL),

      VR_inv_0(NULL),
      VR_inv_1(NULL),
      VR_inv_2(NULL),
      VR_inv_3(NULL),

      W_0(NULL),
      W_1(NULL),
      W_2(NULL),
      W_3(NULL),

      // inverse calculation workspaces
      inv_lwork(-1),
      inv_ipiv(nullptr),
      inv_work(nullptr),

      // eigensolver calculation
      eig_lwork(-1),
      eig_work(nullptr),
      eig_rwork(nullptr) {
  G_inv_0 = new std::complex<scalar_type>[N * N];
  G_inv_1 = new std::complex<scalar_type>[N * N];
  G_inv_2 = new std::complex<scalar_type>[N * N];
  G_inv_3 = new std::complex<scalar_type>[N * N];

  VR_0 = new std::complex<scalar_type>[N * N];
  VR_1 = new std::complex<scalar_type>[N * N];
  VR_2 = new std::complex<scalar_type>[N * N];
  VR_3 = new std::complex<scalar_type>[N * N];

  VR_inv_0 = new std::complex<scalar_type>[N * N];
  VR_inv_1 = new std::complex<scalar_type>[N * N];
  VR_inv_2 = new std::complex<scalar_type>[N * N];
  VR_inv_3 = new std::complex<scalar_type>[N * N];

  W_0 = new std::complex<scalar_type>[N];
  W_1 = new std::complex<scalar_type>[N];
  W_2 = new std::complex<scalar_type>[N];
  W_3 = new std::complex<scalar_type>[N];

  inv_lwork = 128 * std::max(1, N);
  inv_ipiv = new int[N];
  inv_work = new std::complex<scalar_type>[inv_lwork];

  eig_lwork = 128 * std::max(1, 2 * N);
  eig_rwork = new scalar_type[std::max(1, 2 * N)];
  eig_work = new std::complex<scalar_type>[eig_lwork];
}

template <typename scalar_type>
tetrahedron_integration_data<scalar_type>::~tetrahedron_integration_data() {
  delete[] G_inv_0;
  delete[] G_inv_1;
  delete[] G_inv_2;
  delete[] G_inv_3;

  delete[] VR_0;
  delete[] VR_1;
  delete[] VR_2;
  delete[] VR_3;

  delete[] VR_inv_0;
  delete[] VR_inv_1;
  delete[] VR_inv_2;
  delete[] VR_inv_3;

  delete[] W_0;
  delete[] W_1;
  delete[] W_2;
  delete[] W_3;

  delete[] inv_ipiv;
  delete[] inv_work;

  delete[] eig_rwork;
  delete[] eig_work;
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_INTEGRATION_DATA_HPP
