// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the Wannier interpolation kernel using the NFFT library.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_KERNEL_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_KERNEL_HPP

#include <cassert>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <nfft3.h>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Empty template class declaration.
template <typename source_dmn_type, typename target_dmn_type>
class wannier_interpolation_kernel {};

// Template specialization for Wannier interpolation in momentum space.
template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
class wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type> {
  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_cluster_type;

  const static int DIMENSION = k_cluster_type::DIMENSION;

  typedef typename k_cluster_type::dual_type source_r_cluster_type;
  typedef typename k_cluster_type::this_type source_k_cluster_type;

  typedef func::dmn_0<source_r_cluster_type> source_r_dmn_t;
  typedef func::dmn_0<source_k_cluster_type> source_k_dmn_t;

  typedef func::dmn_0<target_dmn_type> target_k_dmn_t;

public:
  wannier_interpolation_kernel();
  ~wannier_interpolation_kernel();

  void reset_input();
  void reset_output();

  void set(std::complex<double>* input_ptr);
  void get(std::complex<double>* output_ptr);

  void execute(std::complex<double>* input_ptr, std::complex<double>* output_ptr);

  std::complex<double>& get_F_r(int i);

private:
  // Checks whether the grid size is larger or equal to the size of the source-k-domain.
  void check_grid_sizes();

  void initialize_centered_r_cluster();

  void initialize_nfft_K_2_R();
  void initialize_nfft_R_2_k();

  void initialize_cut_off();

  void FT_to_centered_function_NFFT();

  void FT_F_K__to__F_R(std::complex<double>* input_ptr);

  void FT_F_R__to__F_k(std::complex<double>* output_ptr);

private:
  struct centered_r_cluster {
    typedef centered_r_cluster this_type;
    typedef std::vector<double> element_type;

    static int get_size() {
      return get_elements().size();
    }

    static std::vector<element_type>& get_elements() {
      static std::vector<element_type> elements(0);
      return elements;
    }
  };

public:
  typedef func::dmn_0<centered_r_cluster> centered_r_cluster_dmn_t;

private:
  static bool INITIALIZED;

  static func::function<int, centered_r_cluster_dmn_t> lies_within_cutoff;

  std::vector<int> grid_size;

  nfft_plan nfft_K_2_R;
  nfft_plan nfft_R_2_k;

  func::function<std::complex<double>, centered_r_cluster_dmn_t> F_R;
};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
bool wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::INITIALIZED = false;

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
func::function<int, typename wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::centered_r_cluster_dmn_t> wannier_interpolation_kernel<
    cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
    target_dmn_type>::lies_within_cutoff("cutoff");

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                             target_dmn_type>::wannier_interpolation_kernel()
    : grid_size(DIMENSION, 32),

      nfft_K_2_R(),
      nfft_R_2_k(),

      F_R("wannier_interpolation_kernel__F_r") {
  for (int i = 0; i < DIMENSION; ++i)
    grid_size[i] = grid_size[i] <= 4 ? 6 : grid_size[i];

  check_grid_sizes();

  if (!INITIALIZED) {
    initialize_centered_r_cluster();

    // Reset functions since their domain has changed.
    F_R.reset();
    lies_within_cutoff.reset();

    initialize_cut_off();

    INITIALIZED = true;
  }

  initialize_nfft_K_2_R();
  initialize_nfft_R_2_k();
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                             target_dmn_type>::~wannier_interpolation_kernel() {
  nfft_finalize(&nfft_K_2_R);
  nfft_finalize(&nfft_R_2_k);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::check_grid_sizes() {
  int size = 1;
  for (int i = 0; i < DIMENSION; ++i)
    size *= grid_size[i];

  if (size < source_k_dmn_t::dmn_size()) {
    std::cout << "\n\n\t INCREASE LDA-grid FOR WANNIER-INTERPOLATION!!! \n\n\n";
    throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
std::complex<double>& wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                                   target_dmn_type>::get_F_r(int i) {
  return F_R(i);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::reset_input() {
  nfft_finalize(&nfft_K_2_R);
  initialize_nfft_K_2_R();
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::reset_output() {
  nfft_finalize(&nfft_R_2_k);
  initialize_nfft_R_2_k();
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::execute(std::complex<double>* input_ptr,
                                                            std::complex<double>* output_ptr) {
  FT_F_K__to__F_R(input_ptr);

  FT_F_R__to__F_k(output_ptr);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::set(std::complex<double>* input_ptr) {
  FT_F_K__to__F_R(input_ptr);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::get(std::complex<double>* output_ptr) {
  FT_F_R__to__F_k(output_ptr);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::initialize_nfft_K_2_R() {
  nfft_init(&nfft_K_2_R, DIMENSION, &grid_size[0], source_k_dmn_t::dmn_size());

  std::vector<std::vector<double>>& collection_k_vecs = source_k_dmn_t::get_elements();

  std::vector<std::vector<double>> collection_k_vecs_affine(source_k_dmn_t::dmn_size(),
                                                            std::vector<double>(DIMENSION, 0));

  for (int j = 0; j < source_k_dmn_t::dmn_size(); j++) {
    // collection_k_vecs_affine[j] =
    // source_k_cluster_type::get_affine_coordinate(collection_k_vecs[j]);
    collection_k_vecs_affine[j] = math::util::coordinates(
        collection_k_vecs[j], source_k_cluster_type::get_super_basis_vectors());

    // math::util::print(collection_k_vecs_affine[j]); cout<<endl;
  }
  // cout<<endl;

  for (int j = 0; j < source_k_dmn_t::dmn_size(); j++) {
    for (int i = 0; i < DIMENSION; i++) {
      while (collection_k_vecs_affine[j][i] < -1. / 2.)
        collection_k_vecs_affine[j][i] += 1.;

      while (collection_k_vecs_affine[j][i] > 1. / 2. - 1.e-6)
        collection_k_vecs_affine[j][i] -= 1.;
    }
  }

  for (int j = 0; j < source_k_dmn_t::dmn_size(); j++) {
    for (int i = 0; i < DIMENSION; i++) {
      nfft_K_2_R.x[j * DIMENSION + i] = collection_k_vecs_affine[j][i];

      assert(nfft_K_2_R.x[j * DIMENSION + i] >= -0.5 - 1.e-6);
      assert(nfft_K_2_R.x[j * DIMENSION + i] < 0.5);
    }
  }

  if (nfft_K_2_R.flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&nfft_K_2_R);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::initialize_nfft_R_2_k() {
  nfft_init(&nfft_R_2_k, DIMENSION, &grid_size[0], target_k_dmn_t::dmn_size());

  std::vector<std::vector<double>>& collection_k_vecs = target_k_dmn_t::get_elements();

  std::vector<std::vector<double>> collection_k_vecs_affine(target_k_dmn_t::dmn_size(),
                                                            std::vector<double>(DIMENSION, 0));

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++) {
    collection_k_vecs_affine[j] = math::util::coordinates(
        collection_k_vecs[j], source_k_cluster_type::get_super_basis_vectors());
  }

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++) {
    for (int i = 0; i < DIMENSION; i++) {
      while (collection_k_vecs_affine[j][i] < -1. / 2.)
        collection_k_vecs_affine[j][i] += 1.;

      while (collection_k_vecs_affine[j][i] > 1. / 2. - 1.e-6)
        collection_k_vecs_affine[j][i] -= 1.;
    }
  }

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++) {
    for (int i = 0; i < DIMENSION; i++) {
      nfft_R_2_k.x[j * DIMENSION + i] = collection_k_vecs_affine[j][i];

      assert(nfft_R_2_k.x[j * DIMENSION + i] >= -0.5 - 1.e-6);
      assert(nfft_R_2_k.x[j * DIMENSION + i] < 0.5);
    }
  }

  if (nfft_R_2_k.flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&nfft_R_2_k);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::FT_F_K__to__F_R(std::complex<double>* input_ptr) {
  assert(source_k_dmn_t::dmn_size() > 0);

  for (int K_ind = 0; K_ind < source_k_dmn_t::dmn_size(); K_ind++) {
    nfft_K_2_R.f[K_ind][0] = real(input_ptr[K_ind]);
    nfft_K_2_R.f[K_ind][1] = imag(input_ptr[K_ind]);
  }

  nfft_adjoint(&nfft_K_2_R);

  for (int R_ind = 0; R_ind < centered_r_cluster_dmn_t::dmn_size(); R_ind++) {
    if (lies_within_cutoff(R_ind) > 0) {
      F_R(R_ind)
          .real(nfft_K_2_R.f_hat[R_ind][0] /
                double(source_k_dmn_t::dmn_size() * lies_within_cutoff(R_ind)));
      F_R(R_ind)
          .imag(nfft_K_2_R.f_hat[R_ind][1] /
                double(source_k_dmn_t::dmn_size() * lies_within_cutoff(R_ind)));
    }
    else
      F_R(R_ind) = 0.;
  }
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::FT_F_R__to__F_k(std::complex<double>* output_ptr) {
  for (int r_ind = 0; r_ind < centered_r_cluster_dmn_t::dmn_size(); r_ind++) {
    nfft_R_2_k.f_hat[r_ind][0] = real(F_R(r_ind));
    nfft_R_2_k.f_hat[r_ind][1] = imag(F_R(r_ind));
  }

  nfft_trafo(&nfft_R_2_k);

  for (int k_ind = 0; k_ind < target_k_dmn_t::dmn_size(); k_ind++) {
    output_ptr[k_ind].real(nfft_R_2_k.f[k_ind][0]);
    output_ptr[k_ind].imag(nfft_R_2_k.f[k_ind][1]);
  }
}

/*!
 *  \brief the centered_r_cluster is in row-major order because of the FFTW, on which NFFT relies!
 */
template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::initialize_centered_r_cluster() {
  centered_r_cluster::get_elements().resize(0);

  switch (DIMENSION) {
    case 1: {
      for (int i0 = -grid_size[0] / 2; i0 < grid_size[0] / 2; i0++) {
        std::vector<double> r_vec(DIMENSION, 0);

        r_vec[0] = i0;

        centered_r_cluster::get_elements().push_back(r_vec);
      }
    } break;

    case 2: {
      for (int i1 = -grid_size[1] / 2; i1 < grid_size[1] / 2; i1++) {
        for (int i0 = -grid_size[0] / 2; i0 < grid_size[0] / 2; i0++) {
          std::vector<double> r_vec(DIMENSION, 0);

          r_vec[0] = i1;
          r_vec[1] = i0;

          centered_r_cluster::get_elements().push_back(r_vec);
        }
      }
    } break;

    case 3: {
      for (int i2 = -grid_size[2] / 2; i2 < grid_size[2] / 2; i2++) {
        for (int i1 = -grid_size[1] / 2; i1 < grid_size[1] / 2; i1++) {
          for (int i0 = -grid_size[0] / 2; i0 < grid_size[0] / 2; i0++) {
            std::vector<double> r_vec(DIMENSION, 0);

            r_vec[0] = i2;
            r_vec[1] = i1;
            r_vec[2] = i0;

            centered_r_cluster::get_elements().push_back(r_vec);
          }
        }
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  assert(int(centered_r_cluster::get_elements().size()) == centered_r_cluster::get_size());
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>,
                                  target_dmn_type>::initialize_cut_off() {
  for (int l = 0; l < lies_within_cutoff.size(); l++)
    lies_within_cutoff(l) = 0;

  for (int R_ind = 0; R_ind < source_r_dmn_t::dmn_size(); R_ind++) {
    std::vector<double> R_vec = source_r_dmn_t::get_elements()[R_ind];

    std::vector<std::vector<double>> r_vecs = cluster_operations::equivalent_vectors(
        R_vec, source_r_cluster_type::get_super_basis_vectors());

    for (size_t l = 0; l < r_vecs.size(); l++) {
      std::vector<double> r_aff =
          math::util::coordinates(r_vecs[l], source_r_cluster_type::get_basis_vectors());

      for (int R_cen_ind = 0; R_cen_ind < centered_r_cluster_dmn_t::dmn_size(); R_cen_ind++)
        if (math::util::distance2(r_aff, centered_r_cluster_dmn_t::get_elements()[R_cen_ind]) < 1.e-3)
          lies_within_cutoff(R_cen_ind) += r_vecs.size();
    }
  }

  if (DIMENSION == 2) {
    for (int l = 0; l < lies_within_cutoff.size(); l++) {
      if (centered_r_cluster_dmn_t::get_elements()[l][0] == -grid_size[1] / 2 or
          centered_r_cluster_dmn_t::get_elements()[l][1] == -grid_size[0] / 2)
        lies_within_cutoff(l) = 0;
    }
  }

  if (DIMENSION == 3) {
    for (int l = 0; l < lies_within_cutoff.size(); l++) {
      if (centered_r_cluster_dmn_t::get_elements()[l][0] == -grid_size[2] / 2 or
          centered_r_cluster_dmn_t::get_elements()[l][1] == -grid_size[1] / 2 or
          centered_r_cluster_dmn_t::get_elements()[l][2] == -grid_size[0] / 2)
        lies_within_cutoff(l) = 0;
    }
  }
}

// Template specialization for source and target domains that include spin-orbital subdomains.
template <typename b_dmn_t, typename source_k_dmn_type, typename target_k_dmn_type>
class wannier_interpolation_kernel<func::dmn_variadic<b_dmn_t, b_dmn_t, source_k_dmn_type>,
                                   func::dmn_variadic<b_dmn_t, b_dmn_t, target_k_dmn_type>> {
  typedef typename source_k_dmn_type::parameter_type source_k_cluster_type;
  typedef typename target_k_dmn_type::parameter_type target_k_cluster_type;

  // typedef typename source_k_cluster_type::base_cluster source_base_cluster_type;
  // typedef typename target_k_cluster_type::base_cluster target_base_cluster_type;

  typedef wannier_interpolation_kernel<source_k_cluster_type, target_k_cluster_type>
      wannier_interpolation_kernel_type;
  // typedef wannier_interpolation_kernel<source_base_cluster_type, target_k_cluster_type>
  // wannier_interpolation_kernel_type;

  typedef typename wannier_interpolation_kernel_type::centered_r_cluster_dmn_t centered_r_dmn_t;
  // typedef func::dmn_0<centered_r_LDA> centered_r_dmn_t;

  typedef func::dmn_variadic<b_dmn_t, b_dmn_t, source_k_dmn_type> input_dmn_t;
  typedef func::dmn_variadic<b_dmn_t, b_dmn_t, target_k_dmn_type> output_dmn_t;

public:
  wannier_interpolation_kernel();

  void reset();

  void reset_functions();

  void set(func::function<std::complex<double>, input_dmn_t>& H_K);
  void get(func::function<std::complex<double>, output_dmn_t>& H_k);

private:
  b_dmn_t dmn;

  wannier_interpolation_kernel_type wannier_kernel_object;

  // func::function<std::complex<double>, centered_r_dmn_t>  in;
  func::function<std::complex<double>, source_k_dmn_type> in;
  func::function<std::complex<double>, target_k_dmn_type> out;

  func::function<std::complex<double>, func::dmn_variadic<b_dmn_t, b_dmn_t, centered_r_dmn_t>> F_R;
};

template <typename b_dmn_t, typename source_k_dmn_type, typename target_k_dmn_type>
wannier_interpolation_kernel<
    func::dmn_variadic<b_dmn_t, b_dmn_t, source_k_dmn_type>,
    func::dmn_variadic<b_dmn_t, b_dmn_t, target_k_dmn_type>>::wannier_interpolation_kernel()
    : wannier_kernel_object(),

      F_R("wannier_interpolation_kernel__F_r") {}

template <typename b_dmn_t, typename source_k_dmn_type, typename target_k_dmn_type>
void wannier_interpolation_kernel<func::dmn_variadic<b_dmn_t, b_dmn_t, source_k_dmn_type>,
                                  func::dmn_variadic<b_dmn_t, b_dmn_t, target_k_dmn_type>>::reset() {
  wannier_kernel_object.reset_output();
}

template <typename b_dmn_t, typename source_k_dmn_type, typename target_k_dmn_type>
void wannier_interpolation_kernel<
    func::dmn_variadic<b_dmn_t, b_dmn_t, source_k_dmn_type>,
    func::dmn_variadic<b_dmn_t, b_dmn_t, target_k_dmn_type>>::reset_functions() {
  in.reset();
  out.reset();

  F_R.reset();
}

template <typename b_dmn_t, typename source_k_dmn_type, typename target_k_dmn_type>
void wannier_interpolation_kernel<func::dmn_variadic<b_dmn_t, b_dmn_t, source_k_dmn_type>,
                                  func::dmn_variadic<b_dmn_t, b_dmn_t, target_k_dmn_type>>::
    set(func::function<std::complex<double>, input_dmn_t>& H_K) {
  for (int i = 0; i < dmn.get_size(); ++i) {
    for (int j = 0; j < dmn.get_size(); ++j) {
      for (int k = 0; k < source_k_cluster_type::get_size(); ++k)
        in(k) = H_K(i, j, k);

      wannier_kernel_object.set(&in(0));

      for (int r = 0; r < centered_r_dmn_t::dmn_size(); ++r)
        F_R(i, j, r) = wannier_kernel_object.get_F_r(r);
    }
  }
}

template <typename b_dmn_t, typename source_k_dmn_type, typename target_k_dmn_type>
void wannier_interpolation_kernel<func::dmn_variadic<b_dmn_t, b_dmn_t, source_k_dmn_type>,
                                  func::dmn_variadic<b_dmn_t, b_dmn_t, target_k_dmn_type>>::
    get(func::function<std::complex<double>, output_dmn_t>& H_k) {
  for (int i = 0; i < dmn.get_size(); ++i) {
    for (int j = 0; j < dmn.get_size(); ++j) {
      for (int r = 0; r < centered_r_dmn_t::dmn_size(); ++r)
        wannier_kernel_object.get_F_r(r) = F_R(i, j, r);

      wannier_kernel_object.get(&out(0));

      for (int k = 0; k < target_k_cluster_type::get_size(); ++k)
        H_k(i, j, k) = out(k);
    }
  }
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_KERNEL_HPP
