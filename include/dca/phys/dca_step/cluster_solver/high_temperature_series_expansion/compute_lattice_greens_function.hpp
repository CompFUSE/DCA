// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the lattice Green's function.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_COMPUTE_LATTICE_GREENS_FUNCTION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_COMPUTE_LATTICE_GREENS_FUNCTION_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace htseries {
// dca::phys::solver::htseries::

template <class parameters_type, class MOMS_type, class k_dmn_t, class w_dmn_t>
class compute_lattice_Greens_function {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using k_HOST =
      func::dmn_0<domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::LATTICE_SP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  using function_type = func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST, w>>;

public:
  compute_lattice_Greens_function(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  function_type& get_G_k_w() {
    return G_k_w;
  }
  function_type& get_G0_k_w() {
    return G0_k_w;
  }

  void execute();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;
  MOMS_type& MOMS;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>> G_k_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>> G0_k_w;
};

template <class parameters_type, class MOMS_type, class k_dmn_t, class w_dmn_t>
compute_lattice_Greens_function<parameters_type, MOMS_type, k_dmn_t, w_dmn_t>::compute_lattice_Greens_function(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),
      MOMS(MOMS_ref),

      G_k_w("G_k_w"),
      G0_k_w("G0_k_w") {}

template <class parameters_type, class MOMS_type, class k_dmn_t, class w_dmn_t>
void compute_lattice_Greens_function<parameters_type, MOMS_type, k_dmn_t, w_dmn_t>::execute() {
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> I_k("I_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inv("G_inv", nu::dmn_size());

  // Allocate the work space for inverse only once.
  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    std::complex<double> i_wm_plus_mu;

    i_wm_plus_mu.real(parameters.get_chemical_potential());
    i_wm_plus_mu.imag(w::get_elements()[w_ind]);

    for (int i = 0; i < nu::dmn_size(); i++)
      I_k(i, i) = i_wm_plus_mu;

    for (int k_ind = 0; k_ind < k_HOST::dmn_size(); k_ind++) {
      {
        for (int j = 0; j < nu::dmn_size(); j++)
          for (int i = 0; i < nu::dmn_size(); i++)
            G_inv(i, j) = I_k(i, j) - MOMS.H_HOST(i, j, k_ind);

        dca::linalg::matrixop::inverse(G_inv, ipiv, work);

        for (int j = 0; j < nu::dmn_size(); j++)
          for (int i = 0; i < nu::dmn_size(); i++)
            G0_k_w(i, j, k_ind, w_ind) = G_inv(i, j);
      }

      {
        for (int j = 0; j < nu::dmn_size(); j++)
          for (int i = 0; i < nu::dmn_size(); i++)
            G_inv(i, j) =
                I_k(i, j) - MOMS.H_HOST(i, j, k_ind) - MOMS.Sigma_lattice(i, j, k_ind, w_ind);

        dca::linalg::matrixop::inverse(G_inv, ipiv, work);

        for (int j = 0; j < nu::dmn_size(); j++)
          for (int i = 0; i < nu::dmn_size(); i++)
            G_k_w(i, j, k_ind, w_ind) = G_inv(i, j);
      }
    }
  }
}

}  // htseries
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_COMPUTE_LATTICE_GREENS_FUNCTION_HPP
