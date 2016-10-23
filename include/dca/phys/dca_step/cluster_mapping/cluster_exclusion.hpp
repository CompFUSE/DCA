// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_CLUSTER_EXCLUSION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_CLUSTER_EXCLUSION_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

#include "comp_library/function_plotting/include_plotting.h"
#include "comp_library/linalg/linalg.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename parameters_type, typename MOMS_type>
class cluster_exclusion {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using DCA_r_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::REAL_SPACE, domains::BRILLOUIN_ZONE>;
  using r_DCA = func::dmn_0<DCA_r_cluster_type>;
  using DCA_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;

public:
  cluster_exclusion(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  void execute();

private:
  void compute_G0_K_w_cluster_excluded();
  void compute_G0_R_t_cluster_excluded();

  void plot_G0_R_t_cluster_excluded();

private:
  parameters_type& parameters;
  MOMS_type& MOMS;
};

template <typename parameters_type, typename MOMS_type>
cluster_exclusion<parameters_type, MOMS_type>::cluster_exclusion(parameters_type& parameters_ref,
                                                                 MOMS_type& MOMS_ref)
    : parameters(parameters_ref), MOMS(MOMS_ref) {}

template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::execute() {
  compute_G0_K_w_cluster_excluded();

  compute_G0_R_t_cluster_excluded();
}

/*
 *   G = G_0 + G_0*S*G = G_0 * (1 + S*G)
 *
 *   G_0 = G*(1 + S*G)^-1
 */
template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::compute_G0_K_w_cluster_excluded() {
  profiler_type profiler(__FUNCTION__, "cluster_exclusion", __LINE__);

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix("G_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> S_matrix("S_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> one_plus_S_G("one_plus_S_G",
                                                                           nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix("G0_matrix", nu::dmn_size());

  // Allocate the work space for inverse only once.
  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int K_ind = 0; K_ind < k_DCA::dmn_size(); K_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_matrix(i, j) = MOMS.G_k_w(i, j, K_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          S_matrix(i, j) = MOMS.Sigma_cluster(i, j, K_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          one_plus_S_G(i, j) = 0;

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          for (int l = 0; l < nu::dmn_size(); l++)
            one_plus_S_G(i, j) += S_matrix(i, l) * G_matrix(l, j);

      for (int i = 0; i < nu::dmn_size(); i++)
        one_plus_S_G(i, i) += 1.;

      dca::linalg::matrixop::inverse(one_plus_S_G, ipiv, work);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = 0;

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          for (int l = 0; l < nu::dmn_size(); l++)
            G0_matrix(i, j) += G_matrix(i, l) * one_plus_S_G(l, j);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          MOMS.G0_k_w_cluster_excluded(i, j, K_ind, w_ind) = G0_matrix(i, j);
    }
  }
}

template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::compute_G0_R_t_cluster_excluded() {
  profiler_type profiler(__FUNCTION__, "cluster_exclusion", __LINE__);

  MOMS.G0_k_w_cluster_excluded -= MOMS.G0_k_w;

  {
    math::transform::FunctionTransform<w, t>::execute(MOMS.G0_k_w_cluster_excluded,
                                                      MOMS.G0_k_t_cluster_excluded);

    MOMS.G0_k_t_cluster_excluded += MOMS.G0_k_t;

    math::transform::FunctionTransform<k_DCA, r_DCA>::execute(MOMS.G0_k_t_cluster_excluded,
                                                              MOMS.G0_r_t_cluster_excluded);
  }

  MOMS.G0_k_w_cluster_excluded += MOMS.G0_k_w;
}

template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::plot_G0_R_t_cluster_excluded() {
  {
    func::function<float, t> tmp("G0_k_t");

    Gnuplot plot_obj("lines");
    for (int R_ind = 0; R_ind < r_DCA::dmn_size(); R_ind++) {
      for (int t_ind = 0; t_ind < t::dmn_size(); t_ind++)
        tmp(t_ind) = MOMS.G0_k_t(0, 0, R_ind, t_ind);

      SHOW::execute(plot_obj, tmp);
    }

    plot_obj.showonscreen();
  }

  {
    func::function<float, t> tmp("G0_k_t_cluster_excluded");

    Gnuplot plot_obj("lines");
    for (int R_ind = 0; R_ind < r_DCA::dmn_size(); R_ind++) {
      for (int t_ind = 0; t_ind < t::dmn_size(); t_ind++)
        tmp(t_ind) = MOMS.G0_k_t_cluster_excluded(0, 0, R_ind, t_ind);

      SHOW::execute(plot_obj, tmp);
    }

    plot_obj.showonscreen();
  }

  {
    func::function<float, t> tmp("G0_r_t");

    Gnuplot plot_obj("lines");
    for (int R_ind = 0; R_ind < r_DCA::dmn_size(); R_ind++) {
      for (int t_ind = 0; t_ind < t::dmn_size(); t_ind++)
        tmp(t_ind) = MOMS.G0_r_t_cluster_excluded(0, 0, R_ind, t_ind);

      SHOW::execute(plot_obj, tmp);
    }

    plot_obj.showonscreen();
  }
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_CLUSTER_EXCLUSION_HPP
