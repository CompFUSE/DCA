// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the deconvolution step of the lattice mapping for single-particle
// functions.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_SP_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_SP_HPP

#include <algorithm>
#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/inference/richardson_lucy_deconvolution.hpp"
#include "dca/phys/dca_step/lattice_mapping/deconvolution/deconvolution_routines.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class deconvolution_sp
    : public deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t> {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

public:
  deconvolution_sp(parameters_type& parameters_ref);

  void execute(
      const func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn_t, w>>& f_source,
      const func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>&
          f_source_interpolated,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target_convoluted,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target);

private:
  void find_shift(func::function<std::complex<double>, func::dmn_variadic<b, b, w>>& shift,
                  const func::function<std::complex<double>,
                                       func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& Sigma_interp);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;
};

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::deconvolution_sp(
    parameters_type& parameters_ref)
    : deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>(parameters_ref),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()) {}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute(
    const func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn_t, w>>& f_source,
    const func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>&
        f_source_interpolated,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target_convoluted,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>& f_target) {
  typedef func::dmn_0<func::dmn<2, int>> z;
  typedef func::dmn_variadic<z, b, b, s, w> p_dmn_t;

  math::inference::RichardsonLucyDeconvolution<source_k_dmn_t, target_k_dmn_t, p_dmn_t> RL_obj(
      this->get_T_source_symmetrized(), this->get_T_symmetrized(),
      parameters.get_deconvolution_tolerance(), parameters.get_deconvolution_iterations());

  func::function<double, func::dmn_variadic<source_k_dmn_t, p_dmn_t>> source("source");
  func::function<double, func::dmn_variadic<target_k_dmn_t, p_dmn_t>> source_interpolated(
      "source_interpolated");
  func::function<double, func::dmn_variadic<target_k_dmn_t, p_dmn_t>> target_convoluted(
      "target_convoluted");
  func::function<double, func::dmn_variadic<target_k_dmn_t, p_dmn_t>> target("target");

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int k_ind = 0; k_ind < source_k_dmn_t::dmn_size(); k_ind++) {
      for (int j = 0; j < b::dmn_size(); j++) {
        for (int i = 0; i < b::dmn_size(); i++) {
          source(k_ind, 0, i, j, 0, w_ind) = real(f_source(i, 0, j, 0, k_ind, w_ind));
          source(k_ind, 1, i, j, 0, w_ind) = imag(f_source(i, 0, j, 0, k_ind, w_ind));
          source(k_ind, 0, i, j, 1, w_ind) = real(f_source(i, 1, j, 1, k_ind, w_ind));
          source(k_ind, 1, i, j, 1, w_ind) = imag(f_source(i, 1, j, 1, k_ind, w_ind));
        }
      }
    }
  }

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int k_ind = 0; k_ind < target_k_dmn_t::dmn_size(); k_ind++) {
      for (int j = 0; j < b::dmn_size(); j++) {
        for (int i = 0; i < b::dmn_size(); i++) {
          source_interpolated(k_ind, 0, i, j, 0, w_ind) =
              real(f_source_interpolated(i, 0, j, 0, k_ind, w_ind));
          source_interpolated(k_ind, 1, i, j, 0, w_ind) =
              imag(f_source_interpolated(i, 0, j, 0, k_ind, w_ind));
          source_interpolated(k_ind, 0, i, j, 1, w_ind) =
              real(f_source_interpolated(i, 1, j, 1, k_ind, w_ind));
          source_interpolated(k_ind, 1, i, j, 1, w_ind) =
              imag(f_source_interpolated(i, 1, j, 1, k_ind, w_ind));
        }
      }
    }
  }

  RL_obj.findTargetFunction(source, source_interpolated, target, target_convoluted,
                            concurrency.id() == concurrency.first());

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int k_ind = 0; k_ind < target_k_dmn_t::dmn_size(); k_ind++) {
      for (int j = 0; j < b::dmn_size(); j++) {
        for (int i = 0; i < b::dmn_size(); i++) {
          f_target_convoluted(i, 0, j, 0, k_ind, w_ind).real(target_convoluted(k_ind, 0, i, j, 0, w_ind));
          f_target_convoluted(i, 0, j, 0, k_ind, w_ind).imag(target_convoluted(k_ind, 1, i, j, 0, w_ind));
          f_target_convoluted(i, 1, j, 1, k_ind, w_ind).real(target_convoluted(k_ind, 0, i, j, 1, w_ind));
          f_target_convoluted(i, 1, j, 1, k_ind, w_ind).imag(target_convoluted(k_ind, 1, i, j, 1, w_ind));
        }
      }
    }
  }

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int k_ind = 0; k_ind < target_k_dmn_t::dmn_size(); k_ind++) {
      for (int j = 0; j < b::dmn_size(); j++) {
        for (int i = 0; i < b::dmn_size(); i++) {
          f_target(i, 0, j, 0, k_ind, w_ind).real(target(k_ind, 0, i, j, 0, w_ind));
          f_target(i, 0, j, 0, k_ind, w_ind).imag(target(k_ind, 1, i, j, 0, w_ind));
          f_target(i, 1, j, 1, k_ind, w_ind).real(target(k_ind, 0, i, j, 1, w_ind));
          f_target(i, 1, j, 1, k_ind, w_ind).imag(target(k_ind, 1, i, j, 1, w_ind));
        }
      }
    }
  }
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::find_shift(
    func::function<std::complex<double>, func::dmn_variadic<b, b, w>>& shift,
    const func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn_t, w>>&
        Sigma_interp) {
  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int j = 0; j < b::dmn_size(); j++) {
      for (int i = 0; i < b::dmn_size(); i++) {
        shift(i, j, w_ind).real(0);
        shift(i, j, w_ind).imag(0);
      }
    }
  }

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int k_ind = 0; k_ind < target_k_dmn_t::dmn_size(); k_ind++) {
      for (int j = 0; j < b::dmn_size(); j++) {
        for (int i = 0; i < b::dmn_size(); i++) {
          if (w_ind < w::dmn_size() / 2) {
            shift(i, j, w_ind)
                .real(std::min(real(shift(i, j, w_ind)), real(Sigma_interp(i, j, k_ind, w_ind))));
            shift(i, j, w_ind)
                .imag(std::min(imag(shift(i, j, w_ind)), imag(Sigma_interp(i, j, k_ind, w_ind))));
          }
          else {
            shift(i, j, w_ind)
                .real(std::max(real(shift(i, j, w_ind)), real(Sigma_interp(i, j, k_ind, w_ind))));
            shift(i, j, w_ind)
                .imag(std::max(imag(shift(i, j, w_ind)), imag(Sigma_interp(i, j, k_ind, w_ind))));
          }
        }
      }
    }
  }

  const double factor = 10.;

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int j = 0; j < b::dmn_size(); j++) {
      for (int i = 0; i < b::dmn_size(); i++) {
        if ((w_ind < w::dmn_size() / 2 && real(shift(i, j, w_ind)) < 0.) ||
            (w_ind >= w::dmn_size() / 2 && real(shift(i, j, w_ind)) > 0.)) {
          shift(i, j, w_ind).real(factor * real(shift(i, j, w_ind)));
        }
        else {
          shift(i, j, w_ind).real(0.);
        }

        if ((w_ind < w::dmn_size() / 2 && imag(shift(i, j, w_ind)) < 0.) ||
            (w_ind >= w::dmn_size() / 2 && imag(shift(i, j, w_ind)) > 0.)) {
          shift(i, j, w_ind).imag(factor * imag(shift(i, j, w_ind)));
        }
        else {
          shift(i, j, w_ind).imag(0.);
        }
      }
    }
  }
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_SP_HPP
