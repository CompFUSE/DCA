// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class symmetrizes Greens functions according to cluster symmetries, matsubara frequencies
// and band-index symmetries.

#ifndef DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_HPP
#define DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_HPP

#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize_single_particle_function.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize_two_particle_function.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <class Parameters>
class Symmetrize : public SymmetrizeSingleParticleFunction<Parameters>,
                   public symmetrize_two_particle_function {
public:
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;  // orbital-spin index
  using NuNuDmn = func::dmn_variadic<NuDmn, NuDmn>;

  template <typename Scalar, typename FDmn0>
  static void execute(func::function<Scalar, FDmn0>& f, bool do_diff = false);

  template <typename Scalar, typename FDmn0>
  static void execute(func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, FDmn0>>& f,
                      bool do_diff = false);

  template <typename Scalar, typename FDmn0, typename f_dmn_1>
  static void execute(func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, FDmn0, f_dmn_1>>& f,
                      bool do_diff = false);

  template <typename Scalar, typename FDmn0, typename f_dmn_1>
  static void execute(func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, FDmn0, f_dmn_1>>& f,
                      func::function<int, NuNuDmn>& H_symmetry, bool do_diff = false);

  template <typename Scalar, typename KDmn, typename WDmn>
  static void execute(
      func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>,
                                                func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename Scalar, typename KDmn, typename WDmn>
  static void execute(
      func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>,
                                                func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>>>& f,
      func::function<int, NuNuDmn>& H_symmetry, std::vector<double> Q, bool do_diff = false);

  template <typename Scalar, typename KDmn, typename WDmn>
  static void execute(
      func::function<Scalar, func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KDmn, KDmn, WDmn, WDmn>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename Scalar, typename KDmn, typename WDmn>
  static void execute(
      func::function<Scalar, func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KDmn, KDmn, WDmn, WDmn>>& f,
      func::function<int, NuNuDmn>& H_symmetry, std::vector<double> Q, bool do_diff = false);

  template <typename Scalar, typename KDmn, typename WDmn>
  static void execute(
      func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>,
                                                func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>, KDmn>>& f,
      bool do_diff = false);
};

template <class Parameters>
template <typename Scalar, typename FDmn0>
void Symmetrize<Parameters>::execute(func::function<Scalar, FDmn0>& f, bool do_diff) {
  SymmetrizeSingleParticleFunction<Parameters>::execute(f, do_diff);
}

template <class Parameters>
template <typename Scalar, typename FDmn0, typename f_dmn_1>
void Symmetrize<Parameters>::execute(
    func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, FDmn0, f_dmn_1>>& f,
    func::function<int, NuNuDmn>& H_symmetry, bool do_diff) {
  SymmetrizeSingleParticleFunction<Parameters>::execute(f, H_symmetry, do_diff);
}

template <class Parameters>
template <typename Scalar, typename FDmn0>
void Symmetrize<Parameters>::execute(func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, FDmn0>>& f,
                                     bool do_diff) {
  SymmetrizeSingleParticleFunction<Parameters>::execute(f, do_diff);
}

template <class Parameters>
template <typename Scalar, typename FDmn0, typename f_dmn_1>
void Symmetrize<Parameters>::execute(
    func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, FDmn0, f_dmn_1>>& f, bool do_diff) {
  SymmetrizeSingleParticleFunction<Parameters>::execute(f, do_diff);
}

template <class Parameters>
template <typename Scalar, typename KDmn, typename WDmn>
void Symmetrize<Parameters>::execute(
    func::function<Scalar, func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KDmn, KDmn, WDmn, WDmn>>& f,
    func::function<int, NuNuDmn>& /*H_symmetry*/, std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <class Parameters>
template <typename Scalar, typename KDmn, typename WDmn>
void Symmetrize<Parameters>::execute(
    func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>,
                                              func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>>>& f,
    func::function<int, NuNuDmn>& /*H_symmetry*/, std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <class Parameters>
template <typename Scalar, typename KDmn, typename WDmn>
void Symmetrize<Parameters>::execute(
    func::function<Scalar, func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KDmn, KDmn, WDmn, WDmn>>& f,
    std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <class Parameters>
template <typename Scalar, typename KDmn, typename WDmn>
void Symmetrize<Parameters>::execute(
    func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>,
                                              func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>>>& f,
    std::vector<double> Q, bool do_diff) {
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template <class Parameters>
template <typename Scalar, typename KDmn, typename WDmn>
void Symmetrize<Parameters>::execute(
    func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>,
                                              func::dmn_variadic<BDmn, BDmn, KDmn, WDmn>, KDmn>>& f,
    bool do_diff) {
  symmetrize_two_particle_function::execute(f, do_diff);
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_HPP
