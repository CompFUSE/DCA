// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_INTERACTION_VERTICES_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_INTERACTION_VERTICES_HPP

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "dca/function/function.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

// Represent a matrix element W(i,j,k,l) of the interaction hamiltonian.
// Each index represent a cluster position and a spin-band index.
struct InteractionElement {
  std::array<ushort, 4> r;
  std::array<ushort, 4> nu;
  double w;
  // TODO: write proper constructor.
  std::vector<ushort> partners_id = std::vector<ushort>();
  // Returns true if r and nu members are equal.
  bool operator==(const InteractionElement& other) const;
};

// Store the interaction terms and allow drawing a random vertex with strength |w|.
class InteractionVertices {
public:
  // Initialize vertices with a density-density Hamiltonian.
  // Precondition: Domain has the shape of dmn_variadic<Nu, Nu, Rdmn>
  template <class Nu, class Rdmn>
  void initializeFromHamiltonian(const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn>>& H_int,
                                 bool double_counted = true);
  template <class Nu, class Rdmn>
  void initializeFromNonDensityHamiltonian(
      const func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>& H_int);

  void insertElement(InteractionElement v);

  void insertElement(const std::vector<double>& vec);

  void reset();
  // In: random number generator object.
  // Returns: first: random vertex sampled with probability proportional to |vertex.w|.
  //          second: first vertex partner's id if it exists, -1 otherwise.
  template <class Rng>
  std::pair<short, short> getInsertionIndices(Rng& rng, bool double_update) const;

  // Returns: the sum of the absolute values of the interaction strengths.
  double integratedInteraction() const {
    return total_weigth_;
  }

  std::size_t size() const {
    return elements_.size();
  }
  const InteractionElement& operator[](const int idx) const {
    assert(idx < int(size()));
    return elements_[idx];
  }

  // Returns the number of possible partners for each non density-density interaction.
  int possiblePartners() const{
      const int partners = elements_.back().partners_id.size();
      assert(partners > 1);
      return partners;
  }

  std::vector<InteractionElement> elements_;

private:
  std::vector<double> cumulative_weigths_;
  double total_weigth_ = 0;
};

template <class Rng>
std::pair<short, short> InteractionVertices::getInsertionIndices(Rng& rng, bool double_update) const {
  const double random = rng() * total_weigth_;
  // search in reverse order.
  const auto it_to_vertex =
      std::upper_bound(cumulative_weigths_.rbegin(), cumulative_weigths_.rend(), random);
  const int index = cumulative_weigths_.rend() - it_to_vertex - 1;
  assert(index >= 0 && index < size());

  if (double_update && elements_[index].partners_id.size()) {  // double insertion
    const auto& partners = elements_[index].partners_id;
    auto partner_id = partners[rng() * partners.size()];
    return std::make_pair(index, partner_id);
  }
  else {
    return std::make_pair(index, -1);
  }
}

template <class Nu, class Rdmn>
void InteractionVertices::initializeFromHamiltonian(
    const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn>>& H_int, bool double_counted) {
  for (ushort nu1 = 0; nu1 < Nu::dmn_size(); nu1++) {
    const ushort nu2_start = double_counted ? nu1 + 1 : 0;
    for (ushort nu2 = nu2_start; nu2 < Nu::dmn_size(); nu2++)
      for (ushort delta_r = 0; delta_r < Rdmn::dmn_size(); delta_r++) {
        const double value = H_int(nu1, nu2, delta_r);
        if (std::abs(value) < 1e-8)
          continue;
        for (ushort r1 = 0; r1 < Rdmn::dmn_size(); r1++) {
          const ushort r2 = Rdmn::parameter_type::add(delta_r, r1);
          insertElement(InteractionElement{{r1, r1, r2, r2}, {nu1, nu1, nu2, nu2}, value});
        }
      }
  }
}

template <class Nu, class Rdmn>
void InteractionVertices::initializeFromNonDensityHamiltonian(
    const func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>& H_int) {
  auto spin = [](ushort nu) { return nu >= Nu::dmn_size() / 2; };
  auto check_spins = [&](ushort nu1, ushort nu2, ushort nu3, ushort nu4) -> bool {
    return spin(nu1) == spin(nu2) && spin(nu3) == spin(nu4);
  };

  for (ushort nu1 = 0; nu1 < Nu::dmn_size(); nu1++)
    for (ushort nu2 = 0; nu2 < Nu::dmn_size(); nu2++)
      for (ushort nu3 = 0; nu3 < Nu::dmn_size(); nu3++)
        for (ushort nu4 = 0; nu4 < Nu::dmn_size(); nu4++)
          for (ushort delta_r = 0; delta_r < Rdmn::dmn_size(); delta_r++) {
            const double value = H_int(nu1, nu2, nu3, nu4, delta_r);
            if (std::abs(value) < 1e-8)
              continue;
            if (not check_spins(nu1, nu2, nu3, nu4))
              throw(std::logic_error(
                  "Format input hamiltonian s.t. pair of creation and annihilation "
                  "operators have same spin."));
            for (ushort r1 = 0; r1 < Rdmn::dmn_size(); r1++) {
              const ushort r2 = Rdmn::parameter_type::add(delta_r, r1);
              insertElement(InteractionElement{{r1, r1, r2, r2}, {nu1, nu2, nu3, nu4}, value});
            }
          }
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_INTERACTION_VERTICES_HPP
