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
  ushort r[4];
  ushort nu[4];
  double w;
  short partner_id = -1;
  // Returns true if r and nu members are equal.
  bool operator==(const InteractionElement& other) const;
};

// Store the interaction terms and allow drawing a random vertex with strength |w|.
class InteractionVertices {
public:
  // Initialize vertices with a density-density Hamiltonian.
  // Precondition: Domain has the shape of dmn_variadic<Nu, Nu, Rdmn>
  template <class Nu, class Rdmn>
  void initializeFromHamiltonian(const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn>>& H_int);
  template <class Nu, class Rdmn>
  void initializeFromNonDensityHamiltonian(
      const func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>& H_int);

  void insertElement(InteractionElement v);

  void insertElement(const std::vector<double>& vec);

  void reset();
  // In: random number generator object.
  // Returns: first: random vertex sampled with probability proportional to |vertex.w|.
  //          second: first vertex partner's id if it exists, -1 otherwise.
  std::pair<short, short> getInsertionIndices(double rand) const;
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

private:
  std::vector<InteractionElement> elements_;
  std::vector<double> cumulative_weigths_;
  double total_weigth_ = 0;
};

template <class Nu, class Rdmn>
void InteractionVertices::initializeFromHamiltonian(
    const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn>>& H_int) {
  // Note: Use this version of the code if double counting in the interaction hamiltonian is removed.
  //
  //  for (ushort nu1 = 0; nu1 < Nu::dmn_size(); nu1++) {
  //    for (ushort nu2 = 0; nu2 < Nu::dmn_size(); nu2++)
  //      for (ushort delta_r = 0; delta_r < Rdmn::dmn_size(); delta_r++) {
  //        const double value = H_int(nu1, nu2, delta_r);
  //        if (std::abs(value) < 1e-8)
  //          continue;
  //        for (ushort r1 = 0; r1 < Rdmn::dmn_size(); r1++) {
  //          const ushort r2 = Rdmn::parameter_type::subtract(delta_r, r1); // delta_r = r1 - r2
  //          insertElement(InteractionElement{{r1, r1, r2, r2}, {nu1, nu1, nu2, nu2}, value, short(-1)});
  //        }
  //      }
  //  }

  // Assume the density-density interaction Hamiltonian function is double counted, i.e.
  // H(b1, b2, r1 - r2) == H(b2, b1, r -r2) and both terms describe the same addendum to the
  // physical Hamiltonian.
  func::function<bool, func::dmn_variadic<Nu, Nu, Rdmn>> already_inserted;
  const int r0 = Rdmn::parameter_type::origin_index();
  for (ushort nu1 = 0; nu1 < Nu::dmn_size(); nu1++) {
    for (ushort nu2 = 0; nu2 < Nu::dmn_size(); nu2++)
      for (ushort delta_r = 0; delta_r < Rdmn::dmn_size(); delta_r++) {
        const double value = H_int(nu1, nu2, delta_r);
        if (std::abs(value) < 1e-8)
          continue;

        const double minus_delta_r = Rdmn::parameter_type::subtract(delta_r, r0);
        if (std::abs(H_int(nu1, nu2, delta_r) - H_int(nu2, nu1, minus_delta_r)) > 1e-8)
          throw(std::logic_error("The interaction double counting is not consistent."));

        if (already_inserted(nu2, nu1, minus_delta_r))  // Avoid double counting.
          continue;

        // Insert
        already_inserted(nu1, nu2, delta_r) = true;
        for (ushort r1 = 0; r1 < Rdmn::dmn_size(); r1++) {
          const ushort r2 = Rdmn::parameter_type::subtract(delta_r, r1);  // delta_r = r1 - r2
          insertElement(InteractionElement{{r1, r1, r2, r2}, {nu1, nu1, nu2, nu2}, value, short(-1)});
        }
      }
  }
}

template <class Nu, class Rdmn>
void InteractionVertices::initializeFromNonDensityHamiltonian(
    const func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>& H_int) {
  auto spin = [](ushort nu) { return nu >= Nu::dmn_size() / 2; };
  auto check_spins = [&](ushort nu1, ushort nu2, ushort nu3, ushort nu4) -> bool {
    return spin(nu1) == spin(nu2) and spin(nu3) == spin(nu4);
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
              insertElement(
                  InteractionElement{{r1, r1, r2, r2}, {nu1, nu2, nu3, nu4}, value, short(-1)});
            }
          }
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_INTERACTION_VERTICES_HPP
