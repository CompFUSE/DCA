// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class represent all possible interaction terms of the hamiltonian, and provide
// functionalities to randomly sample from them and pair non density-density terms together.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_INTERACTION_VERTICES_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_INTERACTION_VERTICES_HPP

#include <algorithm>
#include <array>
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
  std::array<unsigned short, 4> r;
  std::array<unsigned short, 4> nu;
  double w;
  // TODO: write proper constructor.
  std::vector<unsigned short> partners_id = std::vector<unsigned short>();
  // Returns true if r and nu members are equal.
  bool operator==(const InteractionElement& other) const;
};

// Store the interaction terms and allow drawing a random vertex with strength |w|.
class InteractionVertices {
public:
  void initialize(double double_insertion_prob, bool all_sites_partnership);

  // Initialize vertices with a density-density Hamiltonian.
  // Precondition: Domain has the shape of dmn_variadic<Nu, Nu, Rdmn>
  template <class Nu, class Rdmn>
  void initializeFromHamiltonian(const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn>>& H_int);
  template <class Nu, class Rdmn>
  void initializeFromNonDensityHamiltonian(
      const func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>& H_int);

  void insertElement(InteractionElement v);

  void insertElement(const std::vector<double>& vec);

  template <class Nu, class Rdmn, class TDmn>
  void checkForInterbandPropagators(
      const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn, TDmn>>& G_r_t);

  void reset();
  // In: random number generator object.
  // Returns: first: random vertex sampled with probability proportional to |vertex.w|.
  //          second: first vertex partner's id if it exists, -1 otherwise.
  template <class Rng>
  std::pair<short, short> getInsertionIndices(Rng& rng, double double_update_prob) const;

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
  int possiblePartners(unsigned idx) const {
    assert(idx < elements_.size());
    return elements_[idx].partners_id.size();
  }

  std::vector<InteractionElement> elements_;

private:
  enum PartnershipType { NONE, SAME_SITE, ALL_SITES };

  std::vector<double> cumulative_weigths_;
  double total_weigth_ = 0;
  bool interband_propagator_ = false;
  PartnershipType partnership_type_ = NONE;
};

template <class Rng>
std::pair<short, short> InteractionVertices::getInsertionIndices(Rng& rng,
                                                                 double double_update_prob) const {
  const double random = rng() * total_weigth_;
  // search in reverse order.
  const auto it_to_vertex =
      std::upper_bound(cumulative_weigths_.rbegin(), cumulative_weigths_.rend(), random);
  const int index = cumulative_weigths_.rend() - it_to_vertex - 1;
  assert(index >= 0 && index < size());

  assert(double_update_prob >= 0 && double_update_prob <= 1);
  auto do_double = [&]() -> bool {
    if (double_update_prob == 0)
      return 0;
    else if (double_update_prob == 1)
      return 1;
    else
      return rng() < double_update_prob;
  };

  if (elements_[index].partners_id.size() && do_double()) {  // double insertion
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
    const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn>>& H_int) {
  // Note: Use this version of the code if double counting in the interaction hamiltonian is removed.
  //
  //  for (unsigned short nu1 = 0; nu1 < Nu::dmn_size(); nu1++) {
  //    for (unsigned short nu2 = 0; nu2 < Nu::dmn_size(); nu2++)
  //      for (unsigned short delta_r = 0; delta_r < Rdmn::dmn_size(); delta_r++) {
  //        const double value = H_int(nu1, nu2, delta_r);
  //        if (std::abs(value) < 1e-8)
  //          continue;
  //        for (unsigned short r1 = 0; r1 < Rdmn::dmn_size(); r1++) {
  //          const unsigned short r2 = Rdmn::parameter_type::subtract(delta_r, r1); // delta_r = r1 - r2
  //          insertElement(InteractionElement{{r1, r1, r2, r2}, {nu1, nu1, nu2, nu2}, value, short(-1)});
  //        }
  //      }
  //  }

  // Assume the density-density interaction Hamiltonian function is double counted, i.e.
  // H(b1, b2, r1 - r2) == H(b2, b1, r -r2) and both terms describe the same addendum to the
  // physical Hamiltonian.
  func::function<char, func::dmn_variadic<Nu, Nu, Rdmn>> already_inserted;
  const int r0 = Rdmn::parameter_type::origin_index();
  for (unsigned short nu1 = 0; nu1 < Nu::dmn_size(); nu1++) {
    for (unsigned short nu2 = 0; nu2 < Nu::dmn_size(); nu2++)
      for (unsigned short delta_r = 0; delta_r < Rdmn::dmn_size(); delta_r++) {
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
        for (unsigned short r1 = 0; r1 < Rdmn::dmn_size(); r1++) {
          const unsigned short r2 =
              Rdmn::parameter_type::subtract(delta_r, r1);  // delta_r = r1 - r2
          insertElement(InteractionElement{{r1, r1, r2, r2}, {nu1, nu1, nu2, nu2}, value});
        }
      }
  }
}

template <class Nu, class Rdmn>
void InteractionVertices::initializeFromNonDensityHamiltonian(
    const func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>& H_int) {
  auto spin = [](unsigned short nu) { return nu >= Nu::dmn_size() / 2; };
  auto check_spins = [&](unsigned short nu1, unsigned short nu2, unsigned short nu3,
                         unsigned short nu4) -> bool {
    return spin(nu1) == spin(nu2) && spin(nu3) == spin(nu4);
  };

  for (unsigned short nu1 = 0; nu1 < Nu::dmn_size(); nu1++)
    for (unsigned short nu2 = 0; nu2 < Nu::dmn_size(); nu2++)
      for (unsigned short nu3 = 0; nu3 < Nu::dmn_size(); nu3++)
        for (unsigned short nu4 = 0; nu4 < Nu::dmn_size(); nu4++)
          for (unsigned short delta_r = 0; delta_r < Rdmn::dmn_size(); delta_r++) {
            const double value = H_int(nu1, nu2, nu3, nu4, delta_r);
            if (std::abs(value) < 1e-8)
              continue;
            if (not check_spins(nu1, nu2, nu3, nu4))
              throw(std::logic_error(
                  "Format input hamiltonian s.t. pair of creation and annihilation "
                  "operators have same spin."));
            for (unsigned short r1 = 0; r1 < Rdmn::dmn_size(); r1++) {
              const unsigned short r2 = Rdmn::parameter_type::subtract(delta_r, r1);
              insertElement(InteractionElement{{r1, r1, r2, r2}, {nu1, nu2, nu3, nu4}, value});
            }
          }
}

template <class Nu, class RDmn, class TDmn>
void InteractionVertices::checkForInterbandPropagators(
    const func::function<double, func::dmn_variadic<Nu, Nu, RDmn, TDmn>>& G_r_t) {
  interband_propagator_ = false;
  const int t0 = TDmn::dmn_size() / 2;
  const int nb = Nu::dmn_size() / 2;
  const int r0 = RDmn::parameter_type::origin_index();
  for (int r = 0; r < RDmn::dmn_size(); ++r)
    for (int b1 = 0; b1 < nb; ++b1)
      for (int b2 = 0; b2 < nb; ++b2) {
        if (r != r0 && b1 != b2 && std::abs(G_r_t(b1, 0, b2, 0, r, t0)) > 1e-5) {
          interband_propagator_ = true;
          return;
        }
      }
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_INTERACTION_VERTICES_HPP
