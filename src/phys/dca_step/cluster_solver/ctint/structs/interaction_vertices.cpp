// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

void InteractionVertices::initialize(double double_insertion_prob, bool all_sites_partnership) {
  if (double_insertion_prob <= 0) {
    partnership_type_ = NONE;
  }
  else if (all_sites_partnership == false) {
    partnership_type_ = SAME_SITE;
  }
  else {
    partnership_type_ = ALL_SITES;
  }
}

void InteractionVertices::insertElement(InteractionElement v) {
  const short id = size();
  // Find the partner and exchange id.
  const auto& nu = v.nu;
  const auto& r = v.r;
  if (partnership_type_ != NONE && (nu[0] != nu[1] or nu[2] != nu[3] or r[0] != r[1] or
                                    r[2] != r[3])) {  // non density-denity interaction
    const std::array<unsigned short, 4> nu_opposite{nu[1], nu[0], nu[3], nu[2]};
    const std::array<ushort, 4> r_opposite{r[1], r[0], r[3], r[2]};

    for (auto& elem : elements_)
      if (elem.nu == nu_opposite && (partnership_type_ == ALL_SITES || elem.r == r_opposite)) {
        elem.partners_id.push_back(id);
        v.partners_id.push_back(&elem - elements_.data());
      }
  }

  elements_.push_back(v);

  const double w = std::abs(v.w);
  total_weigth_ += w;
  cumulative_weigths_.push_back(0);
  for (auto& cum_weigth : cumulative_weigths_)
    cum_weigth += w;
}

void InteractionVertices::insertElement(const std::vector<double>& vec) {
  assert(vec.size() == 9);

  InteractionElement vert;
  for (int i = 0; i < 4; i++)
    vert.r[i] = static_cast<int>(vec[i]);
  for (int i = 0; i < 4; i++)
    vert.nu[i] = static_cast<int>(vec[i + 4]);
  vert.w = vec[8];
  insertElement(vert);
}

void InteractionVertices::reset() {
  elements_.clear();
  cumulative_weigths_.clear();
  total_weigth_ = 0;
}

bool InteractionElement::operator==(const InteractionElement& other) const {
  for (int i = 0; i < 4; ++i)
    if (r[i] != other.r[i] or nu[i] != other.nu[i])
      return false;
  return true;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
