// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::


std::pair<short, short> InteractionVertices::getInsertionIndices(double random) const {
  random *= total_weigth_;
  // search in reverse order.
  const auto it_to_vertex =
      std::upper_bound(cumulative_weigths_.rbegin(), cumulative_weigths_.rend(), random);
  const int index = cumulative_weigths_.rend() - it_to_vertex - 1;

  assert(index >= 0 && index < size());
  return std::make_pair(index, elements_[index].partner_id);
}

void InteractionVertices::insertElement(InteractionElement v) {
  const short id = size();
 // Find the partner and exchange id.
  const auto& nu = v.nu;
  const auto& r = v.r;
  if (nu[0] != nu[1] or nu[2] != nu[3] or r[0] != r[1] or r[2] != r[3]) {  // non density-denity
    InteractionElement partner{{r[1], r[0], r[3], r[2]}, {nu[1], nu[0], nu[3], nu[2]}, v.w};
    for (auto& elem : elements_)
      if(elem == partner) {
        elem.partner_id = id;
        v.partner_id = &elem - elements_.data();
        break;
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


}  // ctint
}  // solver
}  // phys
}  // dca
