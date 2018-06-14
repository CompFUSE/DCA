// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class provides a (fake) interface to do "collective" gather.
// All the methods are empty.

#ifndef DCA_PARALLEL_NO_CONCURRENCY_SERIAL_COLLECTIVE_GATHER_HPP
#define DCA_PARALLEL_NO_CONCURRENCY_SERIAL_COLLECTIVE_GATHER_HPP

#include <vector>

#include "dca/function/function.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class SerialCollectiveGather {
public:
  template <typename Scalar, class Domain>
  void gather(func::function<Scalar, Domain>& /*f*/) const {}

  template <typename Scalar>
  void gather(std::vector<Scalar>& /*v*/) const {}
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_SERIAL_COLLECTIVE_GATHER_HPP
