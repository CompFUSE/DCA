// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
