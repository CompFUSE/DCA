// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides a (fake) interface to do "collective" sums and averages in serial execution.
// All the summing methods are empty since summing the result of ONE process, of course,  doesn't
// change the result. Only the sum_and_average method does real work by averaging over all the
// measurements.

#ifndef DCA_PARALLEL_NO_CONCURRENCY_SERIAL_COLLECTIVE_SUM_HPP
#define DCA_PARALLEL_NO_CONCURRENCY_SERIAL_COLLECTIVE_SUM_HPP

#include <map>
#include <string>
#include <vector>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"
#include "comp_library/function_library/function.h"

namespace dca {
namespace parallel {
// dca::parallel::

class SerialCollectiveSum {
public:
  template <typename Scalar>
  void sum(Scalar&) const {}
  template <typename Scalar>
  void sum(std::vector<Scalar>&) const {}
  template <typename Scalar>
  void sum(std::map<std::string, std::vector<Scalar>>&) const {}
  template <typename Scalar, class Domain>
  void sum(FUNC_LIB::function<Scalar, Domain>&) const {}
  template <typename Scalar, class Domain>
  void sum(FUNC_LIB::function<Scalar, Domain>& /*f_in*/,
           FUNC_LIB::function<Scalar, Domain>& /*f_out*/) const {}
  template <typename Scalar, class Domain>
  void sum(FUNC_LIB::function<std::vector<Scalar>, Domain>&) const {}
  template <typename Scalar>
  void sum(linalg::Vector<Scalar, linalg::CPU>&) const {}
  template <typename Scalar>
  void sum(dca::linalg::Matrix<Scalar, linalg::CPU>&) const {}

  template <typename T>
  void sum_and_average(T& obj, int measurements) const {
    obj *= 1. / measurements;
  }
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_SERIAL_COLLECTIVE_SUM_HPP
