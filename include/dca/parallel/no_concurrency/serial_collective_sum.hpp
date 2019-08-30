// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
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

#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

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
  void sum(func::function<Scalar, Domain>&) const {}
  template <typename Scalar, class Domain>
  void sum(func::function<Scalar, Domain>& /*f_in*/, func::function<Scalar, Domain>& /*f_out*/) const {
  }
  template <typename Scalar, class Domain>
  void sum(func::function<std::vector<Scalar>, Domain>&) const {}
  template <typename Scalar>
  void sum(linalg::Vector<Scalar, linalg::CPU>&) const {}
  template <typename Scalar>
  void sum(dca::linalg::Matrix<Scalar, linalg::CPU>&) const {}

  template <typename T>
  void sum_and_average(T& obj, int measurements) const {
    obj *= 1. / measurements;
  }

  template <typename T>
  void sum_and_average(T& /*obj*/) const {}

  template <typename T>
  void leaveOneOutAvg(T&) const {}
  template <typename T>
  void leaveOneOutSum(T&) const {}

  template <typename Scalar, class Domain>
  func::function<Scalar, Domain> jackknifeError(func::function<Scalar, Domain>&,
                                                const bool = true) const {
    return func::function<Scalar, Domain>();
  }
  template <typename Scalar, class Domain>
  func::function<std::complex<Scalar>, Domain> jackknifeError(
      func::function<std::complex<Scalar>, Domain>&, const bool = true) const {
    return func::function<std::complex<Scalar>, Domain>();
  }

  template <typename Scalar, typename Domain>
  void average_and_compute_stddev(func::function<Scalar, Domain>& /*f_mean*/,
                                  func::function<Scalar, Domain>& f_stddev) const {
    f_stddev = Scalar(0);
  }

  template <typename ScalarOrComplex, typename Scalar, class Domain, class CovDomain>
  void computeCovarianceAndAvg(func::function<ScalarOrComplex, Domain>& /*f*/,
                               func::function<Scalar, CovDomain>& cov) const {
    cov = Scalar(0);
  }
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_SERIAL_COLLECTIVE_SUM_HPP
