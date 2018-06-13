// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Andrei Plamada (plamada@phys.ethz.ch)
//
//  This class performs a statistical test according to docs/sources/StatisticalTesting.md.

#ifndef DCA_MATH_STATISTICAL_TESTING_STATISTICAL_TESTING_HPP
#define DCA_MATH_STATISTICAL_TESTING_STATISTICAL_TESTING_HPP

#include <algorithm>
#include <string>
#include <vector>

#include "dca/function/function.hpp"

namespace dca {
namespace math {
//  dca::math::

class StatisticalTesting {
public:
  // Creates the test and stores f: the observed average, f_expected, and the covariance of f.
  // In: f, f_expected, covariance
  template <typename Domain>
  StatisticalTesting(const func::function<double, Domain>& f,
                     const func::function<double, Domain>& f_expected,
                     const func::function<double, dca::func::dmn_variadic<Domain, Domain>>& covariance,
                     bool verbose = 0);

  // Only the selected indices are tested. Indices will be reordered and pruned of duplicates.
  // Note: The indices associated with each degree of freedom will be reset. Therefore multiple
  // calls with the same input remove different degrees of freedom.
  // In/Out: indices
  // Precondition: all elements of indices are in [0, get_dof()).
  void selectIndices(std::vector<int>& indices);

  // All the selected indices expect the selected ones are tested. Indices will be reordered and
  // pruned of duplicates.
  // Note: The indices associated with each degree of freedom will be reset. Therefore multiple
  // calls with the same input different degrees of freedom.
  // In/Out: indices
  // Precondition: all elements of indices are in [0, get_dof()).
  void discardIndices(std::vector<int>& indices);

  // Performs the test. The returned pvalue is the probability of obtaining
  // more extreme data then the observation. The stored covariance is destroyed.
  // In: known_covariance : true if the covariance is from a reference run, false if it is computed
  //                        with the empirical data.
  // In: n_samples        : number of gaussian samples that were averaged to obtain f.
  // In: allow_fast       : If true it attempts to invert the covariance in the original base
  //                        and the computation of the normalized samples is skipped.
  // Returns              : The p-value
  double computePValue(bool known_covariance, int n_samples, bool allow_fast = false);

  // Prints to a file the pvalue and the normalized samples, for further analysis.
  void printInfo(const std::string& filename, bool append = false) const;
  // Returns the number of degrees of freedom used in the test.
  int get_dof() const {
    return dof_;
  }

private:
  void computeMahalanobisDistanceSquared();
  void computeFastMahalanobisDistanceSquared();

  std::vector<double> df_;
  std::vector<double> cov_;
  bool verbose_;
  std::vector<double> normalized_samples_;
  int dof_ = -1;
  double distance_ = -1;
  double pvalue_ = -1;
  int samples_ = -1;
};

// Cumulative f distribution
double fCdf(double x, int nu1, int nu2);
// Cumulative chi squared distribution
double chi2Cdf(double x, int k);

// Implementation of templated methods.

template <typename Domain>
StatisticalTesting::StatisticalTesting(
    const func::function<double, Domain>& f, const func::function<double, Domain>& f_expected,
    const func::function<double, dca::func::dmn_variadic<Domain, Domain>>& covariance, bool verbose)
    : df_(f.size()), cov_(covariance.size()), verbose_(verbose) {
  for (auto i = 0; i < Domain::dmn_size(); ++i)
    df_[i] = f(i) - f_expected(i);
  std::copy_n(&covariance(0), covariance.size(), cov_.begin());

  dof_ = f.size();
}

}  // math
}  // dca

#endif  // DCA_MATH_STATISTICAL_TESTING_STATISTICAL_TESTING_HPP
