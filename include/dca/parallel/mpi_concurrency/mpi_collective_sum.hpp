// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class provides an interface to do collective sums and averages with MPI.
// In addition, it can compute the covariance for func::function.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_SUM_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_SUM_HPP

#include <map>
#include <string>
#include <utility>  // std::move, std::swap
#include <vector>

#include <mpi.h>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPICollectiveSum {
public:
  MPICollectiveSum(const MPIProcessorGrouping& grouping) : grouping_(grouping) {}

  template <typename scalar_type>
  void sum(scalar_type& value) const;
  template <typename scalar_type>
  void sum(std::vector<scalar_type>& m) const;
  template <typename scalartype>
  void sum(std::map<std::string, std::vector<scalartype>>& m) const;
  template <typename scalar_type, class domain>
  void sum(func::function<scalar_type, domain>& f) const;
  template <typename scalar_type, class domain>
  void sum(const func::function<scalar_type, domain>& f_in,
           func::function<scalar_type, domain>& f_out) const;
  template <typename scalar_type, class domain>
  void sum(func::function<std::vector<scalar_type>, domain>& f) const;
  template <typename scalar_type>
  void sum(linalg::Vector<scalar_type, linalg::CPU>& vec) const;
  template <typename scalar_type>
  void sum(linalg::Matrix<scalar_type, linalg::CPU>& f) const;

  template <typename some_type>
  void sum_and_average(some_type& obj, int nr_meas_rank = 1) const;
  template <typename some_type>
  void sum_and_average(const some_type& in, some_type& out, int nr_meas_rank = 1) const;

  template <typename scalar_type, class domain>
  void average_and_compute_stddev(func::function<scalar_type, domain>& f_mean,
                                  func::function<scalar_type, domain>& f_stddev) const;
  template <typename scalar_type, class domain>
  void average_and_compute_stddev(func::function<std::complex<scalar_type>, domain>& f_mean,
                                  func::function<std::complex<scalar_type>, domain>& f_stddev) const;

  // Overwrites the values of f with its average across all other ranks.
  // Does nothing if there is only one rank.
  template <typename scalar_type, class domain>
  void leaveOneOutAvg(func::function<scalar_type, domain>& f) const;
  template <typename scalar_type>
  void leaveOneOutAvg(scalar_type& s) const;

  // Compute the jack knife variance as sqrt( (n-1)/n \sum_i (f_i - f_avg)^2 ), where f_i is the
  // function value in the current rank.
  // I/O : f_loo
  // Returns: The jack knife error.
  // Precondition: f_loo contains the estimation with one block of data
  // left out.
  // Postcondition: if overwite_average == true, f_loo contains the average of the L.O.O.
  // estimations.
  template <typename scalar_type, class domain>
  func::function<scalar_type, domain> jackKnifeError(func::function<scalar_type, domain>& f_loo,
                                                     bool overwrite_average = true) const;
  // Same as above. Computes the error for real and imaginary part independently.
  template <typename scalar_type, class domain>
  func::function<std::complex<scalar_type>, domain> jackKnifeError(
      func::function<std::complex<scalar_type>, domain>& f_loo, bool overwrite_average = true) const;

  // Computes the covariance matrix of the measurements of the different mpi ranks.
  // In: f, f_estimated
  // Out: cov
  template <typename Scalar, class Domain>
  void computeCovariance(const func::function<Scalar, Domain>& f,
                         const func::function<Scalar, Domain>& f_estimated,
                         func::function<Scalar, func::dmn_variadic<Domain, Domain>>& cov) const;

  // Computes the covariance matrix of complex measurements of the different mpi ranks.
  // The real part and the imaginary part are treated independently, and cov represents
  // the covariance of the real vector [Re(f[0]), Re(f[1]), ..., Im(f[0]), Im(f[1]), ...].
  // In: f, f_estimated
  // Out: cov
  template <typename Scalar, class Domain, class CovDomain>
  void computeCovariance(const func::function<std::complex<Scalar>, Domain>& f,
                         const func::function<std::complex<Scalar>, Domain>& f_estimated,
                         func::function<Scalar, CovDomain>& cov) const;

private:
  const MPIProcessorGrouping& grouping_;
};

template <typename scalar_type>
void MPICollectiveSum::sum(scalar_type& value) const {
  scalar_type result;

  MPI_Allreduce(&value, &result, MPITypeMap<scalar_type>::factor(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_.get());

  value = result;
}

template <typename scalar_type>
void MPICollectiveSum::sum(std::vector<scalar_type>& m) const {
  std::vector<scalar_type> result(m.size(), scalar_type(0));

  MPI_Allreduce(&(m[0]), &(result[0]), MPITypeMap<scalar_type>::factor() * m.size(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_.get());

  m = std::move(result);
}

template <typename scalar_type>
void MPICollectiveSum::sum(std::map<std::string, std::vector<scalar_type>>& m) const {
  typedef typename std::map<std::string, std::vector<scalar_type>>::iterator iterator_type;

  iterator_type it = m.begin();

  for (; it != m.end(); ++it) {
    std::vector<scalar_type> values((it->second).size());

    for (size_t l = 0; l < (it->second).size(); l++)
      values[l] = (it->second)[l];

    sum(values);

    for (size_t l = 0; l < (it->second).size(); l++)
      (it->second)[l] = values[l];
  }
}

template <typename scalar_type, class domain>
void MPICollectiveSum::sum(func::function<scalar_type, domain>& f) const {
  func::function<scalar_type, domain> f_sum;

  // TODO: Use f.values() instead of &f(0) (same for f_sum).
  MPI_Allreduce(&f(0), &f_sum(0), MPITypeMap<scalar_type>::factor() * f.size(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_.get());

  f = std::move(f_sum);

#ifndef NDEBUG
  for (int i = 0; i < f.size(); ++i) {
    // INTERNAL: Cannot use std::isnan since scalar_type might be std::complex.
    if (f(i) != f(i))
      throw std::logic_error("Summation resulted in not-a-number (NaN) value.");
  }
#endif  // NDEBUG
}

template <typename scalar_type, class domain>
void MPICollectiveSum::sum(const func::function<scalar_type, domain>& f_in,
                           func::function<scalar_type, domain>& f_out) const {
  MPI_Allreduce(&f_in(0), &f_out(0), MPITypeMap<scalar_type>::factor() * f_in.size(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_.get());
}

template <typename scalar_type, class domain>
void MPICollectiveSum::sum(func::function<std::vector<scalar_type>, domain>& f) const {
  int Nr = f(0).size();
  int Nc = f.size();

  linalg::Matrix<scalar_type, linalg::CPU> M("M", std::pair<int, int>(Nr, Nc));

  for (int j = 0; j < Nc; j++)
    for (int i = 0; i < Nr; i++)
      M(i, j) = f(j)[i];

  sum(M);

  for (int j = 0; j < Nc; j++)
    for (int i = 0; i < Nr; i++)
      f(j)[i] = M(i, j);
}

template <typename scalar_type>
void MPICollectiveSum::sum(linalg::Vector<scalar_type, linalg::CPU>& vec) const {
  linalg::Vector<scalar_type, linalg::CPU> vec_sum("vec_sum", vec.size());

  MPI_Allreduce(&vec[0], &vec_sum[0], MPITypeMap<scalar_type>::factor() * vec.size(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_.get());

  vec = vec_sum;

#ifndef NDEBUG
  for (int i = 0; i < vec.size(); ++i)
    if (vec[i] != vec[i])
      throw std::logic_error("Summation resulted in not-a-number (NaN) value.");
#endif  // NDEBUG
}

template <typename scalar_type>
void MPICollectiveSum::sum(linalg::Matrix<scalar_type, linalg::CPU>& f) const {
  linalg::Matrix<scalar_type, linalg::CPU> F("F", f.size(), f.capacity());

  assert(f.capacity().first == F.capacity().first);
  assert(f.capacity().second == F.capacity().second);

  int Nr = f.capacity().first;
  int Nc = f.capacity().second;

  MPI_Allreduce(&f(0, 0), &F(0, 0), MPITypeMap<scalar_type>::factor() * Nr * Nc,
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_.get());

  for (int j = 0; j < F.size().second; j++)
    for (int i = 0; i < F.size().first; i++)
      f(i, j) = F(i, j);
}

template <typename some_type>
void MPICollectiveSum::sum_and_average(some_type& obj, const int nr_meas_rank) const {
  sum(obj);

  const double one_over_N = 1. / (nr_meas_rank * grouping_.get_Nr_threads());

  obj *= one_over_N;
}

template <typename some_type>
void MPICollectiveSum::sum_and_average(const some_type& in, some_type& out,
                                       const int nr_meas_rank) const {
  sum(in, out);

  const double one_over_N = 1. / (nr_meas_rank * grouping_.get_Nr_threads());

  out *= one_over_N;
}

template <typename scalar_type>
void MPICollectiveSum::leaveOneOutAvg(scalar_type& s) const {
  if (grouping_.get_Nr_threads() == 1)  // no jack knife estimate is possible,
    return;

  scalar_type copy = s;
  sum(s);
  s = (s - copy) / (grouping_.get_Nr_threads() - 1);
}

template <typename scalar_type, class domain>
void MPICollectiveSum::leaveOneOutAvg(func::function<scalar_type, domain>& f) const {
  if (grouping_.get_Nr_threads() == 1)  // no jack knife estimate is possible,
    return;

  func::function<scalar_type, domain> f_copy(f);
  sum(f_copy, f);
  const double scale = 1. / (grouping_.get_Nr_threads() - 1);
  for (int i = 0; i < f.size(); i++)
    f(i) = scale * (f(i) - f_copy(i));
}

template <typename scalar_type, class domain>
func::function<scalar_type, domain> MPICollectiveSum::jackKnifeError(
    func::function<scalar_type, domain>& f_loo, bool overwrite_average) const {
  func::function<scalar_type, domain> f_avg(f_loo);
  func::function<scalar_type, domain> err("Jack-Knife-error");

  const int n = grouping_.get_Nr_threads();

  if (n == 1)  // no jack knife estimate is possible,
    return err;

  sum_and_average(f_avg);
  const double scale = double(n - 1) / double(n);
  for (int i = 0; i < f_avg.size(); i++)
    err(i) = (f_loo(i) - f_avg(i)) * (f_loo(i) - f_avg(i));

  sum(err);

  for (int i = 0; i < err.size(); i++)
    err(i) = std::sqrt(scale * err(i));

  if (overwrite_average)
    f_loo = std::move(f_avg);
  return err;
}

template <typename scalar_type, class domain>
func::function<std::complex<scalar_type>, domain> MPICollectiveSum::jackKnifeError(
    func::function<std::complex<scalar_type>, domain>& f_loo, bool overwrite_average) const {
  func::function<std::complex<scalar_type>, domain> f_avg(f_loo);
  func::function<std::complex<scalar_type>, domain> err("Jack-Knife-error");

  sum_and_average(f_avg);
  const int n = grouping_.get_Nr_threads();
  const double scale = double(n - 1) / double(n);
  for (int i = 0; i < f_avg.size(); i++) {
    err(i).real(real(f_loo(i) - f_avg(i)) * real(f_loo(i) - f_avg(i)));
    err(i).imag(imag(f_loo(i) - f_avg(i)) * imag(f_loo(i) - f_avg(i)));
  }
  sum(err);

  for (int i = 0; i < err.size(); i++) {
    err(i).real(std::sqrt(scale * std::real(err(i))));
    err(i).imag(std::sqrt(scale * std::imag(err(i))));
  }

  if (overwrite_average)
    f_loo = std::move(f_avg);
  return err;
}

template <typename scalar_type, class domain>
void MPICollectiveSum::average_and_compute_stddev(func::function<scalar_type, domain>& f_mean,
                                                  func::function<scalar_type, domain>& f_stddev) const {
  scalar_type factor = 1. / grouping_.get_Nr_threads();

  func::function<scalar_type, domain> f_sum("f-sum");
  func::function<scalar_type, domain> f_diff("f-diff");

  sum(f_mean, f_sum);

  f_sum *= factor;

  for (int i = 0; i < f_sum.size(); i++)
    f_diff(i) = (f_mean(i) - f_sum(i)) * (f_mean(i) - f_sum(i));

  for (int i = 0; i < f_sum.size(); i++)
    f_mean(i) = f_sum(i);

  sum(f_diff, f_stddev);

  f_stddev *= factor;

  for (int i = 0; i < f_sum.size(); i++)
    f_stddev(i) = std::sqrt(f_stddev(i));

  f_stddev /= std::sqrt(grouping_.get_Nr_threads());
}

template <typename scalar_type, class domain>
void MPICollectiveSum::average_and_compute_stddev(
    func::function<std::complex<scalar_type>, domain>& f_mean,
    func::function<std::complex<scalar_type>, domain>& f_stddev) const {
  scalar_type factor = 1. / grouping_.get_Nr_threads();

  func::function<std::complex<scalar_type>, domain> f_sum("f-sum");
  func::function<std::complex<scalar_type>, domain> f_diff("f-diff");

  sum(f_mean, f_sum);

  f_sum *= factor;

  for (int i = 0; i < f_sum.size(); i++) {
    f_diff(i).real(real(f_mean(i) - f_sum(i)) * real(f_mean(i) - f_sum(i)));
    f_diff(i).imag(imag(f_mean(i) - f_sum(i)) * imag(f_mean(i) - f_sum(i)));
  }

  for (int i = 0; i < f_sum.size(); i++)
    f_mean(i) = f_sum(i);

  sum(f_diff, f_stddev);

  f_stddev *= factor;

  for (int i = 0; i < f_sum.size(); i++) {
    f_stddev(i).real(std::sqrt(real(f_stddev(i))));
    f_stddev(i).imag(std::sqrt(imag(f_stddev(i))));
  }

  f_stddev /= std::sqrt(grouping_.get_Nr_threads());
}

template <typename Scalar, class Domain>
void MPICollectiveSum::computeCovariance(
    const func::function<Scalar, Domain>& f, const func::function<Scalar, Domain>& f_estimated,
    func::function<Scalar, func::dmn_variadic<Domain, Domain>>& cov) const {
  for (int i = 0; i < f.size(); i++)
    for (int j = 0; j < f.size(); j++)
      cov(i, j) = (f(i) - f_estimated(i)) * (f(j) - f_estimated(j));

  sum_and_average(cov, 1);
}

template <typename Scalar, class Domain, class CovDomain>
void MPICollectiveSum::computeCovariance(const func::function<std::complex<Scalar>, Domain>& f,
                                         const func::function<std::complex<Scalar>, Domain>& f_estimated,
                                         func::function<Scalar, CovDomain>& cov) const {
  assert(4 * f.size() * f.size() == cov.size());

  // Treat real and imaginary parts as independent entries
  for (int i = 0; i < f.size(); i++)
    for (int j = 0; j < f.size(); j++) {
      // Real - Real
      cov(i, j) = (f(i).real() - f_estimated(i).real()) * (f(j).real() - f_estimated(j).real());

      // Imaginary - Imaginary
      cov(i + f.size(), j + f.size()) =
          (f(i).imag() - f_estimated(i).imag()) * (f(j).imag() - f_estimated(j).imag());

      // Real - Imaginary
      cov(i, j + f.size()) =
          (f(i).real() - f_estimated(i).real()) * (f(j).imag() - f_estimated(j).imag());

      // Imaginary - Real
      cov(i + f.size(), j) =
          (f(i).imag() - f_estimated(i).imag()) * (f(j).real() - f_estimated(j).imag());
    }
  sum_and_average(cov, 1);
}

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_SUM_HPP
