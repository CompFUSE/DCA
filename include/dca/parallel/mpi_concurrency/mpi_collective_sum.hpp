// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
#include <memory>
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
  MPICollectiveSum(const std::unique_ptr<const MPIProcessorGrouping>& grouping)
      : grouping_(grouping) {}

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

  // Computes the average of s over all ranks excluding the local value and stores the result back
  // in s.
  // Does nothing, if there is only one rank.
  // In/Out: s
  template <typename Scalar>
  void leaveOneOutAvg(Scalar& s) const;
  // Element-wise implementation for dca::func::function.
  // In/Out: f
  template <typename Scalar, class Domain>
  void leaveOneOutAvg(func::function<Scalar, Domain>& f) const;

  // Computes and returns the element-wise jackknife error
  // \Delta f_{jack} = \sqrt( (n-1)/n \sum_i^n (f_i - f_avg)^2 ),
  // where f_i is the i-th jackknife estimate of f and f_avg is the average of all f_i's.
  // Returns zero for all elements, if there is only one rank.
  // In/Out: f_i
  // In (optional): overwrite (default = true)
  // Preconditions: Each rank holds a unique precomputed jackknife estimate f_i.
  // Postconditions: If overwrite == true, the jackknife estimate f_i is overwritten with the
  //                 average f_avg; else f_i stays unchanged.
  template <typename Scalar, class Domain>
  func::function<Scalar, Domain> jackknifeError(func::function<Scalar, Domain>& f_i,
                                                bool overwrite = true) const;
  // Implementation for std::complex.
  // Real and imaginary parts are treated independently.
  template <typename Scalar, class Domain>
  func::function<std::complex<Scalar>, Domain> jackknifeError(
      func::function<std::complex<Scalar>, Domain>& f_i, bool overwrite = true) const;

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

  // Returns the normalized momenta of the desired orders. The momenta is the averaged across the
  // entries of f.
  // The normalized momenta of order i is defined as \lambda_i = <(f - <f>)^i> / \sigma^i
  // In: f, orders
  // Returns: momenta
  template <typename Scalar, class Domain>
  std::vector<Scalar> avgNormalizedMomenta(const func::function<Scalar, Domain>& f,
                                           const std::vector<int>& orders) const;

private:
  const std::unique_ptr<const MPIProcessorGrouping>& grouping_;
};

template <typename scalar_type>
void MPICollectiveSum::sum(scalar_type& value) const {
  scalar_type result;

  MPI_Allreduce(&value, &result, MPITypeMap<scalar_type>::factor(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_->get());

  value = result;
}

template <typename scalar_type>
void MPICollectiveSum::sum(std::vector<scalar_type>& m) const {
  std::vector<scalar_type> result(m.size(), scalar_type(0));

  MPI_Allreduce(&(m[0]), &(result[0]), MPITypeMap<scalar_type>::factor() * m.size(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_->get());

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

  MPI_Allreduce(f.values(), f_sum.values(), MPITypeMap<scalar_type>::factor() * f.size(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_->get());

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
  MPI_Allreduce(f_in.values(), f_out.values(), MPITypeMap<scalar_type>::factor() * f_in.size(),
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_->get());
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
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_->get());

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
                MPITypeMap<scalar_type>::value(), MPI_SUM, grouping_->get());

  for (int j = 0; j < F.size().second; j++)
    for (int i = 0; i < F.size().first; i++)
      f(i, j) = F(i, j);
}

template <typename some_type>
void MPICollectiveSum::sum_and_average(some_type& obj, const int nr_meas_rank) const {
  sum(obj);

  const double one_over_N = 1. / (nr_meas_rank * grouping_->get_size());

  obj *= one_over_N;
}

template <typename some_type>
void MPICollectiveSum::sum_and_average(const some_type& in, some_type& out,
                                       const int nr_meas_rank) const {
  sum(in, out);

  const double one_over_N = 1. / (nr_meas_rank * grouping_->get_size());

  out *= one_over_N;
}

template <typename Scalar>
void MPICollectiveSum::leaveOneOutAvg(Scalar& s) const {
  if (grouping_->get_size() == 1)
    return;

  const Scalar s_local(s);
  sum(s);
  s = (s - s_local) / (grouping_->get_size() - 1);
}

template <typename Scalar, class Domain>
void MPICollectiveSum::leaveOneOutAvg(func::function<Scalar, Domain>& f) const {
  if (grouping_->get_size() == 1)
    return;

  const func::function<Scalar, Domain> f_local(f);
  sum(f_local, f);

  const double scale = 1. / (grouping_->get_size() - 1);
  for (int i = 0; i < f.size(); ++i)
    f(i) = scale * (f(i) - f_local(i));
}

template <typename Scalar, class Domain>
func::function<Scalar, Domain> MPICollectiveSum::jackknifeError(func::function<Scalar, Domain>& f_i,
                                                                const bool overwrite) const {
  func::function<Scalar, Domain> err("jackknife-error");

  const int n = grouping_->get_size();

  if (n == 1)  // No jackknife procedure possible.
    return err;

  func::function<Scalar, Domain> f_avg("f_avg");
  sum_and_average(f_i, f_avg);

  for (int k = 0; k < f_avg.size(); ++k)
    err(k) = (f_i(k) - f_avg(k)) * (f_i(k) - f_avg(k));

  sum(err);

  const double scale = double(n - 1) / double(n);
  for (int k = 0; k < err.size(); ++k)
    err(k) = std::sqrt(scale * err(k));

  if (overwrite)
    f_i = std::move(f_avg);

  return err;
}

template <typename Scalar, class Domain>
func::function<std::complex<Scalar>, Domain> MPICollectiveSum::jackknifeError(
    func::function<std::complex<Scalar>, Domain>& f_i, const bool overwrite) const {
  func::function<std::complex<Scalar>, Domain> err("jackknife-error");

  const int n = grouping_->get_size();

  if (n == 1)  // No jackknife procedure possible.
    return err;

  func::function<std::complex<Scalar>, Domain> f_avg("f_avg");
  sum_and_average(f_i, f_avg);

  for (int k = 0; k < f_avg.size(); ++k) {
    const auto diff = f_i(k) - f_avg(k);
    err(k).real(diff.real() * diff.real());
    err(k).imag(diff.imag() * diff.imag());
  }

  sum(err);

  const double scale = double(n - 1) / double(n);
  for (int k = 0; k < err.size(); ++k) {
    err(k).real(std::sqrt(scale * err(k).real()));
    err(k).imag(std::sqrt(scale * err(k).imag()));
  }

  if (overwrite)
    f_i = std::move(f_avg);

  return err;
}

template <typename scalar_type, class domain>
void MPICollectiveSum::average_and_compute_stddev(func::function<scalar_type, domain>& f_mean,
                                                  func::function<scalar_type, domain>& f_stddev) const {
  scalar_type factor = 1. / grouping_->get_size();

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

  f_stddev /= std::sqrt(grouping_->get_size());
}

template <typename scalar_type, class domain>
void MPICollectiveSum::average_and_compute_stddev(
    func::function<std::complex<scalar_type>, domain>& f_mean,
    func::function<std::complex<scalar_type>, domain>& f_stddev) const {
  scalar_type factor = 1. / grouping_->get_size();

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

  f_stddev /= std::sqrt(grouping_->get_size());
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

template <typename Scalar, class Domain>
std::vector<Scalar> MPICollectiveSum::avgNormalizedMomenta(const func::function<Scalar, Domain>& f,
                                                           const std::vector<int>& orders) const {
  func::function<Scalar, Domain> avg(f);
  sum_and_average(avg);
  linalg::Matrix<Scalar, linalg::CPU> momenta(std::make_pair(orders.size(), f.size()));
  std::vector<Scalar> var2(f.size());

  for (std::size_t i = 0; i < f.size(); i++) {
    const Scalar diff = f(i) - avg(i);
    var2[i] = diff * diff;
    for (std::size_t j = 0; j < orders.size(); j++)
      momenta(j, i) = std::pow(diff, orders[j]);
  }

  sum(momenta);
  sum(var2);

  // Divide by n and normalize the momenta by sigma^order, then average.
  std::vector<Scalar> momenta_avg(orders.size(), 0);
  const int n = grouping_->get_size();
  for (std::size_t i = 0; i < f.size(); i++) {
    const Scalar var = std::sqrt(var2[i] / n);
    for (std::size_t j = 0; j < orders.size(); j++)
      momenta_avg[j] += std::abs(momenta(j, i)) / (n * std::pow(var, orders[j]));
  }

  for (std::size_t j = 0; j < orders.size(); j++)
    momenta_avg[j] /= f.size();

  return momenta_avg;
}

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_SUM_HPP
