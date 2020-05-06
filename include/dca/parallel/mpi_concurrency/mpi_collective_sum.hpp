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

#include <algorithm>  // std::min
#include <numeric>  // std::partial_sum
#include <map>
#include <string>
#include <utility>  // std::move, std::swap
#include <vector>

#include <mpi.h>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/buffer.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPICollectiveSum : public virtual MPIProcessorGrouping {
public:
  MPICollectiveSum() = default;

  // Wrappers to MPI_Allreduce.
  template <typename scalar_type>
  void sum(scalar_type& value) const;
  template <typename scalar_type>
  void sum(std::vector<scalar_type>& m) const;
  template <typename scalartype>
  void sum(std::map<std::string, std::vector<scalartype>>& m);
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

  // Wrapper to MPI_Reduce.
  template <typename Scalar, class Domain>
  void localSum(func::function<Scalar, Domain>& f, int root_id) const;

  // Wrapper to MPI_Reduce.
  template <typename Scalar, class Domain>
  void gatherv(func::function<Scalar, Domain>& f, int root_id) const;

  // Delay the execution of sum (implemented with MPI_Allreduce) until 'resolveSums' is called,
  // or 'delayedSum' is called with an object of different Scalar type.
  template <typename Scalar>
  void delayedSum(Scalar& obj);
  template <typename Scalar>
  void delayedSum(std::vector<Scalar>& f);
  template <typename Scalar, class domain>
  void delayedSum(func::function<Scalar, domain>& f);

  // Execute all the reductions scheduled with 'delayedSum' or delayed leaveOneOutSum out calls.
  void resolveSums();

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

  // Computes the sum of s over all ranks excluding the local value and stores the result back
  // in s.
  // Does nothing, if there is only one rank.
  // In/Out: s
  // In: delay. If true delay the sum until 'resolveSums' is called.
  template <typename T>
  void leaveOneOutSum(T& s, bool delay = 0);

  // Element-wise implementations for dca::func::function.
  // In/Out: f
  // In: delay. If true delay the sum until 'resolveSums' is called.
  template <typename Scalar, class Domain>
  void leaveOneOutSum(func::function<Scalar, Domain>& f, bool delay = 0);

  // Computes the average of s over all ranks excluding the local value and stores the result back
  // in s.
  // Does nothing, if there is only one rank.
  // In/Out: s
  template <typename T>
  void leaveOneOutAvg(T& s);

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

  // Computes the covariance of f with respect to its sample average. Then moves the average into f.
  // In/Out: f
  // Out: cov
  template <typename ScalarOrComplex, typename Scalar, class Domain, class CovDomain>
  void computeCovarianceAndAvg(func::function<ScalarOrComplex, Domain>& f,
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
  // Compute the sum on process 'rank_id', or all processes if rank_id == -1.
  template <typename T>
  void sum(const T* in, T* out, std::size_t n, int rank_id = -1) const;

  template <typename T>
  void gatherv_helper(const T* in, T* out, std::size_t total_size, int root_id = 0) const;

  template <typename T>
  void delayedSum(T* in, std::size_t n);
  template <typename T>
  void delayedSum(std::complex<T>* in, std::size_t n);

  template <typename T>
  void resolveSumsImplementation();

  // Objects for delayed sums and "leave one out" sums.
  MPI_Datatype current_type_ = 0;
  io::Buffer packed_;
  // Each vector entry stores the properties of a single delayed sum.
  std::vector<char*> pointers_;
  std::vector<std::size_t> sizes_;  // number of scalars involved in each operation, not bytes.
  std::vector<std::int8_t> remove_local_value_;  // use int8_t to avoid vector<bool> issues.
};

template <typename scalar_type>
void MPICollectiveSum::sum(scalar_type& value) const {
  scalar_type result;

  MPI_Allreduce(&value, &result, 1, MPITypeMap<scalar_type>::value(), MPI_SUM,
                MPIProcessorGrouping::get());

  value = result;
}

template <typename scalar_type>
void MPICollectiveSum::sum(std::vector<scalar_type>& m) const {
  std::vector<scalar_type> result(m.size());

  sum(m.data(), result.data(), m.size());

  m = std::move(result);
}

template <typename Scalar>
void MPICollectiveSum::sum(std::map<std::string, std::vector<Scalar>>& m) {
  for (auto it = m.begin(); it != m.end(); ++it) {
    delayedSum((it->second));
  }

  resolveSums();
}

template <typename scalar_type, class domain>
void MPICollectiveSum::sum(func::function<scalar_type, domain>& f) const {
  func::function<scalar_type, domain> f_sum;

  sum(f, f_sum);

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
  sum(f_in.values(), f_out.values(), f_in.size());
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

  sum(vec.ptr(), vec_sum.ptr(), vec.size());

  vec = std::move(vec_sum);

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

  sum(f.ptr(), F.ptr(), Nr * Nc);

  f = std::move(F);
}

template <typename scalar_type, class domain>
void MPICollectiveSum::localSum(func::function<scalar_type, domain>& f, int id) const {
  if (id < 0 || id > get_size())
    throw(std::out_of_range("id out of range."));

  func::function<scalar_type, domain> f_sum;

  sum(f.values(), f_sum.values(), f.size(), id);

  f = std::move(f_sum);
}

template <typename scalar_type, class domain>
void MPICollectiveSum::gatherv(func::function<scalar_type, domain>& f, int id) const {
    if (id < 0 || id > get_size())
        throw(std::out_of_range("id out of range."));

    func::function<scalar_type, domain> f_sum;

    int my_rank, mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(my_rank == 0)
    {
        std::cout << "\n *********performing MPI Gatherv Sum************* \n";
    }

    my_rank == 0 ? gatherv_helper(f.values(), f.values(), f.get_domain().get_size(), id)
        : gatherv_helper(f.values(), f_sum.values(), f.get_domain().get_size(), id);
}

template <typename some_type>
void MPICollectiveSum::sum_and_average(some_type& obj, const int nr_meas_rank) const {
  sum(obj);

  const double one_over_N = 1. / (nr_meas_rank * MPIProcessorGrouping::get_size());

  obj *= one_over_N;
}

template <typename some_type>
void MPICollectiveSum::sum_and_average(const some_type& in, some_type& out,
                                       const int nr_meas_rank) const {
  sum(in, out);

  const double one_over_N = 1. / (nr_meas_rank * MPIProcessorGrouping::get_size());

  out *= one_over_N;
}

template <typename Scalar>
void MPICollectiveSum::leaveOneOutSum(Scalar& s, bool delay) {
  if (MPIProcessorGrouping::get_size() == 1)
    return;

  const Scalar s_local(s);
  if (delay == false) {
    sum(s);
    s = s - s_local;
  }
  else {
    delayedSum(s);
    remove_local_value_.back() = true;
  }
}

template <typename Scalar, class Domain>
void MPICollectiveSum::leaveOneOutSum(func::function<Scalar, Domain>& f, const bool delay) {
  if (MPIProcessorGrouping::get_size() == 1)
    return;

  if (delay == false) {
    const func::function<Scalar, Domain> f_local(f);
    sum(f_local, f);
    for (int i = 0; i < f.size(); ++i)
      f(i) = f(i) - f_local(i);
  }
  else {
    delayedSum(f);
    remove_local_value_.back() = false;
  }
}

template <typename T>
void MPICollectiveSum::leaveOneOutAvg(T& x) {
  if (MPIProcessorGrouping::get_size() == 1)
    return;

  leaveOneOutSum(x);
  x /= MPIProcessorGrouping::get_size() - 1;
}

template <typename Scalar, class Domain>
func::function<Scalar, Domain> MPICollectiveSum::jackknifeError(func::function<Scalar, Domain>& f_i,
                                                                const bool overwrite) const {
  func::function<Scalar, Domain> err("jackknife-error");

  const int n = MPIProcessorGrouping::get_size();

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

  const int n = MPIProcessorGrouping::get_size();

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
  scalar_type factor = 1. / MPIProcessorGrouping::get_size();

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

  f_stddev /= std::sqrt(MPIProcessorGrouping::get_size());
}

template <typename scalar_type, class domain>
void MPICollectiveSum::average_and_compute_stddev(
    func::function<std::complex<scalar_type>, domain>& f_mean,
    func::function<std::complex<scalar_type>, domain>& f_stddev) const {
  scalar_type factor = 1. / MPIProcessorGrouping::get_size();

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

  f_stddev /= std::sqrt(MPIProcessorGrouping::get_size());
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

template <typename ScalarOrComplex, typename Scalar, class Domain, class CovDomain>
void MPICollectiveSum::computeCovarianceAndAvg(func::function<ScalarOrComplex, Domain>& f,
                                               func::function<Scalar, CovDomain>& cov) const {
  auto f_avg = f;
  sum_and_average(f_avg);
  computeCovariance(f, f_avg, cov);
  f = std::move(f_avg);
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
  const int n = MPIProcessorGrouping::get_size();
  for (std::size_t i = 0; i < f.size(); i++) {
    const Scalar var = std::sqrt(var2[i] / n);
    for (std::size_t j = 0; j < orders.size(); j++)
      momenta_avg[j] += std::abs(momenta(j, i)) / (n * std::pow(var, orders[j]));
  }

  for (std::size_t j = 0; j < orders.size(); j++)
    momenta_avg[j] /= f.size();

  return momenta_avg;
}

template <typename T>
void MPICollectiveSum::sum(const T* in, T* out, std::size_t n, int root_id) const {
  // On summit large messages hangs if sizeof(floating point type) * message_size > 2^31-1.
  constexpr std::size_t max_size = dca::util::IsComplex<T>::value
                                       ? 2 * (std::numeric_limits<int>::max() / sizeof(T))
                                       : std::numeric_limits<int>::max() / sizeof(T);

  for (std::size_t start = 0; start < n; start += max_size) {
    const int msg_size = std::min(n - start, max_size);
    if (root_id == -1) {
      MPI_Allreduce(in + start, out + start, msg_size, MPITypeMap<T>::value(), MPI_SUM,
                    MPIProcessorGrouping::get());
    }
    else {
      MPI_Reduce(in + start, out + start, msg_size, MPITypeMap<T>::value(), MPI_SUM, root_id,
                 MPIProcessorGrouping::get());
    }
  }
}

template <typename T>
void MPICollectiveSum::gatherv_helper(const T* in, T* out, std::size_t total_size, int root_id) const {
    // On summit large messages hangs if sizeof(floating point type) * message_size > 2^31-1.
//    constexpr std::size_t max_size = dca::util::IsComplex<T>::value
//                                     ? 2 * (std::numeric_limits<int>::max() / sizeof(T))
//                                     : std::numeric_limits<int>::max() / sizeof(T);

//    for (std::size_t start = 0; start < n; start += max_size) {
//        const int msg_size = std::min(n - start, max_size);
//        std::cout << "\n\n msg_size is " << msg_size
//                    << " max_size is " << max_size
//                    << " dca::util::IsComplex<T>::value is " << dca::util::IsComplex<T>::value
//                    << " (std::numeric_limits<int>::max() / sizeof(T)) is " << std::numeric_limits<int>::max() / sizeof(T)
//                    << "\n\n";
//        if (root_id == -1) {
//            MPI_Allreduce(in + start, out + start, msg_size, MPITypeMap<T>::value(), MPI_SUM,
//                          MPIProcessorGrouping::get());
//        }
//        else {
//            MPI_Reduce(in + start, out + MPIProcessorGrouping::get_id() * msg_size, msg_size, MPITypeMap<T>::value(), MPI_SUM, root_id,
//                       MPIProcessorGrouping::get());

    int mpi_size = MPIProcessorGrouping::get_size();
    int my_rank = MPIProcessorGrouping::get_id();

    uint64_t local_work = total_size / mpi_size;
    uint64_t more_work_before_index;

    std::vector<int> ranks_workload(mpi_size, 0);
    std::vector<int> displs(mpi_size + 1, 0);
    int* p_ranks_workload = ranks_workload.data();
    int* p_displs = displs.data();

    std::fill(ranks_workload.begin(), ranks_workload.end(), local_work);

    bool balanced = (total_size % mpi_size) == 0 ? true : false;

    if(balanced)
    {
        std::partial_sum(ranks_workload.begin(), ranks_workload.end(), displs.begin() + 1, std::plus<int>());
        displs.pop_back();
    }
    else
    {
        more_work_before_index = total_size % mpi_size;
        std::transform(ranks_workload.begin(), ranks_workload.begin() + more_work_before_index-1,
                       ranks_workload.begin(), [](int ele){ return ele+1; });
        std::partial_sum(ranks_workload.begin(), ranks_workload.end(), displs.begin() + 1, std::plus<int>());
        displs.pop_back();
    }


    MPI_Gatherv(in, ranks_workload[my_rank], MPITypeMap<T>::value(), out, p_ranks_workload, p_displs, MPITypeMap<T>::value(),
                root_id, MPIProcessorGrouping::get());
}

template <typename Scalar>
void MPICollectiveSum::delayedSum(Scalar& obj) {
  delayedSum(&obj, 1);
}

template <typename Scalar>
void MPICollectiveSum::delayedSum(std::vector<Scalar>& v) {
  delayedSum(v.data(), v.size());
}

template <typename Scalar, class Domain>
void MPICollectiveSum::delayedSum(func::function<Scalar, Domain>& f) {
  delayedSum(f.values(), f.size());
}

template <typename T>
void MPICollectiveSum::delayedSum(std::complex<T>* in, std::size_t n) {
  delayedSum(reinterpret_cast<T*>(in), n * 2);
}

template <typename T>
void MPICollectiveSum::delayedSum(T* in, std::size_t n) {
  if (current_type_ != 0 && current_type_ != MPITypeMap<T>::value())
    resolveSums();

  current_type_ = MPITypeMap<T>::value();
  pointers_.push_back(reinterpret_cast<char*>(in));
  sizes_.push_back(n);
  remove_local_value_.push_back(false);
  packed_.write(in, n);
}

// TODO: move to .cpp file.
inline void MPICollectiveSum::resolveSums() {
  // Note: we can not use a switch statement as MPI_Datatype is not integer on the Summit system.
  if (current_type_ == MPI_DOUBLE)
    return resolveSumsImplementation<double>();
  else if (current_type_ == MPI_FLOAT)
    return resolveSumsImplementation<float>();
  else if (current_type_ == MPI_UNSIGNED_LONG)
    return resolveSumsImplementation<unsigned long int>();
  else
    throw(std::logic_error("Type not supported."));
}

template <typename T>
inline void MPICollectiveSum::resolveSumsImplementation() {
  io::Buffer output(packed_.size());

  sum(reinterpret_cast<T*>(packed_.data()), reinterpret_cast<T*>(output.data()),
      packed_.size() / sizeof(T));

  for (int i = 0; i < pointers_.size(); ++i) {
    const std::size_t position = output.tellg();
    output.read(pointers_[i], sizes_[i] * sizeof(T));

    if (remove_local_value_[i]) {
      // Subtract local value from sum.
      T* sum = reinterpret_cast<T*>(pointers_[i]);
      const T* local = reinterpret_cast<T*>(packed_.data() + position);
      // sum[0:n] = sum[0:n] - local[0:n]
      std::transform(sum, sum + sizes_[i], local, sum, std::minus<T>());
    }
  }

  current_type_ = 0;
  sizes_.clear();
  remove_local_value_.clear();
  pointers_.clear();
  packed_.clear();
}

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_SUM_HPP
