// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@phys.ethz.ch)
//
// Description

#ifndef DCA_CONCURRENCY_INTERFACES_COLLECTIVE_SUM_INTERFACE_MPI_H
#define DCA_CONCURRENCY_INTERFACES_COLLECTIVE_SUM_INTERFACE_MPI_H

#include "dca/concurrency/interfaces/collective_sum_interface.h"
#include <mpi.h>
#include "dca/concurrency/interfaces/type_map_interface_mpi.h"
#include "dca/concurrency/interfaces/processor_grouping_interface_mpi.h"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <>
class collective_sum_interface<MPI_LIBRARY> {
public:
  collective_sum_interface(processor_grouping<MPI_LIBRARY>& grouping_ref);
  ~collective_sum_interface();

  template <typename scalar_type>
  void sum(scalar_type& value);

  template <typename scalar_type>
  void sum(std::vector<scalar_type>& m);

  template <typename scalartype>
  void sum(std::map<std::string, std::vector<scalartype>>& m);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<scalar_type, domain>& f);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<scalar_type, domain>& f,
           FUNC_LIB::function<scalar_type, domain>& f_target);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& f);

  template <typename scalar_type>
  void sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& f);

  template <typename scalar_type>
  void sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& f);

  template <typename some_type>
  void sum_and_average(some_type& obj, int size);

  template <typename scalar_type, class domain>
  void average_and_compute_stddev(FUNC_LIB::function<scalar_type, domain>& f_mean,
                                  FUNC_LIB::function<scalar_type, domain>& f_stddev, size_t size);

  template <typename scalar_type, class domain>
  void average_and_compute_stddev(FUNC_LIB::function<std::complex<scalar_type>, domain>& f_mean,
                                  FUNC_LIB::function<std::complex<scalar_type>, domain>& f_stddev,
                                  size_t size);

  // Compute the covariance matrix of the measurements of the different nodes.
  // In: f, f_estimated
  // Out: cov
  // TODO: const f, f_estimated
  template <typename Scalar, class domain>
  void computeCovariance(FUNC_LIB::function<Scalar, domain>& f,
                         FUNC_LIB::function<Scalar, domain>& f_estimated,
                         FUNC_LIB::function<Scalar, dmn_variadic<domain, domain>>& cov);

  // Compute the covariance matrix of the measurements of the different nodes.
  // The real part and the imaginary part are treated independently,
  // and cov represents the covariance of the vector [Re(f),Im(f)].
  // In: f, f_estimated
  // Out: cov
  // TODO: const f, f_estimated
  template <typename Scalar, class domain, class cov_domain>
  void computeCovariance(FUNC_LIB::function<std::complex<Scalar>, domain>& f,
                         FUNC_LIB::function<std::complex<Scalar>, domain>& f_estimated,
                         FUNC_LIB::function<Scalar, cov_domain>& cov);

private:
  processor_grouping<MPI_LIBRARY>& grouping;
};

collective_sum_interface<MPI_LIBRARY>::collective_sum_interface(
    processor_grouping<MPI_LIBRARY>& grouping_ref)
    : grouping(grouping_ref) {}

collective_sum_interface<MPI_LIBRARY>::~collective_sum_interface() {}

template <typename scalar_type>
void collective_sum_interface<MPI_LIBRARY>::sum(scalar_type& value) {
  scalar_type result;

  MPI_Allreduce(&value, &result, type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
                type_map_interface<MPI_LIBRARY, scalar_type>::value(), MPI_SUM, grouping.get());

  value = result;
}

template <typename scalar_type>
void collective_sum_interface<MPI_LIBRARY>::sum(std::vector<scalar_type>& m) {
  std::vector<scalar_type> result(m.size(), scalar_type(0));

  MPI_Allreduce(&(m[0]), &(result[0]),
                type_map_interface<MPI_LIBRARY, scalar_type>::factor() * m.size(),
                type_map_interface<MPI_LIBRARY, scalar_type>::value(), MPI_SUM, grouping.get());

  for (size_t i = 0; i < m.size(); i++)
    m[i] = result[i];
}

template <typename scalar_type>
void collective_sum_interface<MPI_LIBRARY>::sum(std::map<std::string, std::vector<scalar_type>>& m) {
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
void collective_sum_interface<MPI_LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& f) {
  FUNC_LIB::function<scalar_type, domain> F;

  MPI_Allreduce(&f(0), &F(0), type_map_interface<MPI_LIBRARY, scalar_type>::factor() * f.size(),
                type_map_interface<MPI_LIBRARY, scalar_type>::value(), MPI_SUM, grouping.get());

  for (int i = 0; i < F.size(); i++)
    f(i) = F(i);

  for (int i = 0; i < F.size(); i++) {
    if (f(i) != f(i)) {
      std::stringstream ss;
      ss << i << "\t" << f.get_name() << "\n";
      std::cout << ss.str();

      throw std::logic_error(__FUNCTION__);
    }
  }
}

template <typename scalar_type, class domain>
void collective_sum_interface<MPI_LIBRARY>::sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& f) {
  int Nr = f(0).size();
  int Nc = f.size();

  LIN_ALG::matrix<scalar_type, LIN_ALG::CPU> M("M", std::pair<int, int>(Nr, Nc));

  for (int j = 0; j < Nc; j++)
    for (int i = 0; i < Nr; i++)
      M(i, j) = f(j)[i];

  sum(M);

  for (int j = 0; j < Nc; j++)
    for (int i = 0; i < Nr; i++)
      f(j)[i] = M(i, j);
}

template <typename scalar_type, class domain>
void collective_sum_interface<MPI_LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& f,
                                                FUNC_LIB::function<scalar_type, domain>& f_target) {
  MPI_Allreduce(&f(0), &f_target(0),
                type_map_interface<MPI_LIBRARY, scalar_type>::factor() * f.size(),
                type_map_interface<MPI_LIBRARY, scalar_type>::value(), MPI_SUM, grouping.get());
}

template <typename scalar_type>
void collective_sum_interface<MPI_LIBRARY>::sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& f) {
  LIN_ALG::vector<scalar_type, LIN_ALG::CPU> F("F", f.get_current_size());

  MPI_Allreduce(&f[0], &F[0],
                type_map_interface<MPI_LIBRARY, scalar_type>::factor() * f.get_current_size(),
                type_map_interface<MPI_LIBRARY, scalar_type>::value(), MPI_SUM, grouping.get());

  for (int i = 0; i < F.get_current_size(); i++)
    f[i] = F[i];

  for (int i = 0; i < F.get_current_size(); i++)
    if (f[i] != f[i])
      throw std::logic_error(__FUNCTION__);
}

template <typename scalar_type>
void collective_sum_interface<MPI_LIBRARY>::sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& f) {
  LIN_ALG::matrix<scalar_type, LIN_ALG::CPU> F("F", f.get_current_size(), f.get_global_size());

  assert(f.get_global_size().first == F.get_global_size().first);
  assert(f.get_global_size().second == F.get_global_size().second);

  int Nr = f.get_global_size().first;
  int Nc = f.get_global_size().second;

  MPI_Allreduce(&f(0, 0), &F(0, 0), type_map_interface<MPI_LIBRARY, scalar_type>::factor() * Nr * Nc,
                type_map_interface<MPI_LIBRARY, scalar_type>::value(), MPI_SUM, grouping.get());

  for (int j = 0; j < F.get_current_size().second; j++)
    for (int i = 0; i < F.get_current_size().first; i++)
      f(i, j) = F(i, j);
}

template <typename some_type>
void collective_sum_interface<MPI_LIBRARY>::sum_and_average(some_type& obj, int size) {
  sum(obj);

  double one_over_N = 1. / (size * grouping.get_Nr_threads());

  obj *= one_over_N;
}

template <typename scalar_type, class domain>
void collective_sum_interface<MPI_LIBRARY>::average_and_compute_stddev(
    FUNC_LIB::function<scalar_type, domain>& f_mean,
    FUNC_LIB::function<scalar_type, domain>& f_stddev, size_t size) {
  scalar_type factor = 1. / (size * grouping.get_Nr_threads());

  FUNC_LIB::function<scalar_type, domain> f_sum("f-sum");
  FUNC_LIB::function<scalar_type, domain> f_diff("f-diff");

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

  f_stddev /= std::sqrt(grouping.get_Nr_threads());
}

template <typename scalar_type, class domain>
void collective_sum_interface<MPI_LIBRARY>::average_and_compute_stddev(
    FUNC_LIB::function<std::complex<scalar_type>, domain>& f_mean,
    FUNC_LIB::function<std::complex<scalar_type>, domain>& f_stddev, size_t size) {
  scalar_type factor = 1. / (size * grouping.get_Nr_threads());

  FUNC_LIB::function<std::complex<scalar_type>, domain> f_sum("f-sum");
  FUNC_LIB::function<std::complex<scalar_type>, domain> f_diff("f-diff");

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

  f_stddev /= std::sqrt(grouping.get_Nr_threads());
}

template <typename Scalar, class domain>
void collective_sum_interface<MPI_LIBRARY>::computeCovariance(
    FUNC_LIB::function<Scalar, domain>& f, FUNC_LIB::function<Scalar, domain>& f_estimated,
    FUNC_LIB::function<Scalar, dmn_variadic<domain, domain>>& cov) {
  for (int i = 0; i < f.size(); i++)
    for (int j = 0; j < f.size(); j++)
      cov(i, j) = (f(i) - f_estimated(i)) * (f(j) - f_estimated(j));

  sum_and_average(cov, 1);
}

template <typename Scalar, class domain, class cov_domain>
void collective_sum_interface<MPI_LIBRARY>::computeCovariance(
    FUNC_LIB::function<std::complex<Scalar>, domain>& f,
    FUNC_LIB::function<std::complex<Scalar>, domain>& f_estimated,
    FUNC_LIB::function<Scalar, cov_domain>& cov) {
  assert(4 * f.size() * f.size() == cov.size());

  // Treat real and imaginary parts as independent entries
  for (int i = 0; i < f.size(); i++)
    for (int j = 0; j < f.size(); j++) {

      // Real - Real
      cov(i, j) = (f(i).real() - f_estimated(i).real()) * (f(j).real() - f_estimated(j).real());

      // Imaginary - Imaginary
      cov(i + f.size(), j + f.size()) =
          (f(i).imag() - f_estimated(i).imag()) * (f(j).imag() - f_estimated(j).imag());

      // Real - Inaginary
      cov(i, j + f.size()) =
          (f(i).real() - f_estimated(i).real()) * (f(j).imag() - f_estimated(j).imag());

      // Imaginary - Real
      cov(i + f.size(), j) =
          (f(i).imag() - f_estimated(i).imag()) * (f(j).real() - f_estimated(j).imag());
    }
  sum_and_average(cov, 1);
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_INTERFACES_COLLECTIVE_SUM_INTERFACE_MPI_H
