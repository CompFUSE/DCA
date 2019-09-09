// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class provides an interface to gather functions with MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_GATHER_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_GATHER_HPP

#include <vector>

#include <mpi.h>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/mpi_concurrency/mpi_gang.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPIGather {
public:
  MPIGather() = default;

  // Gather the function 'f_in' on all processes in 'gang' and copy the result into 'f_out',
  // discarding eventual padding.
  // Precondition: gang.get_size() * f_in.size() >= f_out.size()
  template <class Scalar, class DmnIn, class DmnOut, class Gang>
  void gather(const func::function<Scalar, DmnIn>& f_in, func::function<Scalar, DmnOut>& f_out,
              const Gang& gang) const;

private:
  template <class T, class Gang>
  void gather(const T* in, T* out, int local_n, const Gang& gang) const;
};

template <class Scalar, class DmnIn, class DmnOut, class Gang>
void MPIGather::gather(const func::function<Scalar, DmnIn>& f_in,
                       func::function<Scalar, DmnOut>& f_out, const Gang& gang) const {
  std::vector<Scalar> gathered(f_in.size() * gang.get_size());
  gather(f_in.values(), gathered.data(), f_in.size(), gang);

  if (f_out.size() > gathered.size())
    throw(std::logic_error("Output function is too large."));

  // TODO: move.
  std::copy_n(gathered.data(), f_out.size(), f_out.values());
}

template <class T, class Gang>
void MPIGather::gather(const T* in, T* out, int local_n, const Gang& gang) const {
  MPI_Allgather(in, local_n, MPITypeMap<T>::value(), out, local_n, MPITypeMap<T>::value(),
                gang.get());
}

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_GATHER_HPP
