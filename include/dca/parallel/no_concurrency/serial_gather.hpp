// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class provides a no-op interface for gather functions without MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_SERIAL_GATHER_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_SERIAL_GATHER_HPP

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

class SerialGather {
public:
  SerialGather() = default;

  template <class Scalar, class DmnIn, class DmnOut, class Gang>
  void gather(const func::function<Scalar, DmnIn>& f_in, func::function<Scalar, DmnOut>& f_out,
              const Gang& /*gang*/) const {
    // TODO: move.
    if (f_in.size() != f_out.size())
      throw(std::logic_error("Size mismatch."));
    std::copy_n(f_in.values(), f_in.size(), f_out.values());
  }
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_SERIAL_GATHER_HPP
