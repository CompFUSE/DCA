// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class represent a subgroup of MPIProcessorGrouping.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_SERIAL_GANG_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_SERIAL_GANG_HPP

#include "dca/parallel/no_concurrency/serial_processor_grouping.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class SerialGang {
public:
  SerialGang(const SerialProcessorGrouping& /*group*/, int /*target_size*/) {}
  SerialGang(const SerialGang& other) = delete;

  int get_id() const {
    return 0;
  }
  int get_size() const {
    return 1;
  }

  auto get() const {
    return nullptr;
  }
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_SERIAL_GANG_HPP
