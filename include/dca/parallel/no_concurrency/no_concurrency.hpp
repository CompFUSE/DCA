// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class is the equivalent of MPIConcurrency for serial execution.

#ifndef DCA_PARALLEL_NO_CONCURRENCY_NO_CONCURRENCY_HPP
#define DCA_PARALLEL_NO_CONCURRENCY_NO_CONCURRENCY_HPP

#include <utility>
#include "dca/parallel/no_concurrency/serial_collective_sum.hpp"
#include "dca/parallel/no_concurrency/serial_gang.hpp"
#include "dca/parallel/no_concurrency/serial_gather.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class NoConcurrency : public SerialCollectiveSum, public SerialProcessorGrouping, public SerialGather {
public:
  using Gang = SerialGang;

  NoConcurrency(int /*argc*/, char** /*argv*/){};

  void abort() const;

  int id() const {
    return 0;
  }
  int number_of_processors() const {
    return 1;
  }

  int first() const {
    return 0;
  }
  int last() const {
    return 0;
  }

  template <typename T>
  bool broadcast(const T& /*object*/, int /*root_id*/ = 0) const {
    return true;
  }
  template <typename T>
  bool broadcast_object(const T& /*object*/, int /*root_id*/ = 0) const {
    return true;
  }

  // TODO: Add const to function parameter 'dmn'.
  template <typename Domain>
  std::pair<int, int> get_bounds(Domain& dmn) const {
    return std::make_pair(0, dmn.get_size());
  }

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);
  friend std::ostream& operator<<(std::ostream& some_ostream, const NoConcurrency& this_concurrency);

private:
  constexpr static char parallel_type_str_[] = "NoConcurrency";
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_NO_CONCURRENCY_HPP
