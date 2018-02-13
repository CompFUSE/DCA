// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class is the equivalent of MPIConcurrency for serial execution.

#ifndef DCA_PARALLEL_NO_CONCURRENCY_NO_CONCURRENCY_HPP
#define DCA_PARALLEL_NO_CONCURRENCY_NO_CONCURRENCY_HPP

#include <utility>
#include "dca/parallel/no_concurrency/serial_collective_sum.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class NoConcurrency : public SerialCollectiveSum {
public:
  NoConcurrency(int /*argc*/, char** /*argv*/){};

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
  friend std::ostream& operator << (std::ostream &, const NoConcurrency&);
private:
  constexpr static char concurrency_type_str[] = "No Concurrency";
};

template <typename ReaderOrWriter>
void NoConcurrency::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("concurrency");
    try {
      reader_or_writer.execute("type", concurrency_type_str);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("number_of_processors", number_of_processors());
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("grouping.first", first());
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("grouping.last", last());
    }
    catch (const std::exception& r_e) {
    }    
    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }  
}
  
}  // parallel
}  // dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_NO_CONCURRENCY_HPP
