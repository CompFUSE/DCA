// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class is the equivalent of Pthreading for serial execution.

#ifndef DCA_PARALLEL_NO_THREADING_HPP
#define DCA_PARALLEL_NO_THREADING_HPP

#include "dca/parallel/util/threading_data.hpp"
#include <iostream>
#include <stdexcept>

namespace dca {
namespace parallel {
// dca::parallel::

class NoThreading {
public:
  void execute(int num_threads, void* (*start_routine)(void*), void* arg) {
    for (int id = 0; id < num_threads; id++) {
      data_.id = id;
      data_.num_threads = num_threads;
      data_.arg = arg;
      start_routine(static_cast<void*>(&data_));
    }
  }

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);
  friend std::ostream& operator << (std::ostream&, const NoThreading&);
private:
  constexpr static char concurrency_type_str[] = "No Threading Concurrency";
  ThreadingData data_;
};

template <typename ReaderOrWriter>
void NoThreading::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("concurrency");
    try {
      reader_or_writer.execute("type", concurrency_type_str);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("number_of_threads", data_.num_threads);
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

#endif  // DCA_PARALLEL_NO_THREADING_HPP
