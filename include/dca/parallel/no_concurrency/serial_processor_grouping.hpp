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
// This class provides the (fake) processor grouping for serial execution.

#ifndef DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PROCESSOR_GROUPING_HPP
#define DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PROCESSOR_GROUPING_HPP

namespace dca {
namespace parallel {
// dca::parallel::

class SerialProcessorGrouping {
public:
  SerialProcessorGrouping() : id_(0), size_(1) {}

  int get_id() const {
    return id_;
  }
  int get_size() const {
    return size_;
  }

  int first() const {
    return 0;
  }
  int last() const {
    return size_ - 1;
  }

private:
  const int id_;
  const int size_;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PROCESSOR_GROUPING_HPP
