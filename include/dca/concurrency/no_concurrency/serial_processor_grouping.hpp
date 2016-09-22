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
// This class provides the (fake) processor grouping for serial execution.

#ifndef DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PROCESSOR_GROUPING_HPP
#define DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PROCESSOR_GROUPING_HPP

namespace dca {
namespace concurrency {
// dca::concurrency::

class SerialProcessorGrouping {
public:
  SerialProcessorGrouping() : id_(0), nr_threads_(1) {}

  int get_id() const {
    return id_;
  }
  int get_Nr_threads() const {
    return nr_threads_;
  }

  int first() const {
    return 0;
  }
  int last() const {
    return nr_threads_ - 1;
  }

private:
  const int id_;
  const int nr_threads_;
};

}  // concurrency
}  // dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PROCESSOR_GROUPING_HPP
