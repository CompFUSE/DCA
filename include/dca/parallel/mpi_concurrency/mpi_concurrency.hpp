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
// This class provides a simple interface to the standard Message Passing Interface (MPI).

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_CONCURRENCY_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_CONCURRENCY_HPP

#include <cassert>
#include <iostream>
#include <utility>
#include <stdexcept>

#include <mpi.h>

#include "dca/parallel/mpi_concurrency/mpi_collective_max.hpp"
#include "dca/parallel/mpi_concurrency/mpi_collective_min.hpp"
#include "dca/parallel/mpi_concurrency/mpi_collective_sum.hpp"
#include "dca/parallel/mpi_concurrency/mpi_packing.hpp"
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/util/get_bounds.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPIConcurrency : public MPIPacking,
                       public MPICollectiveMax,
                       public MPICollectiveMin,
                       public MPICollectiveSum {
public:
  MPIConcurrency(int argc, char** argv)
      : MPIPacking(grouping_),
        MPICollectiveMax(grouping_),
        MPICollectiveMin(grouping_),
        MPICollectiveSum(grouping_) {
    // TODO: If there is only one grouping, the MPI_Init call can be moved to the constructor of the
    //       MPIProcessorGrouping class.
    MPI_Init(&argc, &argv);
    grouping_.set();
  }

  ~MPIConcurrency() {
    MPI_Finalize();
  }

  int id() const {
    assert(grouping_.get_id() > -1);
    return grouping_.get_id();
  }
  int number_of_processors() const {
    assert(grouping_.get_Nr_threads() > -1);
    return grouping_.get_Nr_threads();
  }
  int first() const {
    return grouping_.first();
  }
  int last() const {
    return grouping_.last();
  }

  template <typename object_type>
  bool broadcast(object_type& object, int root_id = 0) const;

  template <typename object_type>
  bool broadcast_object(object_type& object, int root_id = 0) const;

  // TODO: Add const to function parameter 'dmn'.
  template <typename domain_type>
  std::pair<int, int> get_bounds(domain_type& dmn) const {
    return util::getBounds(id(), number_of_processors(), dmn);
  }

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);
  friend std::ostream& operator << (std::ostream&, const MPIConcurrency&);  
private:
  constexpr static char concurrency_type_str[] = "MPI Concurrency";
  MPIProcessorGrouping grouping_;
};

template <typename object_type>
bool MPIConcurrency::broadcast(object_type& object, int root_id) const {
  assert(root_id > -1 and root_id < number_of_processors());

  int position = 0;

  if (id() == root_id) {
    int bufferSize = this->get_buffer_size(object);

    MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, grouping_.get());

    char* buffer = new char[bufferSize];

    this->pack(buffer, bufferSize, position, object);

    MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, grouping_.get());

    delete[] buffer;
  }
  else {
    int bufferSize(0);

    MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, grouping_.get());

    char* buffer = new char[bufferSize];

    // receive packed message
    MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, grouping_.get());

    this->unpack(buffer, bufferSize, position, object);

    delete[] buffer;
  }

  return true;
}

template <typename object_type>
bool MPIConcurrency::broadcast_object(object_type& object, int root_id) const {
  assert(root_id > -1 and root_id < number_of_processors());

  int buffer_size = 0;

  if (id() == root_id) {
    buffer_size = object.get_buffer_size(*this);

    MPI_Bcast(&buffer_size, 1, MPI_INT, root_id, grouping_.get());

    char* buffer = new char[buffer_size];

    int off_set = 0;
    object.pack(*this, buffer, buffer_size, off_set);

    MPI_Bcast(buffer, buffer_size, MPI_PACKED, root_id, grouping_.get());

    delete[] buffer;
  }
  else {
    MPI_Bcast(&buffer_size, 1, MPI_INT, root_id, grouping_.get());

    char* buffer = new char[buffer_size];

    MPI_Bcast(buffer, buffer_size, MPI_PACKED, root_id, grouping_.get());

    int off_set = 0;
    object.unpack(*this, buffer, buffer_size, off_set);

    delete[] buffer;
  }

  return true;
}

template <typename ReaderOrWriter>
void MPIConcurrency::readWrite(ReaderOrWriter& reader_or_writer) {
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

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_CONCURRENCY_HPP
