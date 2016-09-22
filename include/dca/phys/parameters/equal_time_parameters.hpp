// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the equal-time parameters.

#ifndef DCA_PHYS_PARAMETERS_EQUAL_TIME_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_EQUAL_TIME_PARAMETERS_HPP

#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class EqualTimeParameters {
public:
  EqualTimeParameters() : do_equal_time_measurements_("false") {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  bool do_equal_time_measurements() const {
    return (do_equal_time_measurements_ == "true");
  }

private:
  std::string do_equal_time_measurements_;
};

template <typename Concurrency>
int EqualTimeParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;
  buffer_size += concurrency.get_buffer_size(do_equal_time_measurements_);
  return buffer_size;
}

template <typename Concurrency>
void EqualTimeParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                               int& position) const {
  concurrency.pack(buffer, buffer_size, position, do_equal_time_measurements_);
}

template <typename Concurrency>
void EqualTimeParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                                 int& position) {
  concurrency.unpack(buffer, buffer_size, position, do_equal_time_measurements_);
}

template <typename ReaderOrWriter>
void EqualTimeParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("equal-time-observables");

    try {
      reader_or_writer.execute("do-equal-time-measurements", do_equal_time_measurements_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_EQUAL_TIME_PARAMETERS_HPP
