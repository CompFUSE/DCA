// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the double-counting parameters.

#ifndef DCA_PHYS_PARAMETERS_DOUBLE_COUNTING_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_DOUBLE_COUNTING_PARAMETERS_HPP

#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class DoubleCountingParameters {
public:
  DoubleCountingParameters() : double_counting_method_("none"), double_counting_correction_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_double_counting_method() const {
    return double_counting_method_;
  }
  double get_double_counting_correction() const {
    return double_counting_correction_;
  }

private:
  std::string double_counting_method_;
  double double_counting_correction_;
};

template <typename Concurrency>
int DoubleCountingParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(double_counting_method_);
  buffer_size += concurrency.get_buffer_size(double_counting_correction_);

  return buffer_size;
}

template <typename Concurrency>
void DoubleCountingParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                                    int& position) const {
  concurrency.pack(buffer, buffer_size, position, double_counting_method_);
  concurrency.pack(buffer, buffer_size, position, double_counting_correction_);
}

template <typename Concurrency>
void DoubleCountingParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                                      int& position) {
  concurrency.unpack(buffer, buffer_size, position, double_counting_method_);
  concurrency.unpack(buffer, buffer_size, position, double_counting_correction_);
}

template <typename ReaderOrWriter>
void DoubleCountingParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("double-counting");

    try {
      reader_or_writer.execute("method", double_counting_method_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("correction", double_counting_correction_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }

  // Check values.
  if (!(double_counting_method_ == "none" ||
        double_counting_method_ == "constant-correction-without-U-correction" ||
        double_counting_method_ == "constant-correction-with-U-correction"))
    throw std::logic_error("Illegal value for double-counting method.");
}

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_DOUBLE_COUNTING_PARAMETERS_HPP
