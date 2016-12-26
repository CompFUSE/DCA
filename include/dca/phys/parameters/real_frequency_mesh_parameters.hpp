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
// This class reads, stores, and writes the parameters that determine the real frequency mesh.

#ifndef DCA_PHYS_PARAMETERS_REAL_FREQUENCY_MESH_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_REAL_FREQUENCY_MESH_PARAMETERS_HPP

#include <stdexcept>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class RealFrequencyMeshParameters {
public:
  RealFrequencyMeshParameters()
      : min_real_frequency_(-10.),
        max_real_frequency_(10.),
        real_frequencies_(128),
        imaginary_damping_(0.01) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_min_real_frequency() const {
    return min_real_frequency_;
  }
  double get_max_real_frequency() const {
    return max_real_frequency_;
  }
  int get_real_frequencies() const {
    return real_frequencies_;
  }
  double get_imaginary_damping() const {
    return imaginary_damping_;
  }

private:
  double min_real_frequency_;
  double max_real_frequency_;
  int real_frequencies_;
  double imaginary_damping_;
};

template <typename Concurrency>
int RealFrequencyMeshParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(min_real_frequency_);
  buffer_size += concurrency.get_buffer_size(max_real_frequency_);
  buffer_size += concurrency.get_buffer_size(real_frequencies_);
  buffer_size += concurrency.get_buffer_size(imaginary_damping_);

  return buffer_size;
}

template <typename Concurrency>
void RealFrequencyMeshParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                                       int& position) const {
  concurrency.pack(buffer, buffer_size, position, min_real_frequency_);
  concurrency.pack(buffer, buffer_size, position, max_real_frequency_);
  concurrency.pack(buffer, buffer_size, position, real_frequencies_);
  concurrency.pack(buffer, buffer_size, position, imaginary_damping_);
}

template <typename Concurrency>
void RealFrequencyMeshParameters::unpack(const Concurrency& concurrency, int* buffer,
                                         int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, min_real_frequency_);
  concurrency.unpack(buffer, buffer_size, position, max_real_frequency_);
  concurrency.unpack(buffer, buffer_size, position, real_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, imaginary_damping_);
}
template <typename ReaderOrWriter>
void RealFrequencyMeshParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("real-frequency-mesh");

    try {
      reader_or_writer.execute("min", min_real_frequency_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("max", max_real_frequency_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("frequencies", real_frequencies_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("imaginary-damping", imaginary_damping_);
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

#endif  // DCA_PHYS_PARAMETERS_REAL_FREQUENCY_MESH_PARAMETERS_HPP
