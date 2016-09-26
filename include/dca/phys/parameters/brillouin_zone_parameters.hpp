// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the Brillouin zone parameters.
//
// TODO: This set of 'parameters' is usually not read from the input file but computed from
//       compute_band_structure. This is why the getter function return non-const references.

#ifndef DCA_PHYS_PARAMETERS_BRILLOUIN_ZONE_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_BRILLOUIN_ZONE_PARAMETERS_HPP

#include <stdexcept>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class BrillouinZoneParameters {
public:
  BrillouinZoneParameters()
      : coordinate_type_("absolute"), coordinate_names_(0), Brillouin_zone_vectors_(0) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  std::string& get_coordinate_type() {
    return coordinate_type_;
  }
  std::vector<std::string>& get_coordinate_names() {
    return coordinate_names_;
  }
  std::vector<std::vector<double>>& get_Brillouin_zone_vectors() {
    return Brillouin_zone_vectors_;
  }

private:
  std::string coordinate_type_;
  std::vector<std::string> coordinate_names_;
  std::vector<std::vector<double>> Brillouin_zone_vectors_;
};

template <typename Concurrency>
int BrillouinZoneParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(coordinate_type_);
  buffer_size += concurrency.get_buffer_size(coordinate_names_);
  buffer_size += concurrency.get_buffer_size(Brillouin_zone_vectors_);

  return buffer_size;
}

template <typename Concurrency>
void BrillouinZoneParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                                   int& position) const {
  concurrency.pack(buffer, buffer_size, position, coordinate_type_);
  concurrency.pack(buffer, buffer_size, position, coordinate_names_);
  concurrency.pack(buffer, buffer_size, position, Brillouin_zone_vectors_);
}

template <typename Concurrency>
void BrillouinZoneParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                                     int& position) {
  concurrency.unpack(buffer, buffer_size, position, coordinate_type_);
  concurrency.unpack(buffer, buffer_size, position, coordinate_names_);
  concurrency.unpack(buffer, buffer_size, position, Brillouin_zone_vectors_);
}

template <typename ReaderOrWriter>
void BrillouinZoneParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("band-structure-cut");

    try {
      reader_or_writer.execute("coordinate-type", coordinate_type_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("Brillouin-zone-names", coordinate_names_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("Brillouin-zone-vectors", Brillouin_zone_vectors_);
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

#endif  // DCA_PHYS_PARAMETERS_BRILLOUIN_ZONE_PARAMETERS_HPP
