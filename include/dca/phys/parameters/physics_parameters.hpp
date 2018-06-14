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
// This class reads, stores, and writes the physics parameters.

#ifndef DCA_PHYS_PARAMETERS_PHYSICS_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_PHYSICS_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class PhysicsParameters {
public:
  PhysicsParameters()
      : beta_(1.), density_(1.), chemical_potential_(0.), adjust_chemical_potential_(true) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_beta() const {
    return beta_;
  }
  double get_density() const {
    return density_;
  }
  // TODO: Need to return a reference since this function is used to set and broadcast the chemical
  //       potential. Also cannot make the function const because of this.
  double& get_chemical_potential() {
    return chemical_potential_;
  }
  bool adjust_chemical_potential() const {
    return adjust_chemical_potential_;
  }

private:
  double beta_;
  double density_;
  double chemical_potential_;
  bool adjust_chemical_potential_;
};

template <typename Concurrency>
int PhysicsParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(beta_);
  buffer_size += concurrency.get_buffer_size(density_);
  buffer_size += concurrency.get_buffer_size(chemical_potential_);
  buffer_size += concurrency.get_buffer_size(adjust_chemical_potential_);

  return buffer_size;
}

template <typename Concurrency>
void PhysicsParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                             int& position) const {
  concurrency.pack(buffer, buffer_size, position, beta_);
  concurrency.pack(buffer, buffer_size, position, density_);
  concurrency.pack(buffer, buffer_size, position, chemical_potential_);
  concurrency.pack(buffer, buffer_size, position, adjust_chemical_potential_);
}

template <typename Concurrency>
void PhysicsParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                               int& position) {
  concurrency.unpack(buffer, buffer_size, position, beta_);
  concurrency.unpack(buffer, buffer_size, position, density_);
  concurrency.unpack(buffer, buffer_size, position, chemical_potential_);
  concurrency.unpack(buffer, buffer_size, position, adjust_chemical_potential_);
}
template <typename ReaderOrWriter>
void PhysicsParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("physics");

    try {
      reader_or_writer.execute("beta", beta_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("density", density_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("chemical-potential", chemical_potential_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("adjust-chemical-potential", adjust_chemical_potential_);
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

#endif  // DCA_PHYS_PARAMETERS_PHYSICS_PARAMETERS_HPP
