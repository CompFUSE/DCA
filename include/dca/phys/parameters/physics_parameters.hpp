// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the physics parameters.

#ifndef DCA_PHYS_PARAMETERS_PHYSICS_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_PHYSICS_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class PhysicsParameters {
public:
  PhysicsParameters()
      : beta_(1.), adjust_chemical_potential_("false"), density_(1.), chemical_potential_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_beta() const {
    return beta_;
  }
  bool adjust_chemical_potential() const {
    return (adjust_chemical_potential_ == "true");
  }
  double get_density() const {
    return density_;
  }
  // TODO: Need to return a reference since this function is used to set and broadcast the chemical
  //       potential. Also cannot make the function const because of this.
  double& get_chemical_potential() {
    return chemical_potential_;
  }

private:
  double beta_;
  std::string adjust_chemical_potential_;
  double density_;
  double chemical_potential_;
};

template <typename Concurrency>
int PhysicsParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(beta_);
  buffer_size += concurrency.get_buffer_size(adjust_chemical_potential_);
  buffer_size += concurrency.get_buffer_size(density_);
  buffer_size += concurrency.get_buffer_size(chemical_potential_);

  return buffer_size;
}

template <typename Concurrency>
void PhysicsParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                             int& position) const {
  concurrency.pack(buffer, buffer_size, position, beta_);
  concurrency.pack(buffer, buffer_size, position, adjust_chemical_potential_);
  concurrency.pack(buffer, buffer_size, position, density_);
  concurrency.pack(buffer, buffer_size, position, chemical_potential_);
}

template <typename Concurrency>
void PhysicsParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                               int& position) {
  concurrency.unpack(buffer, buffer_size, position, beta_);
  concurrency.unpack(buffer, buffer_size, position, adjust_chemical_potential_);
  concurrency.unpack(buffer, buffer_size, position, density_);
  concurrency.unpack(buffer, buffer_size, position, chemical_potential_);
}
template <typename ReaderOrWriter>
void PhysicsParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("physics-parameters");

    try {
      reader_or_writer.execute("beta", beta_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("adjust-chemical-potential", adjust_chemical_potential_);
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

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo physics parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_PHYSICS_PARAMETERS_HPP
