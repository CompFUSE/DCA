// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the Monte Carlo Integration (MCI) parameters.

#ifndef DCA_PHYS_PARAMETERS_MCI_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_MCI_PARAMETERS_HPP

#include <cassert>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#include "dca/phys/error_computation_type.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class MciParameters {
public:
  MciParameters()
      : seed_(default_seed),
        warm_up_sweeps_(20),
        sweeps_per_measurement_(1),
        measurements_(100),
        walkers_(1),
        accumulators_(1),
        shared_walk_and_accumulation_thread_(false),
        adjust_self_energy_for_double_counting_(false),
        error_computation_type_(0) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  int get_seed() const {
    return seed_;
  }
  int get_warm_up_sweeps() const {
    return warm_up_sweeps_;
  }
  int get_sweeps_per_measurement() const {
    return sweeps_per_measurement_;
  }
  int get_measurements() const {
    return measurements_;
  }
  void set_measurements(const int measurements) {
    assert(measurements >= 0);
    measurements_ = measurements;
  }
  int get_walkers() const {
    return walkers_;
  }
  int get_accumulators() const {
    return accumulators_;
  }
  bool shared_walk_and_accumulation_thread() const {
    return shared_walk_and_accumulation_thread_;
  }
  bool adjust_self_energy_for_double_counting() const {
    return adjust_self_energy_for_double_counting_;
  }
  ErrorComputationType get_error_computation_type() const {
    return ErrorComputationType(error_computation_type_);
  }

private:
  void generateRandomSeed() {
    std::random_device rd;
    std::uniform_int_distribution<int> dist(0, std::numeric_limits<int>::max());
    seed_ = dist(rd);
  }

  static constexpr int default_seed = 985456376;

  int seed_;
  int warm_up_sweeps_;
  int sweeps_per_measurement_;
  int measurements_;
  int walkers_;
  int accumulators_;
  bool shared_walk_and_accumulation_thread_;
  bool adjust_self_energy_for_double_counting_;
  int error_computation_type_;
};

template <typename Concurrency>
int MciParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(seed_);
  buffer_size += concurrency.get_buffer_size(warm_up_sweeps_);
  buffer_size += concurrency.get_buffer_size(sweeps_per_measurement_);
  buffer_size += concurrency.get_buffer_size(measurements_);
  buffer_size += concurrency.get_buffer_size(walkers_);
  buffer_size += concurrency.get_buffer_size(accumulators_);
  buffer_size += concurrency.get_buffer_size(shared_walk_and_accumulation_thread_);
  buffer_size += concurrency.get_buffer_size(adjust_self_energy_for_double_counting_);
  buffer_size += concurrency.get_buffer_size(error_computation_type_);

  return buffer_size;
}

template <typename Concurrency>
void MciParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                         int& position) const {
  concurrency.pack(buffer, buffer_size, position, seed_);
  concurrency.pack(buffer, buffer_size, position, warm_up_sweeps_);
  concurrency.pack(buffer, buffer_size, position, sweeps_per_measurement_);
  concurrency.pack(buffer, buffer_size, position, measurements_);
  concurrency.pack(buffer, buffer_size, position, walkers_);
  concurrency.pack(buffer, buffer_size, position, accumulators_);
  concurrency.pack(buffer, buffer_size, position, shared_walk_and_accumulation_thread_);
  concurrency.pack(buffer, buffer_size, position, adjust_self_energy_for_double_counting_);
  concurrency.pack(buffer, buffer_size, position, error_computation_type_);
}

template <typename Concurrency>
void MciParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                           int& position) {
  concurrency.unpack(buffer, buffer_size, position, seed_);
  concurrency.unpack(buffer, buffer_size, position, warm_up_sweeps_);
  concurrency.unpack(buffer, buffer_size, position, sweeps_per_measurement_);
  concurrency.unpack(buffer, buffer_size, position, measurements_);
  concurrency.unpack(buffer, buffer_size, position, walkers_);
  concurrency.unpack(buffer, buffer_size, position, accumulators_);
  concurrency.unpack(buffer, buffer_size, position, shared_walk_and_accumulation_thread_);
  concurrency.unpack(buffer, buffer_size, position, adjust_self_energy_for_double_counting_);
  concurrency.unpack(buffer, buffer_size, position, error_computation_type_);
}

template <typename ReaderOrWriter>
void MciParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("Monte-Carlo-integration");

    if (reader_or_writer.is_reader()) {
      // The input file can contain an integral seed or the seeding option "random".
      try {
        // Try to read a seeding option.
        std::string seed_string;
        reader_or_writer.execute("seed", seed_string);
        if (seed_string == "random")
          generateRandomSeed();
        else {
          std::cerr << "Warning: Invalid seeding option. Using default seed = " << default_seed
                    << "." << std::endl;
          seed_ = default_seed;
        }
      }
      catch (const std::exception& r_e) {
        try {
          // Read the seed as an integer.
          reader_or_writer.execute("seed", seed_);
        }

        catch (const std::exception& r_e2) {
        }
      }
    }

    else {
      // Write the seed.
      try {
        reader_or_writer.execute("seed", seed_);
      }
      catch (const std::exception& r_e) {
      }
    }

    try {
      reader_or_writer.execute("warm-up-sweeps", warm_up_sweeps_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("sweeps-per-measurement", sweeps_per_measurement_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("measurements", measurements_);
    }
    catch (const std::exception& r_e) {
    }

    // Read error computation type.
    std::string error_type = toString(static_cast<ErrorComputationType>(error_computation_type_));
    try {
      reader_or_writer.execute("error-computation-type", error_type);
      error_computation_type_ = static_cast<int>(readErrorComputationType(error_type));
    }
    catch (const std::exception& r_e) {
    }

    // Read arguments for threaded solver.
    try {
      reader_or_writer.open_group("threaded-solver");
      try {
        reader_or_writer.execute("walkers", walkers_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("accumulators", accumulators_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("shared-walk-and-accumulation-thread",
                                 shared_walk_and_accumulation_thread_);
      }
      catch (const std::exception& r_e) {
      }
      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    // TODO: adjust_self_energy_for_double_counting has no effect at the moment. Use default value
    // 'false'.
    // try {
    //   reader_or_writer.execute("adjust-self-energy-for-double-counting",
    //                            adjust_self_energy_for_double_counting_);
    // }
    // catch (const std::exception& r_e) {
    // }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_MCI_PARAMETERS_HPP
