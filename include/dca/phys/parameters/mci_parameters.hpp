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
//
// TODO: Const correctness.

#ifndef DCA_PHYS_PARAMETERS_MCI_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_MCI_PARAMETERS_HPP

#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class MciParameters {
public:
  MciParameters()
      : Sigma_file_("zero"),
        warm_up_sweeps_(20),
        number_of_sweeps_per_measurement_(1.),
        number_of_measurements_(100),
        do_adaptive_double_counting_("false"),
        seed_(default_seed),
        nr_walkers_(1),
        nr_accumulators_(1),
        additional_steps_(1),
        nr_HTS_threads_(1) {}

  template <typename Concurrency>
  int getBufferSize(Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_Sigma_file() const {
    return Sigma_file_;
  };
  int get_warm_up_sweeps() const {
    return warm_up_sweeps_;
  }
  double get_number_of_sweeps_per_measurement() const {
    return number_of_sweeps_per_measurement_;
  }
  int get_number_of_measurements() const {
    return number_of_measurements_;
  }
  bool adjust_self_energy_for_double_counting() const {
    return (do_adaptive_double_counting_ == "true");
  }
  int get_seed() const {
    return seed_;
  }
  int get_nr_walkers() const {
    return nr_walkers_;
  }
  int get_nr_accumulators() const {
    return nr_accumulators_;
  }
  int get_additional_steps() const {
    return additional_steps_;
  }
  int get_nr_HTS_threads() const {
    return nr_HTS_threads_;
  }

private:
  void generateRandomSeed() {
    std::random_device rd;
    std::uniform_int_distribution<int> dist(0, std::numeric_limits<int>::max());
    seed_ = dist(rd);
  }

  static constexpr int default_seed = 985456376;

  std::string Sigma_file_;
  int warm_up_sweeps_;
  double number_of_sweeps_per_measurement_;
  int number_of_measurements_;
  std::string do_adaptive_double_counting_;
  int seed_;
  int nr_walkers_;
  int nr_accumulators_;
  int additional_steps_;
  int nr_HTS_threads_;
};

template <typename Concurrency>
int MciParameters::getBufferSize(Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(Sigma_file_);
  buffer_size += concurrency.get_buffer_size(warm_up_sweeps_);
  buffer_size += concurrency.get_buffer_size(number_of_sweeps_per_measurement_);
  buffer_size += concurrency.get_buffer_size(number_of_measurements_);
  buffer_size += concurrency.get_buffer_size(do_adaptive_double_counting_);
  buffer_size += concurrency.get_buffer_size(seed_);
  buffer_size += concurrency.get_buffer_size(nr_walkers_);
  buffer_size += concurrency.get_buffer_size(nr_accumulators_);
  buffer_size += concurrency.get_buffer_size(additional_steps_);
  buffer_size += concurrency.get_buffer_size(nr_HTS_threads_);

  return buffer_size;
}

template <typename Concurrency>
void MciParameters::pack(Concurrency& concurrency, int* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, Sigma_file_);
  concurrency.pack(buffer, buffer_size, position, warm_up_sweeps_);
  concurrency.pack(buffer, buffer_size, position, number_of_sweeps_per_measurement_);
  concurrency.pack(buffer, buffer_size, position, number_of_measurements_);
  concurrency.pack(buffer, buffer_size, position, do_adaptive_double_counting_);
  concurrency.pack(buffer, buffer_size, position, seed_);
  concurrency.pack(buffer, buffer_size, position, nr_walkers_);
  concurrency.pack(buffer, buffer_size, position, nr_accumulators_);
  concurrency.pack(buffer, buffer_size, position, additional_steps_);
  concurrency.pack(buffer, buffer_size, position, nr_HTS_threads_);
}

template <typename Concurrency>
void MciParameters::unpack(Concurrency& concurrency, int* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, Sigma_file_);
  concurrency.unpack(buffer, buffer_size, position, warm_up_sweeps_);
  concurrency.unpack(buffer, buffer_size, position, number_of_sweeps_per_measurement_);
  concurrency.unpack(buffer, buffer_size, position, number_of_measurements_);
  concurrency.unpack(buffer, buffer_size, position, do_adaptive_double_counting_);
  concurrency.unpack(buffer, buffer_size, position, seed_);
  concurrency.unpack(buffer, buffer_size, position, nr_walkers_);
  concurrency.unpack(buffer, buffer_size, position, nr_accumulators_);
  concurrency.unpack(buffer, buffer_size, position, additional_steps_);
  concurrency.unpack(buffer, buffer_size, position, nr_HTS_threads_);
}

template <typename ReaderOrWriter>
void MciParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("Monte-Carlo-Integration");

    try {
      reader_or_writer.execute("Sigma-file", Sigma_file_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("warm-up-sweeps", warm_up_sweeps_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("sweeps-per-measurement", number_of_sweeps_per_measurement_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("measurements", number_of_measurements_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("adaptive-double-counting", do_adaptive_double_counting_);
    }
    catch (const std::exception& r_e) {
    }

    if (reader_or_writer.is_reader()) {
      // The input file can contain an integral seed or the seeding option "random".
      try {
        // Try to read a seeding option.
        std::string seed_string;
        reader_or_writer.execute("RNG-seed", seed_string);
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
          reader_or_writer.execute("RNG-seed", seed_);
        }

        catch (const std::exception& r_e2) {
        }
      }
    }

    else {
      // Write the RNG seed.
      try {
        reader_or_writer.execute("RNG-seed", seed_);
      }
      catch (const std::exception& r_e) {
      }
    }

    {
      reader_or_writer.open_group("MC-posix-parameters");

      try {
        reader_or_writer.execute("nr-walkers", nr_walkers_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("nr-accumulators", nr_accumulators_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("additional-steps", additional_steps_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("HTS-threads", nr_HTS_threads_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo MCI parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_MCI_PARAMETERS_HPP
