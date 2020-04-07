// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
        sweeps_per_measurement_{1.},
        measurements_{100},
        walkers_(1),
        accumulators_(1),
        shared_walk_and_accumulation_thread_(false),
        // TODO: consider setting default do true.
        fix_meas_per_walker_(false),
        adjust_self_energy_for_double_counting_(false),
        error_computation_type_(ErrorComputationType::NONE),
        store_configuration_(false) {}

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

  const std::vector<double>& get_sweeps_per_measurement() const {
    return sweeps_per_measurement_;
  }

  const std::vector<int>& get_measurements() const {
    return measurements_;
  }
  void set_measurements(const int measurements) {
    assert(measurements >= 0);
    std::fill(measurements_.begin(), measurements_.end(), measurements);
  }

  // Maximum distance (in MC time) considered when computing the correlation between configurations.
  int get_time_correlation_window() const {
    return time_correlation_window_;
  }

  // True if the autocorrelation of G(r = 0, t = 0) is computed.
  bool compute_G_correlation() const {
    return compute_G_correlation_;
  }

  void set_time_correlation_window(int window) {
    time_correlation_window_ = window;
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

  // If true, the number of sweeps performed by each walker is fixed a priory. This avoids possible
  // bias toward faster walkers, at the expanse of load balance.
  bool fix_meas_per_walker() const {
    return fix_meas_per_walker_;
  }
  bool adjust_self_energy_for_double_counting() const {
    return adjust_self_energy_for_double_counting_;
  }
  ErrorComputationType get_error_computation_type() const {
    return error_computation_type_;
  }

  // If true, the MC configuration is stored between DCA iterations, and used to initialize the
  // walker.
  bool store_configuration() const {
    return store_configuration_;
  }
  int stamping_period() const {
    return stamping_period_;
  }

protected:
  // Resize vector arguments to have the same size as the number of iterations.
  void inline solveDcaIterationConflict(int iterations);

private:
  void generateRandomSeed() {
    std::random_device rd;
    std::uniform_int_distribution<int> dist(0, std::numeric_limits<int>::max());
    seed_ = dist(rd);
  }

  static constexpr int default_seed = 985456376;

  int seed_;
  int warm_up_sweeps_;
  std::vector<double> sweeps_per_measurement_;
  std::vector<int> measurements_;
  int measurements_final_iter_ = -1;
  int time_correlation_window_ = 0;
  bool compute_G_correlation_ = true;
  int walkers_;
  int accumulators_;
  bool shared_walk_and_accumulation_thread_;
  bool fix_meas_per_walker_;
  bool adjust_self_energy_for_double_counting_;
  ErrorComputationType error_computation_type_;
  bool store_configuration_;
  int stamping_period_ = 0;
};

template <typename Concurrency>
int MciParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(seed_);
  buffer_size += concurrency.get_buffer_size(warm_up_sweeps_);
  buffer_size += concurrency.get_buffer_size(sweeps_per_measurement_);
  buffer_size += concurrency.get_buffer_size(measurements_);
  buffer_size += concurrency.get_buffer_size(measurements_final_iter_);
  buffer_size += concurrency.get_buffer_size(time_correlation_window_);
  buffer_size += concurrency.get_buffer_size(compute_G_correlation_);
  buffer_size += concurrency.get_buffer_size(walkers_);
  buffer_size += concurrency.get_buffer_size(accumulators_);
  buffer_size += concurrency.get_buffer_size(shared_walk_and_accumulation_thread_);
  buffer_size += concurrency.get_buffer_size(fix_meas_per_walker_);
  buffer_size += concurrency.get_buffer_size(adjust_self_energy_for_double_counting_);
  buffer_size += concurrency.get_buffer_size(error_computation_type_);
  buffer_size += concurrency.get_buffer_size(store_configuration_);
  buffer_size += concurrency.get_buffer_size(stamping_period_);

  return buffer_size;
}

template <typename Concurrency>
void MciParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                         int& position) const {
  concurrency.pack(buffer, buffer_size, position, seed_);
  concurrency.pack(buffer, buffer_size, position, warm_up_sweeps_);
  concurrency.pack(buffer, buffer_size, position, sweeps_per_measurement_);
  concurrency.pack(buffer, buffer_size, position, measurements_);
  concurrency.pack(buffer, buffer_size, position, measurements_final_iter_);
  concurrency.pack(buffer, buffer_size, position, time_correlation_window_);
  concurrency.pack(buffer, buffer_size, position, compute_G_correlation_);
  concurrency.pack(buffer, buffer_size, position, walkers_);
  concurrency.pack(buffer, buffer_size, position, accumulators_);
  concurrency.pack(buffer, buffer_size, position, shared_walk_and_accumulation_thread_);
  concurrency.pack(buffer, buffer_size, position, fix_meas_per_walker_);
  concurrency.pack(buffer, buffer_size, position, adjust_self_energy_for_double_counting_);
  concurrency.pack(buffer, buffer_size, position, error_computation_type_);
  concurrency.pack(buffer, buffer_size, position, store_configuration_);
  concurrency.pack(buffer, buffer_size, position, stamping_period_);
}

template <typename Concurrency>
void MciParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                           int& position) {
  concurrency.unpack(buffer, buffer_size, position, seed_);
  concurrency.unpack(buffer, buffer_size, position, warm_up_sweeps_);
  concurrency.unpack(buffer, buffer_size, position, sweeps_per_measurement_);
  concurrency.unpack(buffer, buffer_size, position, measurements_);
  concurrency.unpack(buffer, buffer_size, position, measurements_final_iter_);
  concurrency.unpack(buffer, buffer_size, position, time_correlation_window_);
  concurrency.unpack(buffer, buffer_size, position, compute_G_correlation_);
  concurrency.unpack(buffer, buffer_size, position, walkers_);
  concurrency.unpack(buffer, buffer_size, position, accumulators_);
  concurrency.unpack(buffer, buffer_size, position, shared_walk_and_accumulation_thread_);
  concurrency.unpack(buffer, buffer_size, position, fix_meas_per_walker_);
  concurrency.unpack(buffer, buffer_size, position, adjust_self_energy_for_double_counting_);
  concurrency.unpack(buffer, buffer_size, position, error_computation_type_);
  concurrency.unpack(buffer, buffer_size, position, store_configuration_);
  concurrency.unpack(buffer, buffer_size, position, stamping_period_);
}

template <typename ReaderOrWriter>
void MciParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  auto try_to_read_write = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (std::exception&) {
    }
  };

  auto try_to_read_write_vector = [&](const std::string& name, auto& vec) {
    try {  // read as a vector.
      reader_or_writer.execute(name, vec);
    }
    catch (std::exception&) {  // read as a scalar.
      vec.resize(1);
      try_to_read_write(name, vec[0]);
    }
  };

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
        // Read the seed as an integer.
        try_to_read_write("seed", seed_);
      }
    }  // is_reader()

    else {
      // Write the seed directly.
      try_to_read_write("seed", seed_);
    }

    // Read error computation type.
    std::string error_type = toString(error_computation_type_);
    try_to_read_write("error-computation-type", error_type);
    error_computation_type_ = stringToErrorComputationType(error_type);

    try_to_read_write("warm-up-sweeps", warm_up_sweeps_);
    try_to_read_write_vector("sweeps-per-measurement", sweeps_per_measurement_);
    try_to_read_write_vector("measurements", measurements_);

    try_to_read_write("time-correlation-window", time_correlation_window_);
    try_to_read_write("compute-G-correlation", compute_G_correlation_);

    try_to_read_write("stamping-period", stamping_period_);
    try_to_read_write("store-configuration", store_configuration_);

    // Read arguments for threaded solver.
    try {
      reader_or_writer.open_group("threaded-solver");

      try_to_read_write("walkers", walkers_);
      try_to_read_write("accumulators", accumulators_);
      try_to_read_write("shared-walk-and-accumulation-thread", shared_walk_and_accumulation_thread_);
      try_to_read_write("fix-meas-per-walker", fix_meas_per_walker_);
      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    // TODO: adjust_self_energy_for_double_counting has no effect at the moment. Use default value
    // 'false'.
    // try_to_read_write("adjust-self-energy-for-double-counting", adjust_self_energy_for_double_counting_);

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }

  // Solve conflicts
  if(!time_correlation_window_)
      compute_G_correlation_ = false;
}

void MciParameters::solveDcaIterationConflict(int iterations) {
  // Solve conflicts between number of iterations and mci parameters.
  auto solve_confilct = [&](auto& param) { param.resize(iterations, param.back()); };

  solve_confilct(measurements_);
  solve_confilct(sweeps_per_measurement_);
}

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_MCI_PARAMETERS_HPP
