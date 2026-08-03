// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)

#ifndef DCA_PHYS_PARAMETERS_DISORDER_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_DISORDER_PARAMETERS_HPP

#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

/** This class handles the disorder parameters
 *  The default values result in no disorder configurations
 *  I'd rather have the code know whether the section was present and
 *  Therefore whether to enable any of the disorder flow at all.
 *  Even though that is a change to the semantics of the parameters
 *  I think it may be better design than defaults that result in no
 *  disorder being used do to zero length loops or branching on
 *  num_disorder_configurations.
 */
class DisorderParameters {
public:
  DisorderParameters()
      : disorder_distribution_("binary"),
        disorder_potential_(0.0),
        disorder_density_(0.0),
        disorder_num_configurations_(0),
        // disorder_max_sites_(1),
        disorder_unique_configs_(false),
        disorder_present_(false) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_disorder_distribution() const {
    return disorder_distribution_;
  }
  double get_disorder_potential() const {
    return disorder_potential_;
  }
  void set_disorder_potential(double potential) {
    disorder_potential_ = potential;
  }
  double get_disorder_density() const {
    return disorder_density_;
  }
  int get_disorder_num_configurations() const {
    return disorder_num_configurations_;
  }
  void set_disorder_num_configurations(int num_configurations) {
    disorder_num_configurations_ = num_configurations;
  }
  // int get_disorder_max_sites() const {
  //   return disorder_max_sites_;
  // }
  bool get_disorder_unique_configs() const {
    return disorder_unique_configs_;
  }
  // Whether a "disorder" section was present in the input. Lets the disorder flow distinguish a
  // deliberate clean run (no section) from a misconfigured one (section present, zero configs).
  bool get_disorder_present() const {
    return disorder_present_;
  }

private:
  std::string disorder_distribution_;
  double disorder_potential_;
  double disorder_density_;
  int disorder_num_configurations_;
  // int disorder_max_sites_;
  bool disorder_unique_configs_;
  bool disorder_present_;
};

template <typename Concurrency>
int DisorderParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(disorder_distribution_);
  buffer_size += concurrency.get_buffer_size(disorder_potential_);
  buffer_size += concurrency.get_buffer_size(disorder_density_);
  buffer_size += concurrency.get_buffer_size(disorder_num_configurations_);
  // buffer_size += concurrency.get_buffer_size(disorder_max_sites_);
  buffer_size += concurrency.get_buffer_size(disorder_unique_configs_);
  buffer_size += concurrency.get_buffer_size(disorder_present_);

  return buffer_size;
}

template <typename Concurrency>
void DisorderParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                              int& position) const {
  concurrency.pack(buffer, buffer_size, position, disorder_distribution_);
  concurrency.pack(buffer, buffer_size, position, disorder_potential_);
  concurrency.pack(buffer, buffer_size, position, disorder_density_);
  concurrency.pack(buffer, buffer_size, position, disorder_num_configurations_);
  // concurrency.pack(buffer, buffer_size, position, disorder_max_sites_);
  concurrency.pack(buffer, buffer_size, position, disorder_unique_configs_);
  concurrency.pack(buffer, buffer_size, position, disorder_present_);
}

template <typename Concurrency>
void DisorderParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                                int& position) {
  concurrency.unpack(buffer, buffer_size, position, disorder_distribution_);
  concurrency.unpack(buffer, buffer_size, position, disorder_potential_);
  concurrency.unpack(buffer, buffer_size, position, disorder_density_);
  concurrency.unpack(buffer, buffer_size, position, disorder_num_configurations_);
  // concurrency.unpack(buffer, buffer_size, position, disorder_max_sites_);
  concurrency.unpack(buffer, buffer_size, position, disorder_unique_configs_);
  concurrency.unpack(buffer, buffer_size, position, disorder_present_);
}

template <typename ReaderOrWriter>
void DisorderParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    disorder_present_ = reader_or_writer.open_group("disorder");

    try {
      reader_or_writer.execute("distribution", disorder_distribution_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("potential", disorder_potential_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("density", disorder_density_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("num-configurations", disorder_num_configurations_);
    }
    catch (const std::exception& r_e) {
    }
    // try {
    //   reader_or_writer.execute("max-sites", disorder_max_sites_);
    // }
    // catch (const std::exception& r_e) {
    // }
    try {
      reader_or_writer.execute("unique-configs", disorder_unique_configs_);
    }
    catch (const std::exception& r_e) {
    }
    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }

  // Check values.
  if (!(disorder_distribution_ == "binary" || disorder_distribution_ == "box"))
    throw std::logic_error("Illegal value for disorder distribution (expected 'binary' or 'box').");
}

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_DISORDER_PARAMETERS_HPP
