// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class reads, stores, and writes parameters for measurements of four-point functions.
//
// TODO: Move the computation of the 'exact' momentum transfer and its index outside of this class
//       since the parameters classes should only read and write and not do any computation.

#ifndef DCA_PHYS_PARAMETERS_FOUR_POINT_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_FOUR_POINT_PARAMETERS_HPP

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

template <int lattice_dimension>
class FourPointParameters {
public:
  using DCA_k_cluster_type =
      domains::cluster_domain<double, lattice_dimension, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  FourPointParameters()
      : four_point_momentum_transfer_input_(lattice_dimension, 0.),
        four_point_frequency_transfer_(0),
        compute_all_transfers_(false) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);
  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const auto& get_four_point_channels() const {
    return four_point_channels_;
  }

  /** somehow ranks other than the head rank get four_point_channels == FourPointType::PARTICLE_HOLE_NONE
   *  when the head node has no four point channels. This works around this bug.  \todo fix this.
   */
  bool isAccumulatingG4() const {
    if (four_point_channels_.size() > 0 &&
        four_point_channels_[0] != FourPointType::PARTICLE_HOLE_NONE)
      return true;
    else
      return false;
  }

  /** problematic set function.
   *  elsewhere in the code the size of FourPointTypes in four_point_channels_ is used to
   *  determinine if tp_accumulation should be done. This creates a special case where
   * PARTICLE_HOLE_NONE should never be added.
   */
  void set_four_point_channel(FourPointType type) {
    if (type != FourPointType::PARTICLE_HOLE_NONE)
      four_point_channels_ = std::vector<FourPointType>{type};
  }
  void set_four_point_channels(const std::vector<FourPointType>& channels) {
    four_point_channels_ = channels;
  }

  const std::vector<double>& get_four_point_momentum_transfer_input() const {
    return four_point_momentum_transfer_input_;
  }

  // Returns the index of the bosonic exchange frequency. If compute_all_transfers() is true, all
  // non-negative frequency transfers up to and included this index are computed. Otherwise only the
  // transfer relative to this index is computed.
  int get_four_point_frequency_transfer() const {
    return four_point_frequency_transfer_;
  }

  // Returns the vector of indexes of bosonic exchange frequencies
  std::vector<int> get_four_point_frequency_transfers() const {
    std::vector<int> transfers;
    if (compute_all_transfers_)
      for (int i = 0; i < four_point_frequency_transfer_; ++i)
        transfers.push_back(i);
    else
      transfers.push_back(four_point_frequency_transfer_);
    return transfers;
  }

  // Returns the 'exact' momentum transfer (q-vector), i.e. the DCA momentum space cluster vector
  // whose distance (L2 norm) to the input momentum transfer is minimal.
  // It assumes that the input q-vectors' distance to the next DCA momentum space cluster vector is
  // smaller than 10^-3.
  const std::vector<double>& get_four_point_momentum_transfer() const {
    static const std::vector<double> q_vec(domains::cluster_operations::find_closest_cluster_vector(
        four_point_momentum_transfer_input_, DCA_k_cluster_type::get_elements(),
        DCA_k_cluster_type::get_super_basis_vectors(), 1.e-3));
    return q_vec;
  }
  // Returns the index of the 'exact' momentum transfer (q-vector) with respect to the elements
  // vector of the DCA momentum space cluster.
  const int& get_four_point_momentum_transfer_index() const {
    static const int q_ind(domains::cluster_operations::index(get_four_point_momentum_transfer(),
                                                              DCA_k_cluster_type::get_elements(),
                                                              DCA_k_cluster_type::SHAPE));
    return q_ind;
  }

  std::size_t get_bin_measurements_every() const {
    return bin_measurements_every_;
  }

  // Returns true if all possible momentum and frequency exchanges are computed, ignoring the values
  // of 'get_four_point_momentum_transfer' and 'get_four_point_momentum_transfer_index'.
  bool compute_all_transfers() const {
    return compute_all_transfers_;
  }

private:
  std::vector<FourPointType> four_point_channels_;
  std::vector<double> four_point_momentum_transfer_input_;
  int four_point_frequency_transfer_ = 0;
  bool compute_all_transfers_ = false;
  std::size_t bin_measurements_every_ = 0;
};

template <int lattice_dimension>
template <typename Concurrency>
int FourPointParameters<lattice_dimension>::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(four_point_channels_);
  buffer_size += concurrency.get_buffer_size(four_point_momentum_transfer_input_);
  buffer_size += concurrency.get_buffer_size(four_point_frequency_transfer_);
  buffer_size += concurrency.get_buffer_size(compute_all_transfers_);
  buffer_size += concurrency.get_buffer_size(bin_measurements_every_);

  return buffer_size;
}

template <int lattice_dimension>
template <typename Concurrency>
void FourPointParameters<lattice_dimension>::pack(const Concurrency& concurrency, char* buffer,
                                                  int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, four_point_channels_);
  concurrency.pack(buffer, buffer_size, position, four_point_momentum_transfer_input_);
  concurrency.pack(buffer, buffer_size, position, four_point_frequency_transfer_);
  concurrency.pack(buffer, buffer_size, position, compute_all_transfers_);
  concurrency.pack(buffer, buffer_size, position, bin_measurements_every_);
}

template <int lattice_dimension>
template <typename Concurrency>
void FourPointParameters<lattice_dimension>::unpack(const Concurrency& concurrency, char* buffer,
                                                    int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, four_point_channels_);
  concurrency.unpack(buffer, buffer_size, position, four_point_momentum_transfer_input_);
  concurrency.unpack(buffer, buffer_size, position, four_point_frequency_transfer_);
  concurrency.unpack(buffer, buffer_size, position, compute_all_transfers_);
  concurrency.unpack(buffer, buffer_size, position, bin_measurements_every_);
}

template <int lattice_dimension>
template <typename ReaderOrWriter>
void FourPointParameters<lattice_dimension>::readWrite(ReaderOrWriter& reader_or_writer) {
  if (reader_or_writer.open_group("four-point")) {
    std::vector<std::string> channel_names;
    const std::string channel_par_name = "channels";
    if (ReaderOrWriter::is_reader) {
      // Support legacy input files specifying a single channel name.
      reader_or_writer.execute(channel_par_name, channel_names);
      for (auto name : channel_names)
        four_point_channels_.push_back(stringToFourPointType(name));
      std::string legacy_channel_name;
      reader_or_writer.execute("type", legacy_channel_name);
      if (!legacy_channel_name.empty())
        four_point_channels_.push_back(stringToFourPointType(legacy_channel_name));

      // Remove duplicates
      std::sort(four_point_channels_.begin(), four_point_channels_.end());
      four_point_channels_.erase(std::unique(four_point_channels_.begin(), four_point_channels_.end()),
                                 four_point_channels_.end());
    }

    else {  // is writer.
      for (auto channel : four_point_channels_)
        channel_names.push_back(toString(channel));
      reader_or_writer.execute(channel_par_name, channel_names);
    }

    reader_or_writer.execute("momentum-transfer", four_point_momentum_transfer_input_);
    reader_or_writer.execute("frequency-transfer", four_point_frequency_transfer_);
    reader_or_writer.execute("compute-all-transfers", compute_all_transfers_);
    reader_or_writer.execute("bin-measurements-every", bin_measurements_every_);

    if (compute_all_transfers_ && four_point_frequency_transfer_ < 0)
      throw(
          std::logic_error("When compute-all-transfers is set, a greater than 0 frequency-transfer "
                           "must be chosen."));
  }
  reader_or_writer.close_group();
}

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_FOUR_POINT_PARAMETERS_HPP
