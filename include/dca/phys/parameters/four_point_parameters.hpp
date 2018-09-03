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
// This class reads, stores, and writes parameters for measurements of four-point functions.
//
// TODO: Move the computation of the 'exact' momentum transfer and its index outside of this class
//       since the parameters classes should only read and write and not do any computation.

#ifndef DCA_PHYS_PARAMETERS_FOUR_POINT_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_FOUR_POINT_PARAMETERS_HPP

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
      : four_point_type_(NONE),
        four_point_momentum_transfer_input_(lattice_dimension, 0.),
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

  FourPointType get_four_point_type() const {
    return static_cast<FourPointType>(four_point_type_);
  }

  void set_four_point_type(FourPointType type) {
    four_point_type_ = type;
  }

  const std::vector<double>& get_four_point_momentum_transfer_input() const {
    return four_point_momentum_transfer_input_;
  }

  // Returns the index of the bosonic exchange frequency if compute_all_transfers() == false, or the
  // index of the maximum transfer if compute_all_transfers() == true.
  int get_four_point_frequency_transfer() const {
    return four_point_frequency_transfer_;
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

  // Returns true if all possible momentum and frequency exchanges are computed, ignoring the values
  // of 'get_four_point_momentum_transfer' and 'get_four_point_momentum_transfer_index'.
  bool compute_all_transfers() const {
    return compute_all_transfers_;
  }

private:
  // There is no utility to communicate enumerations over mpi, so four_point_type_ is stored
  // as an int rather than a FourPointType.
  int four_point_type_;
  std::vector<double> four_point_momentum_transfer_input_;
  int four_point_frequency_transfer_;
  bool compute_all_transfers_;
};

template <int lattice_dimension>
template <typename Concurrency>
int FourPointParameters<lattice_dimension>::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(four_point_type_);
  buffer_size += concurrency.get_buffer_size(four_point_momentum_transfer_input_);
  buffer_size += concurrency.get_buffer_size(four_point_frequency_transfer_);
  buffer_size += concurrency.get_buffer_size(compute_all_transfers_);

  return buffer_size;
}

template <int lattice_dimension>
template <typename Concurrency>
void FourPointParameters<lattice_dimension>::pack(const Concurrency& concurrency, char* buffer,
                                                  int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, four_point_type_);
  concurrency.pack(buffer, buffer_size, position, four_point_momentum_transfer_input_);
  concurrency.pack(buffer, buffer_size, position, four_point_frequency_transfer_);
  concurrency.pack(buffer, buffer_size, position, compute_all_transfers_);
}

template <int lattice_dimension>
template <typename Concurrency>
void FourPointParameters<lattice_dimension>::unpack(const Concurrency& concurrency, char* buffer,
                                                    int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, four_point_type_);
  concurrency.unpack(buffer, buffer_size, position, four_point_momentum_transfer_input_);
  concurrency.unpack(buffer, buffer_size, position, four_point_frequency_transfer_);
  concurrency.unpack(buffer, buffer_size, position, compute_all_transfers_);
}

template <int lattice_dimension>
template <typename ReaderOrWriter>
void FourPointParameters<lattice_dimension>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("four-point");

    std::string four_point_name = toString(static_cast<FourPointType>(four_point_type_));
    try {
      reader_or_writer.execute("type", four_point_name);
      four_point_type_ = stringToFourPointType(four_point_name);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("momentum-transfer", four_point_momentum_transfer_input_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("frequency-transfer", four_point_frequency_transfer_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("compute-all-transfers", compute_all_transfers_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }

  if (compute_all_transfers_ && four_point_frequency_transfer_ < 0)
    throw(std::logic_error(
        "When compute-all-transfers is set, a greater than 0 frequency-transfer must be chosen."));
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_FOUR_POINT_PARAMETERS_HPP
