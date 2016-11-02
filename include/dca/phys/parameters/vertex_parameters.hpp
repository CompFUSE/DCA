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
// This class reads, stores, and writes the vertex parameters.
//
// TODO: - Move the computation of the 'exact' q-channel vector and its index outside of this class
//         since the parameters classes should only read and write and not do any computation.
//       - Fix typo in json key "diagonolize-folded-Gamma-chi_0".

#ifndef DCA_PHYS_PARAMETERS_VERTEX_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_VERTEX_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/vertex_measurement_type.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

template <int lattice_dimension>
class VertexParameters {
public:
  using DCA_k_cluster_type =
      domains::cluster_domain<double, lattice_dimension, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  VertexParameters();

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  VertexMeasurementType get_vertex_measurement_type() const;

  const std::vector<double>& get_q_channel_vec_input() const {
    return q_channel_vec_input_;
  }
  // This function returns the 'exact' q-channel vector, i.e. the k-cluster vector whose distance
  // (L2 norm) to the input q-channel vector is minimal.
  // It assumes that the input q-channel vector's distance to the next k-cluster vector is smaller
  // than 10^-3.
  const std::vector<double>& get_q_channel_vec() const {
    static const std::vector<double> q_channel_vec(
        domains::cluster_operations::find_closest_cluster_vector(
            q_channel_vec_input_, DCA_k_cluster_type::get_elements(),
            DCA_k_cluster_type::get_super_basis_vectors(), 1.e-3));
    return q_channel_vec;
  }
  // Returns the index of the 'exact' q-channel vector.
  const int& get_q_channel_ind() const {
    static const int q_channel_ind(domains::cluster_operations::index(
        get_q_channel_vec(), DCA_k_cluster_type::get_elements(), DCA_k_cluster_type::SHAPE));
    return q_channel_ind;
  }

  int get_w_channel() const {
    return w_channel_;
  }

  double get_singular_value_cut_off() const {
    return singular_value_cut_off_;
  }
  int get_singular_value_index_cut_off() const {
    return singular_value_index_cut_off_;
  }

  bool do_diagonalization_on_folded_Gamma_chi_0() const {
    return (diagonalize_folded_Gamma_chi_0_ == "true");
  }
  double get_BSE_cut_off_radius() const {
    return BSE_cut_off_radius_;
  }

  bool do_deconvolution_of_Gamma() const {
    return (do_deconvolution_of_Gamma_ == "yes");
  }
  bool do_symmetrization_of_Gamma() const {
    return (do_symmetrization_of_Gamma_ == "yes");
  }

  bool compute_chi_0() const {
    return (compute_chi_0_ == "yes");
  }
  bool compute_chi() const {
    return (compute_chi_ == "yes");
  }

  bool compute_eigenvalues() const {
    return (compute_eigenvalues_ == "yes");
  }
  bool compute_P_q_cluster() const {
    return (compute_P_q_cluster_ == "yes");
  }
  bool compute_P_q_lattice() const {
    return (compute_P_q_lattice_ == "yes");
  }

private:
  std::string vertex_measurement_type_;

  std::vector<double> q_channel_vec_input_;
  int w_channel_;

  double singular_value_cut_off_;
  int singular_value_index_cut_off_;

  std::string diagonalize_folded_Gamma_chi_0_;
  double BSE_cut_off_radius_;

  std::string do_deconvolution_of_Gamma_;
  std::string do_symmetrization_of_Gamma_;

  std::string compute_chi_0_;
  std::string compute_chi_;

  std::string compute_eigenvalues_;
  std::string compute_P_q_cluster_;
  std::string compute_P_q_lattice_;
};

template <int lattice_dimension>
VertexParameters<lattice_dimension>::VertexParameters()
    : vertex_measurement_type_("NONE"),

      q_channel_vec_input_(lattice_dimension, 0),
      w_channel_(0),

      singular_value_cut_off_(0.5),
      singular_value_index_cut_off_(128),

      diagonalize_folded_Gamma_chi_0_("false"),
      BSE_cut_off_radius_(10.),

      do_deconvolution_of_Gamma_("yes"),
      do_symmetrization_of_Gamma_("yes"),

      compute_chi_0_("no"),
      compute_chi_("no"),

      compute_eigenvalues_("yes"),
      compute_P_q_cluster_("no"),
      compute_P_q_lattice_("no") {}

template <int lattice_dimension>
template <typename Concurrency>
int VertexParameters<lattice_dimension>::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(vertex_measurement_type_);
  buffer_size += concurrency.get_buffer_size(q_channel_vec_input_);
  buffer_size += concurrency.get_buffer_size(w_channel_);
  buffer_size += concurrency.get_buffer_size(singular_value_cut_off_);
  buffer_size += concurrency.get_buffer_size(singular_value_index_cut_off_);
  buffer_size += concurrency.get_buffer_size(diagonalize_folded_Gamma_chi_0_);
  buffer_size += concurrency.get_buffer_size(BSE_cut_off_radius_);
  buffer_size += concurrency.get_buffer_size(do_deconvolution_of_Gamma_);
  buffer_size += concurrency.get_buffer_size(do_symmetrization_of_Gamma_);
  buffer_size += concurrency.get_buffer_size(compute_chi_0_);
  buffer_size += concurrency.get_buffer_size(compute_chi_);
  buffer_size += concurrency.get_buffer_size(compute_eigenvalues_);
  buffer_size += concurrency.get_buffer_size(compute_P_q_cluster_);
  buffer_size += concurrency.get_buffer_size(compute_P_q_lattice_);

  return buffer_size;
}

template <int lattice_dimension>
template <typename Concurrency>
void VertexParameters<lattice_dimension>::pack(const Concurrency& concurrency, int* buffer,
                                               int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, vertex_measurement_type_);
  concurrency.pack(buffer, buffer_size, position, q_channel_vec_input_);
  concurrency.pack(buffer, buffer_size, position, w_channel_);
  concurrency.pack(buffer, buffer_size, position, singular_value_cut_off_);
  concurrency.pack(buffer, buffer_size, position, singular_value_index_cut_off_);
  concurrency.pack(buffer, buffer_size, position, diagonalize_folded_Gamma_chi_0_);
  concurrency.pack(buffer, buffer_size, position, BSE_cut_off_radius_);
  concurrency.pack(buffer, buffer_size, position, do_deconvolution_of_Gamma_);
  concurrency.pack(buffer, buffer_size, position, do_symmetrization_of_Gamma_);
  concurrency.pack(buffer, buffer_size, position, compute_chi_0_);
  concurrency.pack(buffer, buffer_size, position, compute_chi_);
  concurrency.pack(buffer, buffer_size, position, compute_eigenvalues_);
  concurrency.pack(buffer, buffer_size, position, compute_P_q_cluster_);
  concurrency.pack(buffer, buffer_size, position, compute_P_q_lattice_);
}

template <int lattice_dimension>
template <typename Concurrency>
void VertexParameters<lattice_dimension>::unpack(const Concurrency& concurrency, int* buffer,
                                                 int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, vertex_measurement_type_);
  concurrency.unpack(buffer, buffer_size, position, q_channel_vec_input_);
  concurrency.unpack(buffer, buffer_size, position, w_channel_);
  concurrency.unpack(buffer, buffer_size, position, singular_value_cut_off_);
  concurrency.unpack(buffer, buffer_size, position, singular_value_index_cut_off_);
  concurrency.unpack(buffer, buffer_size, position, diagonalize_folded_Gamma_chi_0_);
  concurrency.unpack(buffer, buffer_size, position, BSE_cut_off_radius_);
  concurrency.unpack(buffer, buffer_size, position, do_deconvolution_of_Gamma_);
  concurrency.unpack(buffer, buffer_size, position, do_symmetrization_of_Gamma_);
  concurrency.unpack(buffer, buffer_size, position, compute_chi_0_);
  concurrency.unpack(buffer, buffer_size, position, compute_chi_);
  concurrency.unpack(buffer, buffer_size, position, compute_eigenvalues_);
  concurrency.unpack(buffer, buffer_size, position, compute_P_q_cluster_);
  concurrency.unpack(buffer, buffer_size, position, compute_P_q_lattice_);
}

template <int lattice_dimension>
template <typename ReaderOrWriter>
void VertexParameters<lattice_dimension>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("vertex-channel");

    try {
      reader_or_writer.execute("vertex-measurement-type", vertex_measurement_type_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("q-channel", q_channel_vec_input_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("w-channel", w_channel_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("deconvolute-Gamma", do_deconvolution_of_Gamma_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("symmetrize-Gamma", do_symmetrization_of_Gamma_);
    }
    catch (const std::exception& r_e) {
    }

    {
      reader_or_writer.open_group("lattice-mapping");

      try {
        reader_or_writer.execute("singular-value-sigma-cut-off", singular_value_cut_off_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("singular-value-index-cut-off", singular_value_index_cut_off_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    {
      reader_or_writer.open_group("lattice-solver");

      try {
        reader_or_writer.execute("diagonolize-folded-Gamma-chi_0", diagonalize_folded_Gamma_chi_0_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("BSE-cut-off-radius", BSE_cut_off_radius_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    {
      reader_or_writer.open_group("options");

      try {
        reader_or_writer.execute("compute-chi", compute_chi_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("compute-chi_0", compute_chi_0_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("compute-eigenvalues", compute_eigenvalues_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("compute-Pq-cluster", compute_P_q_cluster_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("compute-Pq-lattice", compute_P_q_lattice_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo vertex parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template <int lattice_dimension>
VertexMeasurementType VertexParameters<lattice_dimension>::get_vertex_measurement_type() const {
  if (vertex_measurement_type_ == "PARTICLE_HOLE_TRANSVERSE")
    return PARTICLE_HOLE_TRANSVERSE;

  else if (vertex_measurement_type_ == "PARTICLE_HOLE_MAGNETIC")
    return PARTICLE_HOLE_MAGNETIC;

  else if (vertex_measurement_type_ == "PARTICLE_HOLE_CHARGE")
    return PARTICLE_HOLE_CHARGE;

  else if (vertex_measurement_type_ == "PARTICLE_PARTICLE_SUPERCONDUCTING")
    return PARTICLE_PARTICLE_SUPERCONDUCTING;

  else if (vertex_measurement_type_ == "NONE")
    return NONE;

  else {
    std::cout << vertex_measurement_type_ << " is not of the type: "
              << "PARTICLE_HOLE_TRANSVERSE"
              << ", "
              << "PARTICLE_HOLE_MAGNETIC"
              << ", "
              << "PARTICLE_HOLE_CHARGE"
              << ", "
              << "PARTICLE_PARTICLE_SUPERCONDUCTING"
              << ", "
              << "NONE"
              << ", " << std::endl;

    throw std::logic_error(__FUNCTION__);
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_VERTEX_PARAMETERS_HPP
