// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains all functions on the real frequency axis.

#ifndef DCA_PHYS_DCA_DATA_DCA_DATA_REAL_FREQ_HPP
#define DCA_PHYS_DCA_DATA_DCA_DATA_REAL_FREQ_HPP

#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <class parameters_type>
class DcaDataRealFreq {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  DcaDataRealFreq(parameters_type& parameters_ref);

  void read(std::string filename);
  template <typename Reader>
  void read(Reader& reader);

  void write(std::string filename);
  template <typename Writer>
  void write(Writer& writer);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

public:
  func::function<double, w_REAL> A_w;
  func::function<double, w_REAL> A_w_stddev;

  func::function<double, func::dmn_variadic<b, s, w_REAL>> A_nu_w;
  func::function<double, func::dmn_variadic<b, s, w_REAL>> A_nu_w_stddev;

  func::function<double, w_REAL> A0_w;
  func::function<double, func::dmn_variadic<b, s, w_REAL>> A0_nu_w;

  func::function<double, w_REAL> E_w;
  func::function<double, w_REAL> E0_w;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w_REAL>> Sigma;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w_REAL>> Sigma_stddev;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w_REAL>> G0_k_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w_REAL>> G0_r_w;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w_REAL>> G_k_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w_REAL>> G_k_w_stddev;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w_REAL>> G_r_w;
};

template <class parameters_type>
DcaDataRealFreq<parameters_type>::DcaDataRealFreq(parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      A_w("spectral-density"),
      A_w_stddev("spectral-density-stddev"),

      A_nu_w("spectral-density-per-orbital"),
      A_nu_w_stddev("spectral-density-per-orbital-stddev"),

      A0_w("free-spectral-density"),
      A0_nu_w("free-spectral-density-per-orbital"),

      E_w("E_w"),
      E0_w("E0_w"),

      Sigma("self-energy-real-axis"),
      Sigma_stddev("self-energy-real-axis-stddev"),

      G0_k_w("cluster-greens-function-G0_k_w-real-axis"),
      G0_r_w("cluster-greens-function-G0_r_w-real-axis"),

      G_k_w("cluster-greens-function-G_k_w-real-axis"),
      G_k_w_stddev("cluster-greens-function-G_k_w-real-axis-stddev"),

      G_r_w("cluster-greens-function-G_r_w-real-axis") {}

template <class parameters_type>
void DcaDataRealFreq<parameters_type>::read(std::string filename) {
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\t starts reading \n\n";

    const std::string& output_format = parameters.get_output_format();

    if (output_format == "JSON") {
      dca::io::JSONReader reader;
      reader.open_file(filename);
      this->read(reader);
      reader.close_file();
    }

    else if (output_format == "HDF5") {
      dca::io::HDF5Reader reader;
      reader.open_file(filename);
      this->read(reader);
      reader.close_file();
    }

    else
      throw std::logic_error(__FUNCTION__);
  }

  concurrency.broadcast_object(A_w);
  concurrency.broadcast_object(A0_w);

  concurrency.broadcast_object(E_w);
  concurrency.broadcast_object(E0_w);

  concurrency.broadcast_object(Sigma);

  concurrency.broadcast_object(G_k_w);
  concurrency.broadcast_object(G0_k_w);

  concurrency.broadcast_object(G_r_w);
  concurrency.broadcast_object(G0_r_w);
}

template <class parameters_type>
template <typename Reader>
void DcaDataRealFreq<parameters_type>::read(Reader& reader) {
  reader.open_group("spectral-functions");

  reader.execute(A_w);
  reader.execute(A0_w);

  reader.execute(E_w);
  reader.execute(E0_w);

  reader.execute(Sigma);

  reader.execute(G0_k_w);
  reader.execute(G_k_w);

  reader.execute(G_r_w);
  reader.execute(G0_r_w);

  reader.close_group();
}

template <class parameters_type>
void DcaDataRealFreq<parameters_type>::write(std::string filename) {
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\t starts writing \n\n";

    const std::string& output_format = parameters.get_output_format();

    if (output_format == "JSON") {
      dca::io::JSONWriter writer;
      writer.open_file(filename);
      this->write(writer);
      writer.close_file();
    }

    else if (output_format == "HDF5") {
      dca::io::HDF5Writer writer;
      writer.open_file(filename);
      this->write(writer);
      writer.close_file();
    }

    else
      throw std::logic_error(__FUNCTION__);
  }
}

template <class parameters_type>
template <typename Writer>
void DcaDataRealFreq<parameters_type>::write(Writer& writer) {
  writer.open_group("spectral-functions");

  writer.execute(A_w);
  writer.execute(A_w_stddev);

  writer.execute(A_nu_w);
  writer.execute(A_nu_w_stddev);

  writer.execute(A0_w);
  writer.execute(A0_nu_w);

  writer.execute(E_w);
  writer.execute(E0_w);

  writer.execute(Sigma);
  writer.execute(Sigma_stddev);

  writer.execute(G_k_w);
  writer.execute(G_k_w_stddev);

  writer.execute(G0_k_w);

  writer.execute(G_r_w);
  writer.execute(G0_r_w);

  writer.close_group();
}

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_DATA_DCA_DATA_REAL_FREQ_HPP
