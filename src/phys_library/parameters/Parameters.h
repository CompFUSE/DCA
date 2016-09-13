// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// Description

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMTERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMTERS_H

#include <iostream>
#include <string>
#include <vector>

#include "dca/phys/parameters/brillouin_zone_parameters.hpp"
#include "dca/phys/parameters/cpe_parameters.hpp"
#include "dca/phys/parameters/dca_parameters.hpp"
#include "dca/phys/parameters/double_counting_parameters.hpp"
#include "dca/phys/parameters/equal_time_parameters.hpp"
#include "dca/phys/parameters/ed_solver_parameters.hpp"
#include "dca/phys/parameters/filename_parameters.hpp"
#include "dca/phys/parameters/function_parameters.hpp"
#include "dca/phys/parameters/mc_solver_parameters.hpp"
#include "dca/phys/parameters/mci_parameters.hpp"
#include "dca/phys/parameters/model_parameters.hpp"
#include "dca/phys/parameters/physics_parameters.hpp"
#include "dca/phys/parameters/vertex_parameters.hpp"
#include "dca/util/print_type.hpp"
#include "comp_library/function_library/domains/special_domains/dmn_0.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/cluster/cluster_domain_family.h"
#include "phys_library/domains/cluster/cluster_domain_initializer.h"
#include "phys_library/domains/cluster/cluster_domain_symmetry_initializer.h"
#include "phys_library/domains/Quantum_domain/DCA_iteration_domain.h"
#include "phys_library/domains/Quantum_domain/numerical_error_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_real_axis.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_imag_axis.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"
#include "phys_library/domains/time_and_frequency/time_domain_left_oriented.h"
#include "phys_library/domains/time_and_frequency/vertex_time_domain.h"

// TODO: Remove when the Parameters class itself is put inside the namespace dca::phys::params.
using dca::phys::params::BrillouinZoneParameters;
using dca::phys::params::CpeParameters;
using dca::phys::params::DcaParameters;
using dca::phys::params::DoubleCountingParameters;
using dca::phys::params::EdSolverParameters;
using dca::phys::params::EqualTimeParameters;
using dca::phys::params::FilenameParameters;
using dca::phys::params::FunctionParameters;
using dca::phys::params::McSolverParameters;
using dca::phys::params::MciParameters;
using dca::phys::params::ModelParameters;
using dca::phys::params::PhysicsParameters;
using dca::phys::params::VertexParameters;

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
class Parameters : public FilenameParameters,
                   public PhysicsParameters,
                   public ModelParameters<model_t>,
                   public DcaParameters,

                   public FunctionParameters,
                   public MciParameters,

                   public McSolverParameters<solver_name>,
                   public EdSolverParameters,

                   public VertexParameters<model_t::lattice_type::DIMENSION>,
                   public EqualTimeParameters,
                   public CpeParameters,

                   public BrillouinZoneParameters,
                   public DoubleCountingParameters {
public:
  using concurrency_type = concurrency_t;
  using profiler_type = Profiler;
  using random_number_generator = rng_t;
  using model_type = model_t;
  using lattice_type = typename model_t::lattice_type;

  // Time and frequency domains
  using t = dmn_0<time_domain>;
  using w = dmn_0<frequency_domain>;
  using w_VERTEX_EXTENDED = dmn_0<DCA::vertex_frequency_domain<DCA::EXTENDED>>;
  using w_VERTEX_EXTENDED_POS = dmn_0<DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>>;

  // DCA cluster domains
  using r_DCA =
      dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, CLUSTER, REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_DCA =
      dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  // Host cluster domains
  using r_HOST =
      dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, LATTICE_SP, REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST = dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, LATTICE_SP,
                                      MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  // Host vertex cluster domains
  using r_HOST_VERTEX =
      dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, LATTICE_TP, REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST_VERTEX = dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, LATTICE_TP,
                                             MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  // LDA cluster domains
  using r_LDA = dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, LATTICE_SP,
                                     REAL_SPACE, PARALLELLEPIPEDUM>>;
  using k_LDA = dmn_0<cluster_domain<double, model_t::lattice_type::DIMENSION, LATTICE_SP,
                                     MOMENTUM_SPACE, PARALLELLEPIPEDUM>>;

  using DCA_cluster_family_type =
      cluster_domain_family<double, model_t::lattice_type::DIMENSION, CLUSTER, BRILLOUIN_ZONE>;
  using HOST_sp_cluster_family_type =
      cluster_domain_family<double, model_t::lattice_type::DIMENSION, LATTICE_SP, BRILLOUIN_ZONE>;
  using HOST_tp_cluster_family_type =
      cluster_domain_family<double, model_t::lattice_type::DIMENSION, LATTICE_TP, BRILLOUIN_ZONE>;

  constexpr static int lattice_dimension = model_t::lattice_type::DIMENSION;

#ifdef DCA_WITH_SINGLE_PRECISION_MEASUREMENTS
  typedef float MC_measurement_scalar_type;
#else
  typedef double MC_measurement_scalar_type;
#endif  // DCA_WITH_SINGLE_PRECISION_MEASUREMENTS

#ifdef DCA_WITH_REDUCED_VERTEX_FUNCTION
  typedef w_VERTEX_EXTENDED_POS G4_w1_dmn_t;
  typedef w_VERTEX_EXTENDED G4_w2_dmn_t;
#else
  typedef w_VERTEX_EXTENDED G4_w1_dmn_t;
  typedef w_VERTEX_EXTENDED G4_w2_dmn_t;
#endif  // DCA_WITH_REDUCED_VERTEX_FUNCTION

public:
  Parameters(std::string version_stamp, concurrency_type& concurrency_obj);

  template <typename Reader>
  void read(Reader& reader);

  template <typename Writer>
  void write(Writer& writer);

  template <typename Reader>
  void read_input_and_broadcast(std::string file_name);

  template <typename read_write_t>
  void read_write(read_write_t& read_write_obj);

  void update_model();
  void update_domains();

  int get_buffer_size(concurrency_type& concurrency);
  void pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);
  void unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  // get_functions
  concurrency_type& get_concurrency();

private:
  static std::string make_python_readable(std::string tmp);

private:
  std::string version_stamp;

  std::string date_str;
  std::string time_str;
  std::string compiler_str;

  concurrency_type& concurrency_obj;
};

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::Parameters(
    std::string version_stamp_str, concurrency_t& concurrency_object)
    : FilenameParameters(),
      PhysicsParameters(),
      ModelParameters<model_t>(),
      DcaParameters(model_t::DIMENSION),

      FunctionParameters(model_t::DIMENSION),
      MciParameters(),

      McSolverParameters<solver_name>(),
      VertexParameters<model_t::DIMENSION>(),
      EqualTimeParameters(),
      CpeParameters(),

      BrillouinZoneParameters(),

      version_stamp(make_python_readable(version_stamp_str)),

      date_str(__DATE__),
      time_str(__TIME__),

      compiler_str("????"),

      concurrency_obj(concurrency_object) {
#ifdef _CRAYC
  compiler_str = _RELEASE_STRING;
#endif

#ifdef __GNUC__
  compiler_str = __VERSION__;
#endif
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
template <typename Writer>
void Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::write(Writer& writer) {
  {
    writer.open_group("parameters");

    this->read_write(writer);

    writer.close_group();
  }

  {
    writer.open_group("domains");

    DCA_cluster_family_type::write(writer);
    HOST_sp_cluster_family_type::write(writer);
    HOST_tp_cluster_family_type::write(writer);

    t::parameter_type::write(writer);
    w::parameter_type::write(writer);

    DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>::write(writer);

    DCA_iteration_domain::write(writer);

    if (VertexParameters<model_t::DIMENSION>::get_vertex_measurement_type() != NONE) {
      DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN>::write(writer);
      DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN>::write(writer);
      DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN_POSITIVE>::write(writer);
      DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN_POSITIVE>::write(writer);

      DCA::vertex_frequency_domain<DCA::COMPACT>::write(writer);
      DCA::vertex_frequency_domain<DCA::EXTENDED>::write(writer);

      DCA::vertex_frequency_domain<DCA::COMPACT_POSITIVE>::write(writer);
      DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>::write(writer);
    }

    if (do_CPE()) {
      frequency_domain_real_axis::write(writer);
      frequency_domain_imag_axis::write(writer);
    }

#ifdef DCA_WITH_QMC_BIT
    numerical_error_domain::write(writer);
#endif  // DCA_WITH_QMC_BIT

    writer.close_group();
  }
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
template <typename Reader>
void Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::read_input_and_broadcast(
    std::string filename) {
  if (concurrency_obj.id() == concurrency_obj.first()) {
    Reader read_obj;
    read_obj.open_file(filename);
    this->read_write(read_obj);
    read_obj.close_file();
  }

  concurrency_obj.broadcast_object(*this);
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
void Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::update_model() {
  model_t::initialize(*this);
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
void Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::update_domains() {
  DCA_iteration_domain::initialize(*this);
  electron_band_domain::initialize(*this, model_t::BANDS, model_t::get_flavors(),
                                   model_t::get_a_vectors());

  // time && frequency-domains
  time_domain::initialize(*this);
  time_domain_left_oriented::initialize(*this);
  frequency_domain::initialize(*this);
  frequency_domain_real_axis::initialize(*this);
  frequency_domain_imag_axis::initialize(*this);

  DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN>::initialize(*this);
  DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN>::initialize(*this);
  DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN_POSITIVE>::initialize(*this);
  DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN_POSITIVE>::initialize(*this);

  DCA::vertex_frequency_domain<DCA::COMPACT>::initialize(*this);
  DCA::vertex_frequency_domain<DCA::EXTENDED>::initialize(*this);

  DCA::vertex_frequency_domain<DCA::COMPACT_POSITIVE>::initialize(*this);
  DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>::initialize(*this);

  DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>::initialize(*this);

  // DCA
  cluster_domain_initializer<r_DCA>::execute(model_t::get_r_DCA_basis(),
                                             DcaParameters::get_DCA_cluster());

  cluster_domain_symmetry_initializer<r_DCA, typename model_t::lattice_type::DCA_point_group>::execute();

  if (concurrency_obj.id() == concurrency_obj.last())
    k_DCA::parameter_type::print(std::cout);

  // host
  cluster_domain_initializer<r_HOST>::execute(
      model_t::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      FunctionParameters::get_sp_cluster());

  cluster_domain_symmetry_initializer<r_HOST, typename model_t::lattice_type::DCA_point_group>::execute();

  if (concurrency_obj.id() == concurrency_obj.last())
    k_HOST::parameter_type::print(std::cout);

  // host
  cluster_domain_initializer<r_HOST_VERTEX>::execute(
      model_t::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      FunctionParameters::get_tp_cluster());

  cluster_domain_symmetry_initializer<r_HOST_VERTEX,
                                      typename model_t::lattice_type::DCA_point_group>::execute();

  if (concurrency_obj.id() == concurrency_obj.last())
    k_HOST_VERTEX::parameter_type::print(std::cout);

  // LDA
  cluster_domain_initializer<r_LDA>::execute(
      model_t::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      FunctionParameters::get_H_k_grid_size());

  if (concurrency_obj.id() == concurrency_obj.last())
    k_LDA::parameter_type::print(std::cout);
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
template <typename read_write_t>
void Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::read_write(
    read_write_t& read_write_obj) {
  if (read_write_obj.is_writer()) {
    read_write_obj.execute("date", date_str);
    read_write_obj.execute("time", time_str);
    read_write_obj.execute("compiler", compiler_str);

    // TODO: Remove this redundant object that is only needed since the execute function of reader
    // types expects an l-value.
    std::string rng_type_str = dca::util::Type<rng_t>::print();
    read_write_obj.execute("random-number-generator", rng_type_str);
  }

  FilenameParameters::readWrite(read_write_obj);

  PhysicsParameters::readWrite(read_write_obj);
  ModelParameters<model_t>::readWrite(read_write_obj);
  DcaParameters::readWrite(read_write_obj);

  MciParameters::readWrite(read_write_obj);
  McSolverParameters<solver_name>::readWrite(read_write_obj);
  EdSolverParameters::readWrite(read_write_obj);

  FunctionParameters::readWrite(read_write_obj);
  VertexParameters<model_t::DIMENSION>::readWrite(read_write_obj);
  EqualTimeParameters::readWrite(read_write_obj);

  CpeParameters::readWrite(read_write_obj);
  DoubleCountingParameters::readWrite(read_write_obj);
  BrillouinZoneParameters::readWrite(read_write_obj);
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
int Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::get_buffer_size(
    concurrency_t& concurrency) {
  int buffer_size = 0;

  buffer_size += FilenameParameters::getBufferSize(concurrency);

  buffer_size += PhysicsParameters::getBufferSize(concurrency);
  buffer_size += ModelParameters<model_t>::getBufferSize(concurrency);
  buffer_size += DcaParameters::getBufferSize(concurrency);

  buffer_size += MciParameters::getBufferSize(concurrency);
  buffer_size += McSolverParameters<solver_name>::getBufferSize(concurrency);
  buffer_size += EdSolverParameters::getBufferSize(concurrency);

  buffer_size += FunctionParameters::getBufferSize(concurrency);
  buffer_size += VertexParameters<model_t::DIMENSION>::getBufferSize(concurrency);
  buffer_size += EqualTimeParameters::getBufferSize(concurrency);
  buffer_size += CpeParameters::getBufferSize(concurrency);

  buffer_size += BrillouinZoneParameters::getBufferSize(concurrency);
  buffer_size += DoubleCountingParameters::getBufferSize(concurrency);

  return buffer_size;
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
void Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::pack(
    concurrency_t& concurrency, int* buffer, int buffer_size, int& position) {
  FilenameParameters::pack(concurrency, buffer, buffer_size, position);

  PhysicsParameters::pack(concurrency, buffer, buffer_size, position);
  ModelParameters<model_t>::pack(concurrency, buffer, buffer_size, position);
  DcaParameters::pack(concurrency, buffer, buffer_size, position);

  MciParameters::pack(concurrency, buffer, buffer_size, position);
  McSolverParameters<solver_name>::pack(concurrency, buffer, buffer_size, position);
  EdSolverParameters::pack(concurrency, buffer, buffer_size, position);

  FunctionParameters::pack(concurrency, buffer, buffer_size, position);
  VertexParameters<model_t::DIMENSION>::pack(concurrency, buffer, buffer_size, position);
  EqualTimeParameters::pack(concurrency, buffer, buffer_size, position);
  CpeParameters::pack(concurrency, buffer, buffer_size, position);

  BrillouinZoneParameters::pack(concurrency, buffer, buffer_size, position);
  DoubleCountingParameters::pack(concurrency, buffer, buffer_size, position);
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
void Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::unpack(
    concurrency_t& concurrency, int* buffer, int buffer_size, int& position) {
  FilenameParameters::unpack(concurrency, buffer, buffer_size, position);

  PhysicsParameters::unpack(concurrency, buffer, buffer_size, position);
  ModelParameters<model_t>::unpack(concurrency, buffer, buffer_size, position);
  DcaParameters::unpack(concurrency, buffer, buffer_size, position);

  MciParameters::unpack(concurrency, buffer, buffer_size, position);
  McSolverParameters<solver_name>::unpack(concurrency, buffer, buffer_size, position);
  EdSolverParameters::unpack(concurrency, buffer, buffer_size, position);

  FunctionParameters::unpack(concurrency, buffer, buffer_size, position);
  VertexParameters<model_t::DIMENSION>::unpack(concurrency, buffer, buffer_size, position);
  EqualTimeParameters::unpack(concurrency, buffer, buffer_size, position);
  CpeParameters::unpack(concurrency, buffer, buffer_size, position);

  BrillouinZoneParameters::unpack(concurrency, buffer, buffer_size, position);
  DoubleCountingParameters::unpack(concurrency, buffer, buffer_size, position);
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
concurrency_t& Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::get_concurrency() {
  return concurrency_obj;
}

template <typename concurrency_t, typename Profiler, typename model_t, typename rng_t,
          DCA::ClusterSolverName solver_name>
std::string Parameters<concurrency_t, Profiler, model_t, rng_t, solver_name>::make_python_readable(
    std::string str) {
  {
    std::string tmp("\n\n");
    while (true) {
      int pos = str.find(tmp);
      if (pos == -1)
        break;
      str.replace(pos, tmp.size(), std::string(", "));
    }
  }

  {
    std::string tmp("\n");
    while (true) {
      int pos = str.find(tmp);
      if (pos == -1)
        break;
      str.replace(pos, tmp.size(), std::string(", "));
    }
  }

  {
    std::string tmp("\t");
    while (true) {
      int pos = str.find(tmp);
      if (pos == -1)
        break;
      str.replace(pos, tmp.size(), std::string(" "));
    }
  }

  return str;
}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMTERS_H
