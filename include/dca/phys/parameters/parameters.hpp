// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class manages the simulation parameters.

#ifndef DCA_PHYS_PARAMETERS_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_PARAMETERS_HPP

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

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
class Parameters : public FilenameParameters,
                   public PhysicsParameters,
                   public ModelParameters<Model>,
                   public DcaParameters,
                   public FunctionParameters,
                   public MciParameters,
                   public McSolverParameters<solver_name>,
                   public EdSolverParameters,
                   public VertexParameters<Model::lattice_type::DIMENSION>,
                   public EqualTimeParameters,
                   public CpeParameters,
                   public BrillouinZoneParameters,
                   public DoubleCountingParameters {
public:
  using concurrency_type = Concurrency;
  using profiler_type = Profiler;
  using random_number_generator = RandomNumberGenerator;
  using model_type = Model;
  using lattice_type = typename Model::lattice_type;

  // Time and frequency domains
  using t = dmn_0<time_domain>;
  using w = dmn_0<frequency_domain>;
  using w_VERTEX_EXTENDED = dmn_0<DCA::vertex_frequency_domain<DCA::EXTENDED>>;
  using w_VERTEX_EXTENDED_POS = dmn_0<DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>>;

  // DCA cluster domains
  using r_DCA =
      dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, CLUSTER, REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_DCA =
      dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  // Host cluster domains
  using r_HOST =
      dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, LATTICE_SP, REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST = dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, LATTICE_SP,
                                      MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  // Host vertex cluster domains
  using r_HOST_VERTEX =
      dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, LATTICE_TP, REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST_VERTEX = dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, LATTICE_TP,
                                             MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  // LDA cluster domains
  using r_LDA =
      dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, LATTICE_SP, REAL_SPACE, PARALLELLEPIPEDUM>>;
  using k_LDA = dmn_0<cluster_domain<double, Model::lattice_type::DIMENSION, LATTICE_SP,
                                     MOMENTUM_SPACE, PARALLELLEPIPEDUM>>;

  using DCA_cluster_family_type =
      cluster_domain_family<double, Model::lattice_type::DIMENSION, CLUSTER, BRILLOUIN_ZONE>;
  using HOST_sp_cluster_family_type =
      cluster_domain_family<double, Model::lattice_type::DIMENSION, LATTICE_SP, BRILLOUIN_ZONE>;
  using HOST_tp_cluster_family_type =
      cluster_domain_family<double, Model::lattice_type::DIMENSION, LATTICE_TP, BRILLOUIN_ZONE>;

  constexpr static int lattice_dimension = Model::lattice_type::DIMENSION;

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

  Parameters(const std::string& version_stamp, concurrency_type& concurrency);

  template <typename Writer>
  void write(Writer& writer);
  template <typename Reader>
  void read_input_and_broadcast(const std::string& file_name);

  void update_model();
  void update_domains();

  int get_buffer_size(concurrency_type& concurrency) const;
  void pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position) const;
  void unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  concurrency_type& get_concurrency() {
    return concurrency_;
  }

private:
  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  std::string make_python_readable(std::string tmp);

  std::string version_stamp_;

  std::string date_;
  std::string time_;
  std::string compiler_;

  concurrency_type& concurrency_;
};

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::Parameters(
    const std::string& version_stamp, concurrency_type& concurrency)
    : FilenameParameters(),
      PhysicsParameters(),
      ModelParameters<Model>(),
      DcaParameters(Model::DIMENSION),
      FunctionParameters(Model::DIMENSION),
      MciParameters(),
      McSolverParameters<solver_name>(),
      VertexParameters<Model::DIMENSION>(),
      EqualTimeParameters(),
      CpeParameters(),
      BrillouinZoneParameters(),

      version_stamp_(make_python_readable(version_stamp)),

      date_(__DATE__),
      time_(__TIME__),
      compiler_("????"),

      concurrency_(concurrency) {
#ifdef _CRAYC
  compiler_ = _RELEASE_STRING;
#endif

#ifdef __GNUC__
  compiler_ = __VERSION__;
#endif
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
template <typename Writer>
void Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::write(Writer& writer) {
  writer.open_group("parameters");
  this->readWrite(writer);
  writer.close_group();

  writer.open_group("domains");

  DCA_cluster_family_type::write(writer);
  HOST_sp_cluster_family_type::write(writer);
  HOST_tp_cluster_family_type::write(writer);

  t::parameter_type::write(writer);
  w::parameter_type::write(writer);

  DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>::write(writer);

  DCA_iteration_domain::write(writer);

  if (VertexParameters<Model::DIMENSION>::get_vertex_measurement_type() != NONE) {
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

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
template <typename Reader>
void Parameters<Concurrency, Profiler, Model, RandomNumberGenerator,
                solver_name>::read_input_and_broadcast(const std::string& filename) {
  if (concurrency_.id() == concurrency_.first()) {
    Reader read_obj;
    read_obj.open_file(filename);
    this->readWrite(read_obj);
    read_obj.close_file();
  }

  concurrency_.broadcast_object(*this);
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
void Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::update_model() {
  Model::initialize(*this);
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
void Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::update_domains() {
  DCA_iteration_domain::initialize(*this);
  electron_band_domain::initialize(*this, Model::BANDS, Model::get_flavors(), Model::get_a_vectors());

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
  cluster_domain_initializer<r_DCA>::execute(Model::get_r_DCA_basis(),
                                             DcaParameters::get_DCA_cluster());

  cluster_domain_symmetry_initializer<r_DCA, typename Model::lattice_type::DCA_point_group>::execute();

  if (concurrency_.id() == concurrency_.last())
    k_DCA::parameter_type::print(std::cout);

  // host
  cluster_domain_initializer<r_HOST>::execute(
      Model::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      FunctionParameters::get_sp_cluster());

  cluster_domain_symmetry_initializer<r_HOST, typename Model::lattice_type::DCA_point_group>::execute();

  if (concurrency_.id() == concurrency_.last())
    k_HOST::parameter_type::print(std::cout);

  // host
  cluster_domain_initializer<r_HOST_VERTEX>::execute(
      Model::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      FunctionParameters::get_tp_cluster());

  cluster_domain_symmetry_initializer<r_HOST_VERTEX,
                                      typename Model::lattice_type::DCA_point_group>::execute();

  if (concurrency_.id() == concurrency_.last())
    k_HOST_VERTEX::parameter_type::print(std::cout);

  // LDA
  cluster_domain_initializer<r_LDA>::execute(
      Model::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      FunctionParameters::get_H_k_grid_size());

  if (concurrency_.id() == concurrency_.last())
    k_LDA::parameter_type::print(std::cout);
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
int Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::get_buffer_size(
    Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += FilenameParameters::getBufferSize(concurrency);
  buffer_size += PhysicsParameters::getBufferSize(concurrency);
  buffer_size += ModelParameters<Model>::getBufferSize(concurrency);
  buffer_size += DcaParameters::getBufferSize(concurrency);
  buffer_size += MciParameters::getBufferSize(concurrency);
  buffer_size += McSolverParameters<solver_name>::getBufferSize(concurrency);
  buffer_size += EdSolverParameters::getBufferSize(concurrency);
  buffer_size += FunctionParameters::getBufferSize(concurrency);
  buffer_size += VertexParameters<Model::DIMENSION>::getBufferSize(concurrency);
  buffer_size += EqualTimeParameters::getBufferSize(concurrency);
  buffer_size += CpeParameters::getBufferSize(concurrency);
  buffer_size += BrillouinZoneParameters::getBufferSize(concurrency);
  buffer_size += DoubleCountingParameters::getBufferSize(concurrency);

  return buffer_size;
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
void Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::pack(
    Concurrency& concurrency, int* buffer, int buffer_size, int& position) const {
  FilenameParameters::pack(concurrency, buffer, buffer_size, position);
  PhysicsParameters::pack(concurrency, buffer, buffer_size, position);
  ModelParameters<Model>::pack(concurrency, buffer, buffer_size, position);
  DcaParameters::pack(concurrency, buffer, buffer_size, position);
  MciParameters::pack(concurrency, buffer, buffer_size, position);
  McSolverParameters<solver_name>::pack(concurrency, buffer, buffer_size, position);
  EdSolverParameters::pack(concurrency, buffer, buffer_size, position);
  FunctionParameters::pack(concurrency, buffer, buffer_size, position);
  VertexParameters<Model::DIMENSION>::pack(concurrency, buffer, buffer_size, position);
  EqualTimeParameters::pack(concurrency, buffer, buffer_size, position);
  CpeParameters::pack(concurrency, buffer, buffer_size, position);
  BrillouinZoneParameters::pack(concurrency, buffer, buffer_size, position);
  DoubleCountingParameters::pack(concurrency, buffer, buffer_size, position);
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
void Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::unpack(
    Concurrency& concurrency, int* buffer, int buffer_size, int& position) {
  FilenameParameters::unpack(concurrency, buffer, buffer_size, position);
  PhysicsParameters::unpack(concurrency, buffer, buffer_size, position);
  ModelParameters<Model>::unpack(concurrency, buffer, buffer_size, position);
  DcaParameters::unpack(concurrency, buffer, buffer_size, position);
  MciParameters::unpack(concurrency, buffer, buffer_size, position);
  McSolverParameters<solver_name>::unpack(concurrency, buffer, buffer_size, position);
  EdSolverParameters::unpack(concurrency, buffer, buffer_size, position);
  FunctionParameters::unpack(concurrency, buffer, buffer_size, position);
  VertexParameters<Model::DIMENSION>::unpack(concurrency, buffer, buffer_size, position);
  EqualTimeParameters::unpack(concurrency, buffer, buffer_size, position);
  CpeParameters::unpack(concurrency, buffer, buffer_size, position);
  BrillouinZoneParameters::unpack(concurrency, buffer, buffer_size, position);
  DoubleCountingParameters::unpack(concurrency, buffer, buffer_size, position);
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
template <typename ReaderOrWriter>
void Parameters<Concurrency, Profiler, Model, RandomNumberGenerator, solver_name>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  if (reader_or_writer.is_writer()) {
    reader_or_writer.execute("date", date_);
    reader_or_writer.execute("time", time_);
    reader_or_writer.execute("compiler", compiler_);

    // TODO: Remove this redundant object that is only needed since the execute function of reader
    //       types expects an l-value.
    std::string RandomNumberGeneratorype_str = dca::util::Type<RandomNumberGenerator>::print();
    reader_or_writer.execute("random-number-generator", RandomNumberGeneratorype_str);
  }

  FilenameParameters::readWrite(reader_or_writer);
  PhysicsParameters::readWrite(reader_or_writer);
  ModelParameters<Model>::readWrite(reader_or_writer);
  DcaParameters::readWrite(reader_or_writer);
  MciParameters::readWrite(reader_or_writer);
  McSolverParameters<solver_name>::readWrite(reader_or_writer);
  EdSolverParameters::readWrite(reader_or_writer);
  FunctionParameters::readWrite(reader_or_writer);
  VertexParameters<Model::DIMENSION>::readWrite(reader_or_writer);
  EqualTimeParameters::readWrite(reader_or_writer);
  CpeParameters::readWrite(reader_or_writer);
  DoubleCountingParameters::readWrite(reader_or_writer);
  BrillouinZoneParameters::readWrite(reader_or_writer);
}

template <typename Concurrency, typename Profiler, typename Model, typename RandomNumberGenerator,
          DCA::ClusterSolverName solver_name>
std::string Parameters<Concurrency, Profiler, Model, RandomNumberGenerator,
                       solver_name>::make_python_readable(std::string str) {
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

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_PARAMETERS_HPP
