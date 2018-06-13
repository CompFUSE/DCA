// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class manages the simulation parameters.
//
// TODO: Const correctness.

#ifndef DCA_PHYS_PARAMETERS_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_PARAMETERS_HPP

#include <iostream>
#include <string>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/phys/parameters/analysis_parameters.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/parameters/dca_parameters.hpp"
#include "dca/phys/parameters/domains_parameters.hpp"
#include "dca/phys/parameters/double_counting_parameters.hpp"
#include "dca/phys/parameters/ed_solver_parameters.hpp"
#include "dca/phys/parameters/four_point_parameters.hpp"
#include "dca/phys/parameters/mc_solver_parameters.hpp"
#include "dca/phys/parameters/mci_parameters.hpp"
#include "dca/phys/parameters/model_parameters.hpp"
#include "dca/phys/parameters/output_parameters.hpp"
#include "dca/phys/parameters/physics_parameters.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_family.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/cluster_domain_symmetry_initializer.hpp"
#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"
#include "dca/phys/domains/quantum/dca_iteration_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/numerical_error_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_time_domain.hpp"
#include "dca/util/print_type.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
class Parameters : public AnalysisParameters,
                   public DcaParameters,
                   public DomainsParameters,
                   public DoubleCountingParameters,
                   public EdSolverParameters,
                   public FourPointParameters<Model::lattice_type::DIMENSION>,
                   public McSolverParameters<solver_name>,
                   public MciParameters,
                   public ModelParameters<Model>,
                   public OutputParameters,
                   public PhysicsParameters {
public:
  using concurrency_type = Concurrency;
  using ThreadingType = Threading;
  using profiler_type = Profiler;
  using random_number_generator = RandomNumberGenerator;
  using model_type = Model;
  using lattice_type = typename Model::lattice_type;

  // Time and frequency domains
  using TDmn = func::dmn_0<domains::time_domain>;
  using WDmn = func::dmn_0<domains::frequency_domain>;
  using WTpExtDmn = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;
  using WTpExtPosDmn = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED_POSITIVE>>;

  constexpr static int lattice_dimension = Model::lattice_type::DIMENSION;
  using CDA = ClusterDomainAliases<lattice_dimension>;
  // DCA cluster domains
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  // Host cluster domains
  using RSpHostDmn = typename CDA::RSpHostDmn;
  using KSpHostDmn = typename CDA::KSpHostDmn;

  // Host vertex cluster domains
  using RTpHostDmn = typename CDA::RTpHostDmn;
  using KTpHostDmn = typename CDA::KTpHostDmn;

  using DcaClusterFamily = typename CDA::DcaClusterFamily;
  using HostSpClusterFamily = typename CDA::HostSpClusterFamily;
  using HostTpClusterFamily = typename CDA::HostTpClusterFamily;

  constexpr static int bands = Model::lattice_type::BANDS;

#ifdef DCA_WITH_SINGLE_PRECISION_MEASUREMENTS
  typedef float MC_measurement_scalar_type;
#else
  typedef double MC_measurement_scalar_type;
#endif  // DCA_WITH_SINGLE_PRECISION_MEASUREMENTS

  Parameters(const std::string& version_stamp, concurrency_type& concurrency);

  template <typename Writer>
  void write(Writer& writer);
  template <typename Reader>
  void read_input_and_broadcast(const std::string& file_name);

  void update_model();
  void update_domains();

  int get_buffer_size(const concurrency_type& concurrency) const;
  void pack(const concurrency_type& concurrency, char* buffer, int buffer_size, int& position) const;
  void unpack(const concurrency_type& concurrency, char* buffer, int buffer_size, int& position);

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

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator, solver_name>::Parameters(
    const std::string& version_stamp, concurrency_type& concurrency)
    : AnalysisParameters(),
      DcaParameters(),
      DomainsParameters(Model::DIMENSION),
      DoubleCountingParameters(),
      EdSolverParameters(),
      FourPointParameters<Model::DIMENSION>(),
      McSolverParameters<solver_name>(),
      MciParameters(),
      ModelParameters<Model>(),
      OutputParameters(),
      PhysicsParameters(),

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

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
template <typename Writer>
void Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator, solver_name>::write(
    Writer& writer) {
  writer.open_group("parameters");
  this->readWrite(writer);
  writer.close_group();

  writer.open_group("domains");

  DcaClusterFamily::write(writer);
  HostSpClusterFamily::write(writer);
  HostTpClusterFamily::write(writer);

  TDmn::parameter_type::write(writer);
  WDmn::parameter_type::write(writer);

  domains::vertex_frequency_domain<domains::EXTENDED_BOSONIC>::write(writer);

  domains::DCA_iteration_domain::write(writer);

  if (FourPointParameters<Model::DIMENSION>::get_four_point_type() != NONE) {
    domains::vertex_time_domain<domains::SP_TIME_DOMAIN>::write(writer);
    domains::vertex_time_domain<domains::TP_TIME_DOMAIN>::write(writer);
    domains::vertex_time_domain<domains::SP_TIME_DOMAIN_POSITIVE>::write(writer);
    domains::vertex_time_domain<domains::TP_TIME_DOMAIN_POSITIVE>::write(writer);

    domains::vertex_frequency_domain<domains::COMPACT>::write(writer);
    domains::vertex_frequency_domain<domains::EXTENDED>::write(writer);

    domains::vertex_frequency_domain<domains::COMPACT_POSITIVE>::write(writer);
    domains::vertex_frequency_domain<domains::EXTENDED_POSITIVE>::write(writer);

    domains::FrequencyExchangeDomain::write(writer);
    domains::MomentumExchangeDomain::write(writer);
  }

  domains::frequency_domain_real_axis::write(writer);

#ifdef DCA_WITH_QMC_BIT
  domains::numerical_error_domain::write(writer);
#endif  // DCA_WITH_QMC_BIT

  writer.close_group();
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
template <typename Reader>
void Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator,
                solver_name>::read_input_and_broadcast(const std::string& filename) {
  if (concurrency_.id() == concurrency_.first()) {
    Reader read_obj;
    read_obj.open_file(filename);
    this->readWrite(read_obj);
    read_obj.close_file();
  }

  concurrency_.broadcast_object(*this);
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
void Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator,
                solver_name>::update_model() {
  Model::initialize(*this);
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
void Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator,
                solver_name>::update_domains() {
  domains::DCA_iteration_domain::initialize(*this);
  domains::electron_band_domain::initialize(*this, Model::BANDS, Model::get_flavors(),
                                            Model::get_a_vectors());

  // time and frequency-domains
  domains::time_domain::initialize(*this);
  domains::time_domain_left_oriented::initialize(*this);
  domains::frequency_domain::initialize(*this);
  domains::frequency_domain_real_axis::initialize(*this);

  domains::vertex_time_domain<domains::SP_TIME_DOMAIN>::initialize(*this);
  domains::vertex_time_domain<domains::TP_TIME_DOMAIN>::initialize(*this);
  domains::vertex_time_domain<domains::SP_TIME_DOMAIN_POSITIVE>::initialize(*this);
  domains::vertex_time_domain<domains::TP_TIME_DOMAIN_POSITIVE>::initialize(*this);

  domains::vertex_frequency_domain<domains::COMPACT>::initialize(*this);
  domains::vertex_frequency_domain<domains::EXTENDED>::initialize(*this);

  domains::vertex_frequency_domain<domains::COMPACT_POSITIVE>::initialize(*this);
  domains::vertex_frequency_domain<domains::EXTENDED_POSITIVE>::initialize(*this);

  domains::vertex_frequency_domain<domains::EXTENDED_BOSONIC>::initialize(*this);

  domains::FrequencyExchangeDomain::initialize(*this);

  // DCA cluster
  domains::cluster_domain_initializer<RClusterDmn>::execute(Model::get_r_DCA_basis(),
                                                            DomainsParameters::get_cluster());
  domains::cluster_domain_symmetry_initializer<
      RClusterDmn, typename Model::lattice_type::DCA_point_group>::execute();

  if (concurrency_.id() == concurrency_.first())
    KClusterDmn::parameter_type::print(std::cout);

  // Host grid for single-particle functions ((sp-)lattice)
  domains::cluster_domain_initializer<RSpHostDmn>::execute(Model::get_r_DCA_basis(),
                                                           DomainsParameters::get_sp_host());
  domains::cluster_domain_symmetry_initializer<
      RSpHostDmn, typename Model::lattice_type::DCA_point_group>::execute();

  domains::MomentumExchangeDomain::initialize(*this);

  if (concurrency_.id() == concurrency_.first())
    KSpHostDmn::parameter_type::print(std::cout);

  // Host grid for two-particle functions (tp-lattice)
  if (do_dca_plus()) {
    domains::cluster_domain_initializer<RTpHostDmn>::execute(Model::get_r_DCA_basis(),
                                                             DomainsParameters::get_tp_host());
  }
  // Set equal to DCA cluster, if standard DCA is used.
  // In this way, we can keep the BseLatticeSolver general.
  else {
    domains::cluster_domain_initializer<RTpHostDmn>::execute(Model::get_r_DCA_basis(),
                                                             DomainsParameters::get_cluster());
  }
  domains::cluster_domain_symmetry_initializer<
      RTpHostDmn, typename Model::lattice_type::DCA_point_group>::execute();

  if (concurrency_.id() == concurrency_.first())
    KTpHostDmn::parameter_type::print(std::cout);
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
int Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator,
               solver_name>::get_buffer_size(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += AnalysisParameters::getBufferSize(concurrency);
  buffer_size += DcaParameters::getBufferSize(concurrency);
  buffer_size += DomainsParameters::getBufferSize(concurrency);
  buffer_size += DoubleCountingParameters::getBufferSize(concurrency);
  buffer_size += EdSolverParameters::getBufferSize(concurrency);
  buffer_size += FourPointParameters<Model::DIMENSION>::getBufferSize(concurrency);
  buffer_size += McSolverParameters<solver_name>::getBufferSize(concurrency);
  buffer_size += MciParameters::getBufferSize(concurrency);
  buffer_size += ModelParameters<Model>::getBufferSize(concurrency);
  buffer_size += OutputParameters::getBufferSize(concurrency);
  buffer_size += PhysicsParameters::getBufferSize(concurrency);

  return buffer_size;
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
void Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator, solver_name>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  AnalysisParameters::pack(concurrency, buffer, buffer_size, position);
  DcaParameters::pack(concurrency, buffer, buffer_size, position);
  DomainsParameters::pack(concurrency, buffer, buffer_size, position);
  DoubleCountingParameters::pack(concurrency, buffer, buffer_size, position);
  EdSolverParameters::pack(concurrency, buffer, buffer_size, position);
  FourPointParameters<Model::DIMENSION>::pack(concurrency, buffer, buffer_size, position);
  McSolverParameters<solver_name>::pack(concurrency, buffer, buffer_size, position);
  MciParameters::pack(concurrency, buffer, buffer_size, position);
  ModelParameters<Model>::pack(concurrency, buffer, buffer_size, position);
  OutputParameters::pack(concurrency, buffer, buffer_size, position);
  PhysicsParameters::pack(concurrency, buffer, buffer_size, position);
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
void Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator, solver_name>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  AnalysisParameters::unpack(concurrency, buffer, buffer_size, position);
  DcaParameters::unpack(concurrency, buffer, buffer_size, position);
  DomainsParameters::unpack(concurrency, buffer, buffer_size, position);
  DoubleCountingParameters::unpack(concurrency, buffer, buffer_size, position);
  EdSolverParameters::unpack(concurrency, buffer, buffer_size, position);
  FourPointParameters<Model::DIMENSION>::unpack(concurrency, buffer, buffer_size, position);
  McSolverParameters<solver_name>::unpack(concurrency, buffer, buffer_size, position);
  MciParameters::unpack(concurrency, buffer, buffer_size, position);
  ModelParameters<Model>::unpack(concurrency, buffer, buffer_size, position);
  OutputParameters::unpack(concurrency, buffer, buffer_size, position);
  PhysicsParameters::unpack(concurrency, buffer, buffer_size, position);
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
template <typename ReaderOrWriter>
void Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator, solver_name>::readWrite(
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

  AnalysisParameters::readWrite(reader_or_writer);
  DcaParameters::readWrite(reader_or_writer);
  DomainsParameters::readWrite(reader_or_writer);
  DoubleCountingParameters::readWrite(reader_or_writer);
  EdSolverParameters::readWrite(reader_or_writer);
  FourPointParameters<Model::DIMENSION>::readWrite(reader_or_writer);
  McSolverParameters<solver_name>::readWrite(reader_or_writer);
  MciParameters::readWrite(reader_or_writer);
  ModelParameters<Model>::readWrite(reader_or_writer);
  OutputParameters::readWrite(reader_or_writer);
  PhysicsParameters::readWrite(reader_or_writer);
}

template <typename Concurrency, typename Threading, typename Profiler, typename Model,
          typename RandomNumberGenerator, solver::ClusterSolverName solver_name>
std::string Parameters<Concurrency, Threading, Profiler, Model, RandomNumberGenerator,
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
