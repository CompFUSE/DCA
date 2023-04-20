// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains physical and computational data of the DCA(+) loop.

#ifndef DCA_PHYS_DCA_LOOP_DCA_LOOP_DATA_HPP
#define DCA_PHYS_DCA_LOOP_DCA_LOOP_DATA_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/filesystem.hpp"
#include "dca/io/reader.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/dca_iteration_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#ifdef DCA_HAVE_ADIOS2
#include "dca/io/adios2/adios2_writer.hpp"
#endif

namespace dca {
namespace phys {
// dca::phys::

template <typename Parameters>
class DcaLoopData {
public:
  using Concurrency = typename Parameters::concurrency_type;
  using expansion_dmn_t = func::dmn_0<func::dmn<32, int>>;
  using DCA_iteration_domain_type = func::dmn_0<domains::DCA_iteration_domain>;

  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;
  
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, Parameters::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  DcaLoopData();

  template <typename WRITER>
  void write(WRITER& writer, Concurrency& concurrency);

  // Attempts to read the loop functions from 'filename'. If successful returns
  // the last completed iteration from the input file, otherwise it returns -1.
  int readData(const std::string& filename, const std::string& format,
               const Concurrency& concurrency);
#ifdef DCA_HAVE_ADIOS2
  int readData(const std::string& filename, const std::string& format,
               const Concurrency& concurrency, adios2::ADIOS& adios);
#endif
  template <class READER>
  void readLoopDataCommon(READER& reader, const std::string& filename, const std::string& format);
  func::function<double, DCA_iteration_domain_type> Gflop_per_mpi_task;
  func::function<double, DCA_iteration_domain_type> times_per_mpi_task;

  func::function<double, DCA_iteration_domain_type> thermalization_per_mpi_task;
  func::function<double, DCA_iteration_domain_type> MC_integration_per_mpi_task;

  func::function<double, DCA_iteration_domain_type> Gflops_per_mpi_task;
  func::function<double, DCA_iteration_domain_type> max_Gflops_per_mpi_task;

  func::function<Scalar, DCA_iteration_domain_type> sign;

  func::function<double, DCA_iteration_domain_type> L2_Sigma_difference;

  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> Sigma_zero_moment;
  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> standard_deviation;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, expansion_dmn_t, DCA_iteration_domain_type>>
      sigma_lambda;

  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> n_k;
  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> A_k;
  func::function<double, func::dmn_variadic<nu, DCA_iteration_domain_type>> orbital_occupancies;

  func::function<double, DCA_iteration_domain_type> density;
  func::function<double, DCA_iteration_domain_type> chemical_potential;
  func::function<double, DCA_iteration_domain_type> average_expansion_order;

  int last_completed_iteration = -1;
};

template <typename Parameters>
DcaLoopData<Parameters>::DcaLoopData()
    : Gflop_per_mpi_task("Gflop_per_mpi_task"),
      times_per_mpi_task("times_per_mpi_task"),

      thermalization_per_mpi_task("thermalization_per_mpi_task"),
      MC_integration_per_mpi_task("MC_integration_per_mpi_task"),

      Gflops_per_mpi_task("Gflops_per_mpi_task"),
      max_Gflops_per_mpi_task("max_Gflops_per_mpi_task"),

      sign("sign"),
      L2_Sigma_difference("L2_Sigma_difference"),

      Sigma_zero_moment("Sigma_zero_moment"),
      standard_deviation("standard_deviation"),

      sigma_lambda("sigma_lambda"),

      n_k("n_k"),
      A_k("A_k"),
      orbital_occupancies("orbital-occupancies"),

      density("density"),
      chemical_potential("chemical-potential"),
      average_expansion_order("expansion_order") {}

template <typename Parameters>
template <typename Writer>
void DcaLoopData<Parameters>::write(Writer& writer, Concurrency& concurrency) {
  // This is suboptimal really the writer should report whether the file is open.
  if (concurrency.id() == concurrency.first()) {
    if (writer.isOpen()) {
      writer.open_group("DCA-loop-functions");
      {
        writer.execute(Gflop_per_mpi_task);
        writer.execute(times_per_mpi_task);
        writer.execute(Gflops_per_mpi_task);

        writer.execute(sign);

        writer.execute(L2_Sigma_difference);
        writer.execute(standard_deviation);

        writer.execute(chemical_potential);
        writer.execute(density);
        writer.execute(average_expansion_order);

        writer.execute(Sigma_zero_moment);

        writer.execute(n_k);
        writer.execute(A_k);
        writer.execute(orbital_occupancies);
        writer.execute("completed-iteration", last_completed_iteration);
      }

      writer.close_group();
    }
    else {
      throw std::logic_error(
          "Dca_loop_data write called on writer with closed file.  This is likely a developer "
          "error.");
    }
  }
}

template <typename Parameters>
int DcaLoopData<Parameters>::readData(const std::string& filename, const std::string& format,
                                          const Concurrency& concurrency) {
  if (concurrency.id() == concurrency.first() && filesystem::exists(filename)) {
    io::Reader reader(concurrency, format, false);
    readLoopDataCommon(reader, filename, format);
  }
  concurrency.broadcast(last_completed_iteration);
  return last_completed_iteration;
}

#ifdef DCA_HAVE_ADIOS2
template <typename Parameters>
int DcaLoopData<Parameters>::readData(const std::string& filename, const std::string& format,
                                          const Concurrency& concurrency, adios2::ADIOS& adios) {
  std::cout << "Reading dca_loop_data with ADIOS2\n";
  if (concurrency.id() == concurrency.first() && filesystem::exists(filename)) {
    io::Reader reader(adios, concurrency, format, false);
    auto& adios2_reader = std::get<io::ADIOS2Reader<Concurrency>>(reader.getUnderlying());
    std::size_t step_count = adios2_reader.getStepCount();
    for (std::size_t i = 0; i < step_count; ++i) {
      adios2_reader.begin_step();
      adios2_reader.end_step();
    }
    readLoopDataCommon(reader, filename, format);
  }
  concurrency.broadcast(last_completed_iteration);
  return last_completed_iteration;
}
#endif

template <typename Parameters>
template <class READER>
void DcaLoopData<Parameters>::readLoopDataCommon(READER& reader, const std::string& filename,
                                                     const std::string& format [[maybe_unused]]) {
  reader.open_file(filename);
  reader.open_group("DCA-loop-functions");

  reader.execute("completed-iteration", last_completed_iteration);

  reader.execute(Gflop_per_mpi_task);
  reader.execute(times_per_mpi_task);
  reader.execute(Gflops_per_mpi_task);

  reader.execute(sign);

  reader.execute(L2_Sigma_difference);
  reader.execute(standard_deviation);

  reader.execute(chemical_potential);
  reader.execute(density);
  reader.execute(average_expansion_order);

  reader.execute(Sigma_zero_moment);

  reader.execute(n_k);
  reader.execute(A_k);
  reader.execute(orbital_occupancies);

  reader.close_file();
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_LOOP_DCA_LOOP_DATA_HPP
