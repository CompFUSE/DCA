// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class organizes the MC walker in the CT-AUX QMC.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_WALKER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_WALKER_HPP

#include <cassert>
#include <cstdint>  // uint64_t
#include <cstdlib>  // std::size_t
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/config/threading.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/util/gpu_event.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/math/util/phase.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_vertex_move_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ct_aux_hs_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ctaux_walker_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/read_write_config.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/ct_aux_walker_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/n_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/shrink_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/walker_bit.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
class CtauxWalker : public WalkerBIT<Parameters, Data>, public CtauxWalkerData<device_t, Parameters> {
public:
  using parameters_type = Parameters;  // typedef for derived class.
  using Configuration = CT_AUX_HS_configuration<Parameters>;
  using Rng = typename Parameters::random_number_generator;

  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename CtauxTypedefs<Parameters, Data>::concurrency_type;
  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;

  constexpr static dca::linalg::DeviceType device = device_t;
  static constexpr bool is_complex = dca::util::IsComplex<Scalar>::value;

public:
  /** constructor
   *  param[in] id    thread id?
   */
  CtauxWalker(Parameters& parameters_ref, Data& MOMS_ref, Rng& rng_ref, int id);

  void initialize(int iteration);

  bool is_thermalized() const {
    return thermalized_;
  }

  void markThermalized();

  // Does one sweep, if the walker is not yet thermalized (warm-up).
  // Otherwise, does multiple sweeps according to the input parameter "sweeps-per-measurement".
  void doSweep();

  // Does one submatrix step of at most single_spin_updates_todo single spin updates.
  // Precondition: single_spin_updates_todo > 0
  // Postcondition: single_spin_updates_todo has been updated according to the executed submatrix
  //                update.
  // In/Out: single_spin_updates_todo
  void doStep(int& single_spin_updates_todo);

  /// Factored out of read_Gamma_matrices to make this computation testable.
  void calcExpVandPushVertexIndices(e_spin_states e_spin);

  /// Do significant just in time computation prior to moving data off of "device"
  void read_Gamma_matrices(e_spin_states e_spin);

  // Computes M(N).
  // Out: Ms.
  // Returns: pointer to the event marking the end of the computation.
  template <typename AccumType>
  const linalg::util::GpuEvent* computeM(std::array<linalg::Matrix<AccumType, device_t>, 2>& Ms);

  auto& get_configuration() {
    return configuration_;
  }
  const auto& get_configuration() const {
    return configuration_;
  }

  auto get_matrix_configuration() const {
    return configuration_.get_matrix_configuration();
  }

  unsigned long get_steps() const {
    return n_steps_;
  }

  auto get_sign() const;

  int get_thread_id();

  double get_Gflop();

  template <class stream_type>
  void to_JSON(stream_type& /*ss*/) {}

  void readConfig(dca::io::Buffer& buff);
  dca::io::Buffer dumpConfig() const;

  // Writes the current progress, the number of interacting spins and the total configuration_ size
  // to stdout.
  // TODO: Before this method can be made const, CT_AUX_HS_configuration and vertex_pair need to be
  //       made const correct.
  void updateShell(const int done, const int total) /*const*/;

  // Writes a summary of the walker's Markov chain updates and visited configurations to stdout.
  void printSummary() const;

  std::size_t deviceFingerprint() const {
    if (device_t == linalg::GPU)
      return G0_tools_obj.deviceFingerprint() + N_tools_obj.deviceFingerprint() +
             SHRINK_tools_obj.deviceFingerprint() +
             CtauxWalkerData<device_t, Parameters>::deviceFingerprint();
    else
      return 0;
  }

  Real get_MC_log_weight() const {
    return mc_log_weight_;
  }

  void addNonInteractingSpinsToMatrices();

  void clean_up_the_configuration();

  /// Actually doesn't seem to "compute the Gamma" but adds spins and calls solve_Gamma_blocked
  void compute_Gamma_matrices();

  template <dca::linalg::DeviceType dev_t = device_t>
  void actually_download_from_device();

  template <dca::linalg::DeviceType dev_t = device_t>
  std::enable_if_t<dev_t == device_t && device_t != dca::linalg::CPU, void> download_from_device();

  template <dca::linalg::DeviceType dev_t = device_t>
  std::enable_if_t<dev_t == device_t && device_t == dca::linalg::CPU, void> download_from_device();

  void generate_delayed_spins(int& single_spin_updates_todo);

  void update_N_matrix_with_Gamma_matrix();

  template <dca::linalg::DeviceType dev_t = device_t>
  std::enable_if_t<dev_t == device_t && device_t != dca::linalg::CPU, void> upload_to_device();

  template <dca::linalg::DeviceType dev_t = device_t>
  std::enable_if_t<dev_t == device_t && device_t == dca::linalg::CPU, void> upload_to_device();

  // for testing
  auto& getNUp() {
    return N_up;
  }
  auto& getNDown() {
    return N_dn;
  }

private:
  // Generates delayed single spin updates.
  // Returns the total number of proposed single spin updates including "static" steps.
  // Version that aborts when a Bennett spin is proposed for removal.
  int generateDelayedSpinsAbortAtBennett(int single_spin_updates_todo);

  // Version that neglects Bennett updates.
  int generateDelayedSpinsNeglectBennett(int single_spin_updates_todo);

  void finalizeDelayedSpins();

  void add_delayed_spin(int& delayed_index, int& Gamma_up_size, int& Gamma_dn_size);

  void add_delayed_spins_to_the_configuration();

  void remove_non_accepted_and_bennett_spins_from_Gamma(int& Gamma_up_size, int& Gamma_dn_size);

  void apply_bennett_on_Gamma_matrices(int& Gamma_up_size, int& Gamma_dn_size);

  void neutralize_delayed_spin(int& delayed_index, int& Gamma_up_size, int& Gamma_dn_size);

  void add_delayed_spins();

  HS_vertex_move_type get_new_HS_move();

  // INTERNAL: Unused.
  int get_new_vertex_index(HS_vertex_move_type HS_current_move);

  // INTERNAL: Unused.
  HS_spin_states_type get_new_spin_value(HS_vertex_move_type HS_current_move);

  auto calculate_acceptance_ratio(Scalar ratio, HS_vertex_move_type HS_current_move,
                                  Real QMC_factor) const;

  bool assert_exp_delta_V_value(HS_field_sign HS_field, int random_vertex_ind,
                                HS_spin_states_type new_HS_spin_value, Scalar exp_delta_V);

  void recomputeMCWeight();

  using WalkerBIT<Parameters, Data>::check_G0_matrices;
  using WalkerBIT<Parameters, Data>::check_N_matrices;
  using WalkerBIT<Parameters, Data>::check_G_matrices;

protected:
  const Parameters& parameters_;
  Data& data_;
  const Concurrency& concurrency_;

  int thread_id;
  int stream_id;

  CV<Parameters> CV_obj;

  // So actually GPU CT_AUX_WALKER_TOOLS are not currently used.
  CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar> ctaux_tools;

  Rng& rng;
  Configuration configuration_;

  G0Interpolation<device_t, Parameters> G0_tools_obj;
  N_TOOLS<device_t, Parameters> N_tools_obj;
  G_TOOLS<device_t, Parameters> G_tools_obj;

  SHRINK_TOOLS<Profiler, device_t, Scalar> SHRINK_tools_obj;

  struct delayed_spin_struct {
    int delayed_spin_index;

    HS_vertex_move_type HS_current_move;
    int random_vertex_ind;
    HS_spin_states_type new_HS_spin_value;

    Real QMC_factor;

    bool is_accepted_move;
    bool is_a_bennett_spin;

    e_spin_states e_spin_HS_field_DN;
    e_spin_states e_spin_HS_field_UP;

    int configuration_e_spin_index_HS_field_DN;
    int configuration_e_spin_index_HS_field_UP;

    int Gamma_index_HS_field_DN;
    int Gamma_index_HS_field_UP;

    Scalar exp_V_HS_field_DN;
    Scalar exp_V_HS_field_UP;

    Scalar exp_delta_V_HS_field_DN;
    Scalar exp_delta_V_HS_field_UP;

    Scalar exp_minus_delta_V_HS_field_UP;
    Scalar exp_minus_delta_V_HS_field_DN;
  };

  dca::linalg::Matrix<Scalar, dca::linalg::CPU> Gamma_up_CPU;
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> Gamma_dn_CPU;

  int warm_up_sweeps_done_;
  util::Accumulator<std::size_t> warm_up_expansion_order_;
  util::Accumulator<std::size_t> num_delayed_spins_;

  using CtauxWalkerData<device_t, Parameters>::N_up;
  using CtauxWalkerData<device_t, Parameters>::N_dn;

  using CtauxWalkerData<device_t, Parameters>::G0_up;
  using CtauxWalkerData<device_t, Parameters>::G0_dn;

  using CtauxWalkerData<device_t, Parameters>::Gamma_up;
  using CtauxWalkerData<device_t, Parameters>::Gamma_dn;

  using CtauxWalkerData<device_t, Parameters>::G_up;
  using CtauxWalkerData<device_t, Parameters>::G_dn;

  std::vector<Scalar> exp_V_CPU;
  dca::linalg::Vector<Scalar, device_t> exp_V;
  dca::linalg::Vector<Scalar, device_t> exp_delta_V;
  dca::linalg::Vector<int, device_t> vertex_indixes;

private:
  Real Gamma_up_diag_max;
  Real Gamma_up_diag_min;
  Real Gamma_dn_diag_max;
  Real Gamma_dn_diag_min;

  dca::linalg::Matrix<Scalar, dca::linalg::CPU> stored_Gamma_up_CPU;
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> stored_Gamma_dn_CPU;

  std::vector<int> random_vertex_vector;
  std::vector<HS_vertex_move_type> HS_current_move_vector;
  std::vector<HS_spin_states_type> new_HS_spin_value_vector;

  std::vector<int> vertex_indixes_CPU;
  std::vector<Scalar> exp_delta_V_CPU;

  std::vector<delayed_spin_struct> delayed_spins;
  std::vector<delayed_spin_struct> bennett_spins;

  int number_of_interacting_spins;

  int number_of_creations;
  int number_of_annihilations;

  bool annihilation_proposal_aborted_;
  uint64_t aborted_vertex_id_;

  bool thermalized_;
  bool Bennett;

  math::Phase<Scalar> phase_;

  Real mc_log_weight_ = 0.;
  const Real mc_log_weight_constant_;

  int currently_proposed_creations_ = 0;
  int currently_proposed_annihilations_ = 0;

  //  std::array<linalg::Matrix<Scalar, device_t>, 2> M_;
  std::array<linalg::Vector<Scalar, linalg::CPU>, 2> exp_v_minus_one_;
  std::array<linalg::Vector<Scalar, device_t>, 2> exp_v_minus_one_dev_;
  std::array<linalg::util::GpuEvent, 2> m_computed_events_;

  double sweeps_per_measurement_ = 1.;
  unsigned long n_steps_ = 0;
  linalg::util::GpuEvent sync_streams_event_;
  bool config_initialized_;
};

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
CtauxWalker<device_t, Parameters, Data>::CtauxWalker(Parameters& parameters_ref, Data& MOMS_ref,
                                                     Rng& rng_ref, int id)
    : WalkerBIT<Parameters, Data>(parameters_ref, MOMS_ref, id),
      CtauxWalkerData<device_t, Parameters>(parameters_ref, id),
      parameters_(parameters_ref),

      data_(MOMS_ref),
      concurrency_(parameters_.get_concurrency()),

      thread_id(id),
      stream_id(0),

      CV_obj(parameters_ref),
      ctaux_tools(CtauxWalkerData<device_t, Parameters>::MAX_VERTEX_SINGLETS *
                  parameters_.get_max_submatrix_size()),

      rng(rng_ref),

      configuration_(parameters_, rng),

      G0_tools_obj(thread_id, parameters_),
      N_tools_obj(thread_id, parameters_, CV_obj),
      G_tools_obj(thread_id, parameters_, CV_obj),

      SHRINK_tools_obj(thread_id),

      Gamma_up_CPU("Gamma_up_CPU", Gamma_up.size(), Gamma_up.capacity()),
      Gamma_dn_CPU("Gamma_dn_CPU", Gamma_dn.size(), Gamma_dn.capacity()),

      warm_up_sweeps_done_(0),
      warm_up_expansion_order_(),
      num_delayed_spins_(),

      Gamma_up_diag_max(1),
      Gamma_up_diag_min(1),
      Gamma_dn_diag_max(1),
      Gamma_dn_diag_min(1),

      stored_Gamma_up_CPU("stored_Gamma_up_CPU", Gamma_up.size(), Gamma_up.capacity()),
      stored_Gamma_dn_CPU("stored_Gamma_dn_CPU", Gamma_dn.size(), Gamma_dn.capacity()),

      random_vertex_vector(0),
      HS_current_move_vector(0),
      new_HS_spin_value_vector(0),

      annihilation_proposal_aborted_(false),
      aborted_vertex_id_(0),

      thermalized_(false),
      Bennett(false),
      mc_log_weight_constant_(
          std::log(parameters_.get_expansion_parameter_K() / (2. * parameters_.get_beta()))),

      config_initialized_(false) {
  if (concurrency_.id() == 0 and thread_id == 0) {
    std::cout << "\n\n"
              << "\t\t"
              << "CT-AUX walker test"
              << "\n\n";
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::printSummary() const {
  // std::defaultfloat is only supported by GCC 5 or later.
  std::cout.unsetf(std::ios_base::floatfield);
  std::cout << "\n"
            << "Walker: process ID = " << concurrency_.id() << ", thread ID = " << thread_id << "\n"
            << "-------------------------------------------\n";

  if (warm_up_expansion_order_.count())
    std::cout << "estimate for sweep size: " << warm_up_expansion_order_.mean() << "\n";
  if (num_delayed_spins_.count())
    std::cout << "average number of delayed spins: " << num_delayed_spins_.mean() << "\n";

  std::cout << "# creations / # annihilations: "
            << static_cast<Real>(number_of_creations) / static_cast<Real>(number_of_annihilations)
            << "\n"
            << std::endl;

  std::cout << std::scientific;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
auto CtauxWalker<device_t, Parameters, Data>::get_sign() const {
  return phase_.getSign();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
int CtauxWalker<device_t, Parameters, Data>::get_thread_id() {
  return thread_id;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
double CtauxWalker<device_t, Parameters, Data>::get_Gflop() {
  Real Gflop = 0.;

  Gflop += N_tools_obj.get_Gflop();
  Gflop += G_tools_obj.get_Gflop();

  return Gflop;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::markThermalized() {
  thermalized_ = true;
  recomputeMCWeight();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::initialize(int iteration) {
  WalkerBIT<Parameters, Data>::initialize();
  sweeps_per_measurement_ = parameters_.get_sweeps_per_measurement().at(iteration);

  number_of_creations = 0;
  number_of_annihilations = 0;

  annihilation_proposal_aborted_ = false;
  // aborted_vertex_id_ = 0;

  phase_.reset();

  CV_obj.initialize(data_);

  if (!config_initialized_)
    configuration_.initialize();
  // configuration_.print();

  thermalized_ = false;

  // TODO: Reset accumulators of warm-up expansion order and number of delayed spins, and set
  //       warm_up_sweeps_done_ to zero?

  {
    // std::cout << "\n\n\t G0-TOOLS \n\n";

    G0_tools_obj.initialize(data_);
    G0_tools_obj.build_G0_matrix(configuration_, G0_up, e_UP);
    G0_tools_obj.build_G0_matrix(configuration_, G0_dn, e_DN);

#ifdef DCA_WITH_QMC_BIT
    check_G0_matrices(configuration_, G0_up, G0_dn);
#endif  // DCA_WITH_QMC_BIT
  }

  {
    // std::cout << "\n\n\t N-TOOLS \n\n";

    N_tools_obj.build_N_matrix(configuration_, N_up, G0_up, e_UP);
    N_tools_obj.build_N_matrix(configuration_, N_dn, G0_dn, e_DN);

#ifdef DCA_WITH_QMC_BIT
    check_N_matrices(configuration_, G0_up, G0_dn, N_up, N_dn);
#endif  // DCA_WITH_QMC_BIT
  }

  {
    // std::cout << "\n\n\t G-TOOLS \n\n";

    G_tools_obj.build_G_matrix(configuration_, N_up, G0_up, G_up, e_UP);
    G_tools_obj.build_G_matrix(configuration_, N_dn, G0_dn, G_dn, e_DN);

#ifdef DCA_WITH_QMC_BIT
    check_G_matrices(configuration_, G0_up, G0_dn, N_up, N_dn, G_up, G_dn);
#endif  // DCA_WITH_QMC_BIT
  }

  recomputeMCWeight();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::doSweep() {
  Profiler profiler("do_sweep", "CT-AUX walker", __LINE__, thread_id);
  const double sweeps_per_measurement{thermalized_ ? sweeps_per_measurement_ : 1.};

  // Do at least one single spin update per sweep.
  const int single_spin_updates_per_sweep{warm_up_expansion_order_.count() > 0 &&
                                                  warm_up_expansion_order_.mean() > 1.
                                              ? static_cast<int>(warm_up_expansion_order_.mean())
                                              : 1};

  // Reset the warm-up expansion order accumulator after half the warm-up sweeps to get a better
  // estimate for the expansion order of the thermalized system.
  if (warm_up_sweeps_done_ == parameters_.get_warm_up_sweeps() / 2)
    warm_up_expansion_order_.reset();

  int single_spin_updates_todo{single_spin_updates_per_sweep *
                               static_cast<int>(sweeps_per_measurement)};

  while (single_spin_updates_todo > 0) {
    doStep(single_spin_updates_todo);
  }

  assert(single_spin_updates_todo == 0);

  if (!thermalized_)
    ++warm_up_sweeps_done_;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::doStep(int& single_spin_updates_todo) {
  configuration_.prepare_configuration();

  generate_delayed_spins(single_spin_updates_todo);

  addNonInteractingSpinsToMatrices();

  download_from_device();

  // This is strange since CTAUX_TOOLS::computeGamma gets called as a side effect of the download_from_device!
  compute_Gamma_matrices();

  upload_to_device();

  update_N_matrix_with_Gamma_matrix();

  clean_up_the_configuration();

  if (!thermalized_)
    warm_up_expansion_order_.addSample(configuration_.get_number_of_interacting_HS_spins());
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <dca::linalg::DeviceType dev_t>
void CtauxWalker<device_t, Parameters, Data>::actually_download_from_device() {
  if constexpr (device_t == dca::linalg::GPU) {
    Gamma_up_CPU.setAsync(Gamma_up, thread_id, stream_id);
    Gamma_dn_CPU.setAsync(Gamma_dn, thread_id, stream_id);
    linalg::util::syncStream(thread_id, stream_id);
  }
  else if constexpr (device_t == dca::linalg::CPU) {
    Gamma_up_CPU.swap(Gamma_up);
    Gamma_dn_CPU.swap(Gamma_dn);
  }
}

// In case Gamma_up and Gamma_down do not reside in the CPU memory, copy them.
template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <dca::linalg::DeviceType dev_t>
std::enable_if_t<dev_t == device_t && device_t != dca::linalg::CPU, void> CtauxWalker<
    device_t, Parameters, Data>::download_from_device() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  calcExpVandPushVertexIndices(e_UP);
  read_Gamma_matrices(e_UP);
  calcExpVandPushVertexIndices(e_DN);
  read_Gamma_matrices(e_DN);

  actually_download_from_device<device_t>();
}

// In case Gamma_up and Gamma_down reside in the CPU memory, avoid the copies using swap.
template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <dca::linalg::DeviceType dev_t>
std::enable_if_t<dev_t == device_t && device_t == dca::linalg::CPU, void> CtauxWalker<
    device_t, Parameters, Data>::download_from_device() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  assert(Gamma_up_CPU.capacity() == Gamma_up.capacity());
  assert(Gamma_dn_CPU.capacity() == Gamma_dn.capacity());

  calcExpVandPushVertexIndices(e_UP);
  read_Gamma_matrices(e_UP);
  calcExpVandPushVertexIndices(e_DN);
  read_Gamma_matrices(e_DN);

  actually_download_from_device<device_t>();
}

// In case Gamma_up and Gamma_down do not reside in the CPU memory, copy them.
template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <dca::linalg::DeviceType dev_t>
std::enable_if_t<dev_t == device_t && device_t != dca::linalg::CPU, void> CtauxWalker<
    device_t, Parameters, Data>::upload_to_device() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  Gamma_up.setAsync(Gamma_up_CPU, thread_id, stream_id);
  Gamma_dn.setAsync(Gamma_dn_CPU, thread_id, stream_id);

  linalg::util::syncStream(thread_id, stream_id);
}

// In case Gamma_up and Gamma_down reside in the CPU memory, avoid the copies using swap.
template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <dca::linalg::DeviceType dev_t>
std::enable_if_t<dev_t == device_t && device_t == dca::linalg::CPU, void> CtauxWalker<
    device_t, Parameters, Data>::upload_to_device() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  assert(Gamma_up_CPU.capacity() == Gamma_up.capacity());
  assert(Gamma_dn_CPU.capacity() == Gamma_dn.capacity());

  Gamma_up.swap(Gamma_up_CPU);
  Gamma_dn.swap(Gamma_dn_CPU);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::addNonInteractingSpinsToMatrices() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  Gamma_up.resizeNoCopy(0);
  Gamma_dn.resizeNoCopy(0);

  {  // update G0 for new shuffled vertices
    Profiler p("G0-matrix (update)", "CT-AUX walker", __LINE__, thread_id);

    G0_tools_obj.update_G0_matrix(configuration_, G0_up, e_UP);
    G0_tools_obj.update_G0_matrix(configuration_, G0_dn, e_DN);

#ifdef DCA_WITH_QMC_BIT
    check_G0_matrices(configuration_, G0_up, G0_dn);
#endif  // DCA_WITH_QMC_BIT
  }
  {  // update N for new shuffled vertices
    Profiler p("N-matrix (update)", "CT-AUX walker", __LINE__, thread_id);

    N_tools_obj.update_N_matrix(configuration_, G0_up, N_up, e_UP);
    N_tools_obj.update_N_matrix(configuration_, G0_dn, N_dn, e_DN);

#ifdef DCA_WITH_QMC_BIT
    check_N_matrices(configuration_, G0_up, G0_dn, N_up, N_dn);
#endif  // DCA_WITH_QMC_BIT
  }

  {  // update N for new shuffled vertices
    Profiler p("G-matrix (update)", "CT-AUX walker", __LINE__, thread_id);

    G_tools_obj.build_G_matrix(configuration_, N_up, G0_up, G_up, e_UP);
    G_tools_obj.build_G_matrix(configuration_, N_dn, G0_dn, G_dn, e_DN);

#ifdef DCA_WITH_QMC_BIT
    check_G_matrices(configuration_, G0_up, G0_dn, N_up, N_dn, G_up, G_dn);
#endif  // DCA_WITH_QMC_BIT
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::generate_delayed_spins(int& single_spin_updates_todo) {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  assert(single_spin_updates_todo > 0);

  const auto single_spin_updates_proposed =
      parameters_.neglect_bennett_updates()
          ? generateDelayedSpinsNeglectBennett(single_spin_updates_todo)
          : generateDelayedSpinsAbortAtBennett(single_spin_updates_todo);

  n_steps_ += single_spin_updates_proposed;

  single_spin_updates_todo -= single_spin_updates_proposed;
  assert(single_spin_updates_todo >= 0);

  if (thermalized_)
    num_delayed_spins_.addSample(delayed_spins.size());

  finalizeDelayedSpins();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
int CtauxWalker<device_t, Parameters, Data>::generateDelayedSpinsAbortAtBennett(
    const int single_spin_updates_todo) {
  assert(single_spin_updates_todo > 0);

  const auto max_num_delayed_spins = parameters_.get_max_submatrix_size();

  delayed_spins.resize(0);

  int num_statics = 0;
  int single_spin_updates_proposed = 0;

  currently_proposed_annihilations_ = 0;
  currently_proposed_creations_ = 0;

  // Do the aborted annihilation proposal.
  if (annihilation_proposal_aborted_) {
    delayed_spin_struct delayed_spin;
    delayed_spin.is_accepted_move = false;
    delayed_spin.is_a_bennett_spin = false;

    const std::size_t aborted_vertex_index = configuration_.find(aborted_vertex_id_);

    // The vertex that caused the annihilation proposal to be aborted is still in the
    // configuration_. Propose removal of this vertex again.
    if (aborted_vertex_index < configuration_.size()) {
      delayed_spin.HS_current_move = ANNIHILATION;
      delayed_spin.random_vertex_ind = aborted_vertex_index;
      delayed_spin.new_HS_spin_value = HS_ZERO;

      delayed_spins.push_back(delayed_spin);
      ++currently_proposed_annihilations_;
    }

    // Propose removal of a different vertex or do a static step if the configuration_ is empty.
    else {
      delayed_spin.HS_current_move =
          configuration_.get_number_of_interacting_HS_spins() == 0 ? STATIC : ANNIHILATION;

      if (delayed_spin.HS_current_move == ANNIHILATION) {
        delayed_spin.random_vertex_ind = configuration_.get_random_interacting_vertex();
        delayed_spin.new_HS_spin_value = HS_ZERO;

        delayed_spins.push_back(delayed_spin);
        ++currently_proposed_annihilations_;
      }

      else {
        ++num_statics;
      }
    }

    ++single_spin_updates_proposed;
    annihilation_proposal_aborted_ = false;
  }

  // Generate more delayed spins.
  while (!annihilation_proposal_aborted_ && single_spin_updates_proposed < single_spin_updates_todo &&
         delayed_spins.size() < max_num_delayed_spins) {
    delayed_spin_struct delayed_spin;
    delayed_spin.is_accepted_move = false;
    delayed_spin.is_a_bennett_spin = false;
    delayed_spin.HS_current_move = get_new_HS_move();

    if (delayed_spin.HS_current_move == ANNIHILATION) {
      delayed_spin.random_vertex_ind = configuration_.get_random_interacting_vertex();
      delayed_spin.new_HS_spin_value = HS_ZERO;

      // Check whether the spin has already been changed, i.e. it has already been proposed for
      // removal, if it is an interacting spin, or it is a "virtual" interacting spin
      // (= non-interacting spin proposed for insertion).
      // If we encounter such a spin, we stop. The aborted annihilation proposal is then done in the
      // next submatrix step.
      for (const auto& other_spin : delayed_spins) {
        if (delayed_spin.random_vertex_ind == other_spin.random_vertex_ind) {
          annihilation_proposal_aborted_ = true;
          aborted_vertex_id_ = configuration_[delayed_spin.random_vertex_ind].get_id();
          break;
        }
      }

      if (!annihilation_proposal_aborted_) {
        delayed_spins.push_back(delayed_spin);
        ++currently_proposed_annihilations_;
        ++single_spin_updates_proposed;
      }
    }

    else if (delayed_spin.HS_current_move == CREATION) {
      delayed_spin.random_vertex_ind = configuration_.size();
      configuration_.insert_random_noninteracting_vertex(true);
      delayed_spin.new_HS_spin_value = rng() > 0.5 ? HS_UP : HS_DN;

      delayed_spins.push_back(delayed_spin);
      ++currently_proposed_creations_;
      ++single_spin_updates_proposed;
    }

    else {
      assert(delayed_spin.HS_current_move == STATIC);
      ++num_statics;
      ++single_spin_updates_proposed;
    }
  }

  // We need to unmark all "virtual" interacting spins, that we have temporarily marked as
  // annihilatable in CT_AUX_HS_configuration::get_random_noninteracting_vertex().
  // TODO: Eliminate the need to mark and unmark these spins.
  for (const auto& spin : delayed_spins)
    if (spin.HS_current_move == CREATION)
      configuration_.unmarkAsAnnihilatable(spin.random_vertex_ind);

  assert(single_spin_updates_proposed ==
         currently_proposed_creations_ + currently_proposed_annihilations_ + num_statics);

  return single_spin_updates_proposed;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
int CtauxWalker<device_t, Parameters, Data>::generateDelayedSpinsNeglectBennett(
    const int single_spin_updates_todo) {
  assert(single_spin_updates_todo > 0);

  const auto max_num_delayed_spins = parameters_.get_max_submatrix_size();
  const auto num_interacting_spins_initial = configuration_.get_number_of_interacting_HS_spins();

  delayed_spins.resize(0);

  int single_spin_updates_proposed = 0;
  int num_statics = 0;

  currently_proposed_annihilations_ = 0;
  currently_proposed_creations_ = 0;

  while ((num_interacting_spins_initial == 0 ||
          currently_proposed_annihilations_ < num_interacting_spins_initial) &&
         single_spin_updates_proposed < single_spin_updates_todo &&
         delayed_spins.size() < max_num_delayed_spins) {
    delayed_spin_struct delayed_spin;
    delayed_spin.is_accepted_move = false;
    delayed_spin.is_a_bennett_spin = false;
    delayed_spin.HS_current_move = get_new_HS_move();

    if (delayed_spin.HS_current_move == ANNIHILATION) {
      delayed_spin.new_HS_spin_value = HS_ZERO;

      bool has_already_been_chosen = true;

      while (has_already_been_chosen) {
        delayed_spin.random_vertex_ind = configuration_.get_random_interacting_vertex();
        has_already_been_chosen = false;

        for (const auto& other_spin : delayed_spins)
          if (delayed_spin.random_vertex_ind == other_spin.random_vertex_ind) {
            has_already_been_chosen = true;
            break;
          }
      }

      delayed_spins.push_back(delayed_spin);
      ++currently_proposed_annihilations_;
      ++single_spin_updates_proposed;
    }

    else if (delayed_spin.HS_current_move == CREATION) {
      delayed_spin.random_vertex_ind = configuration_.size();
      configuration_.insert_random_noninteracting_vertex(false);
      delayed_spin.new_HS_spin_value = rng() > 0.5 ? HS_UP : HS_DN;

      delayed_spins.push_back(delayed_spin);
      ++currently_proposed_creations_;
      ++single_spin_updates_proposed;
    }

    else {
      assert(delayed_spin.HS_current_move == STATIC);
      ++num_statics;
      ++single_spin_updates_proposed;
    }
  }

  assert(single_spin_updates_proposed ==
         currently_proposed_creations_ + currently_proposed_annihilations_ + num_statics);

  return single_spin_updates_proposed;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::finalizeDelayedSpins() {
  int Gamma_dn_size = 0;
  int Gamma_up_size = 0;

  for (std::size_t i = 0; i < delayed_spins.size(); ++i) {
    delayed_spins[i].delayed_spin_index = i;

    const int random_vertex_ind = delayed_spins[i].random_vertex_ind;
    const HS_spin_states_type new_HS_spin_value = delayed_spins[i].new_HS_spin_value;

    delayed_spins[i].QMC_factor =
        CV_obj.get_QMC_factor(configuration_[random_vertex_ind], new_HS_spin_value);

    delayed_spins[i].e_spin_HS_field_DN = configuration_[random_vertex_ind].get_e_spins().first;
    delayed_spins[i].e_spin_HS_field_UP = configuration_[random_vertex_ind].get_e_spins().second;

    delayed_spins[i].configuration_e_spin_index_HS_field_DN =
        configuration_[random_vertex_ind].get_configuration_e_spin_indices().first;
    delayed_spins[i].configuration_e_spin_index_HS_field_UP =
        configuration_[random_vertex_ind].get_configuration_e_spin_indices().second;

    if (delayed_spins[i].e_spin_HS_field_DN == e_UP) {
      delayed_spins[i].Gamma_index_HS_field_DN = Gamma_up_size++;

      const vertex_singleton& v_j =
          configuration_.get(e_UP)[delayed_spins[i].configuration_e_spin_index_HS_field_DN];

      delayed_spins[i].exp_V_HS_field_DN = CV_obj.exp_V(v_j);
      delayed_spins[i].exp_delta_V_HS_field_DN = CV_obj.exp_delta_V(v_j, new_HS_spin_value);
      delayed_spins[i].exp_minus_delta_V_HS_field_DN =
          CV_obj.exp_minus_delta_V(v_j, new_HS_spin_value);
    }

    else {
      delayed_spins[i].Gamma_index_HS_field_DN = Gamma_dn_size++;

      const vertex_singleton& v_j =
          configuration_.get(e_DN)[delayed_spins[i].configuration_e_spin_index_HS_field_DN];

      delayed_spins[i].exp_V_HS_field_DN = CV_obj.exp_V(v_j);
      delayed_spins[i].exp_delta_V_HS_field_DN = CV_obj.exp_delta_V(v_j, new_HS_spin_value);
      delayed_spins[i].exp_minus_delta_V_HS_field_DN =
          CV_obj.exp_minus_delta_V(v_j, new_HS_spin_value);
    }

    if (delayed_spins[i].e_spin_HS_field_UP == e_UP) {
      delayed_spins[i].Gamma_index_HS_field_UP = Gamma_up_size++;

      const vertex_singleton& v_j =
          configuration_.get(e_UP)[delayed_spins[i].configuration_e_spin_index_HS_field_UP];

      delayed_spins[i].exp_V_HS_field_UP = CV_obj.exp_V(v_j);
      delayed_spins[i].exp_delta_V_HS_field_UP = CV_obj.exp_delta_V(v_j, new_HS_spin_value);
      delayed_spins[i].exp_minus_delta_V_HS_field_UP =
          CV_obj.exp_minus_delta_V(v_j, new_HS_spin_value);
    }

    else {
      delayed_spins[i].Gamma_index_HS_field_UP = Gamma_dn_size++;

      const vertex_singleton& v_j =
          configuration_.get(e_DN)[delayed_spins[i].configuration_e_spin_index_HS_field_UP];

      delayed_spins[i].exp_V_HS_field_UP = CV_obj.exp_V(v_j);
      delayed_spins[i].exp_delta_V_HS_field_UP = CV_obj.exp_delta_V(v_j, new_HS_spin_value);
      delayed_spins[i].exp_minus_delta_V_HS_field_UP =
          CV_obj.exp_minus_delta_V(v_j, new_HS_spin_value);
    }
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::calcExpVandPushVertexIndices(e_spin_states e_spin) {
  exp_V_CPU.resize(0);
  exp_delta_V_CPU.resize(0);
  vertex_indixes_CPU.resize(0);

  for (size_t i = 0; i < delayed_spins.size(); ++i) {
    if (delayed_spins[i].e_spin_HS_field_DN == e_spin) {
      exp_V_CPU.push_back(delayed_spins[i].exp_V_HS_field_DN);
      exp_delta_V_CPU.push_back(delayed_spins[i].exp_delta_V_HS_field_DN);
      vertex_indixes_CPU.push_back(delayed_spins[i].configuration_e_spin_index_HS_field_DN);
    }

    if (delayed_spins[i].e_spin_HS_field_UP == e_spin) {
      exp_V_CPU.push_back(delayed_spins[i].exp_V_HS_field_UP);
      exp_delta_V_CPU.push_back(delayed_spins[i].exp_delta_V_HS_field_UP);
      vertex_indixes_CPU.push_back(delayed_spins[i].configuration_e_spin_index_HS_field_UP);
    }
  }

  vertex_indixes.setAsync(vertex_indixes_CPU, thread_id, stream_id);
  exp_V.setAsync(exp_V_CPU, thread_id, stream_id);
  exp_delta_V.setAsync(exp_delta_V_CPU, thread_id, stream_id);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::read_Gamma_matrices(e_spin_states e_spin) {
  // std::cout << __FUNCTION__ << "\n";

  // Profiler profiler(concurrency_, __FUNCTION__, "CT-AUX walker", __LINE__, thread_id);
  switch (e_spin) {
    case e_DN:
      linalg::util::syncStream(thread_id, stream_id);
      CT_AUX_WALKER_TOOLS<device_t, Scalar>::compute_Gamma(
          Gamma_dn, N_dn, G_dn, vertex_indixes, exp_V, exp_delta_V, thread_id, stream_id);
      // assume we've no guarantee this will be allowed to finish before the async copy starts
      linalg::util::syncStream(thread_id, stream_id);
      break;

    case e_UP:
      linalg::util::syncStream(thread_id, stream_id);
      CT_AUX_WALKER_TOOLS<device_t, Scalar>::compute_Gamma(
          Gamma_up, N_up, G_up, vertex_indixes, exp_V, exp_delta_V, thread_id, stream_id);
      // assume we've no guarantee this will be allowed to finish before the async copy starts
      linalg::util::syncStream(thread_id, stream_id);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::compute_Gamma_matrices() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  bennett_spins.resize(0);

  int Gamma_up_size = 0;
  int Gamma_dn_size = 0;

  Gamma_up_diag_max = 1;
  Gamma_up_diag_min = 1;
  Gamma_dn_diag_max = 1;
  Gamma_dn_diag_min = 1;

  number_of_interacting_spins = configuration_.get_number_of_interacting_HS_spins();

  for (int delayed_index = 0; delayed_index < int(delayed_spins.size()); delayed_index++) {
    if (delayed_spins[delayed_index].HS_current_move == ANNIHILATION) {
      Real alpha = 0;
      if (number_of_interacting_spins > 0)
        alpha = Real(bennett_spins.size()) / Real(number_of_interacting_spins);

      // INTERNAL: Does this turn off the Bennett updates?
      // TODO: Clean this up by e.g. using a flag 'Bennett_updates'.
      if (false && rng() < alpha) {
        apply_bennett_on_Gamma_matrices(Gamma_up_size, Gamma_dn_size);

        neutralize_delayed_spin(delayed_index, Gamma_up_size, Gamma_dn_size);
      }
      else {
        add_delayed_spin(delayed_index, Gamma_up_size, Gamma_dn_size);
      }
    }
    else {
      add_delayed_spin(delayed_index, Gamma_up_size, Gamma_dn_size);
    }
  }

  add_delayed_spins_to_the_configuration();

  remove_non_accepted_and_bennett_spins_from_Gamma(Gamma_up_size, Gamma_dn_size);

  // #ifdef DCA_WITH_QMC_BIT
  //   if (concurrency_.id() == 0 and thread_id == 0)
  //     std::cout << "\n\t Gamma-update check : \n\n";

  //   GAMMA_tools_obj.check_Gamma_LU(Gamma_up_CPU, N_up, G_up, configuration_, e_UP);
  //   GAMMA_tools_obj.check_Gamma_LU(Gamma_dn_CPU, N_dn, G_dn, configuration_, e_DN);
  // #endif  // DCA_WITH_QMC_BIT
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::neutralize_delayed_spin(int& delayed_index,
                                                                      int& Gamma_up_size,
                                                                      int& Gamma_dn_size) {
  // std::cout << __FUNCTION__ << "\n";

  delayed_spins[delayed_index].is_accepted_move = false;

  if (delayed_spins[delayed_index].e_spin_HS_field_DN == e_UP) {
    Gamma_up_size += 1;
    CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(
        Gamma_up_CPU, delayed_spins[delayed_index].Gamma_index_HS_field_DN);
  }
  else {
    Gamma_dn_size += 1;
    CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(
        Gamma_dn_CPU, delayed_spins[delayed_index].Gamma_index_HS_field_DN);
  }

  if (delayed_spins[delayed_index].e_spin_HS_field_UP == e_UP) {
    Gamma_up_size += 1;
    CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(
        Gamma_up_CPU, delayed_spins[delayed_index].Gamma_index_HS_field_UP);
  }
  else {
    Gamma_dn_size += 1;
    CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(
        Gamma_dn_CPU, delayed_spins[delayed_index].Gamma_index_HS_field_UP);
  }

  //     delayed_index += 1;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::add_delayed_spins_to_the_configuration() {
  // std::cout << __FUNCTION__ << "\n";

  for (size_t i = 0; i < delayed_spins.size(); ++i) {
    int configuration_index = delayed_spins[i].random_vertex_ind;

    if (delayed_spins[i].is_accepted_move) {
      if (not delayed_spins[i].is_a_bennett_spin) {
        configuration_.add_delayed_HS_spin(configuration_index, delayed_spins[i].new_HS_spin_value);
      }
      else {
        configuration_[configuration_index].set_annihilatable(false);
      }
    }
  }

  assert(number_of_interacting_spins == configuration_.get_number_of_interacting_HS_spins());
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::remove_non_accepted_and_bennett_spins_from_Gamma(
    int& Gamma_up_size, int& Gamma_dn_size) {
  // std::cout << __FUNCTION__ << "\n";

  for (int i = delayed_spins.size() - 1; i > -1; i--) {
    if ((not delayed_spins[i].is_accepted_move) or delayed_spins[i].is_a_bennett_spin) {
      e_spin_states e_spin_HS_field_DN = delayed_spins[i].e_spin_HS_field_DN;
      e_spin_states e_spin_HS_field_UP = delayed_spins[i].e_spin_HS_field_UP;

      int Gamma_index_HS_field_UP = delayed_spins[i].Gamma_index_HS_field_UP;
      int Gamma_index_HS_field_DN = delayed_spins[i].Gamma_index_HS_field_DN;

      if (e_spin_HS_field_UP == e_UP) {
        Gamma_up_size -= 1;
        dca::linalg::matrixop::removeRowAndCol(Gamma_up_CPU, Gamma_index_HS_field_UP);
      }
      else {
        Gamma_dn_size -= 1;
        dca::linalg::matrixop::removeRowAndCol(Gamma_dn_CPU, Gamma_index_HS_field_UP);
      }

      if (e_spin_HS_field_DN == e_UP) {
        Gamma_up_size -= 1;
        dca::linalg::matrixop::removeRowAndCol(Gamma_up_CPU, Gamma_index_HS_field_DN);
      }
      else {
        Gamma_dn_size -= 1;
        dca::linalg::matrixop::removeRowAndCol(Gamma_dn_CPU, Gamma_index_HS_field_DN);
      }
    }
  }
  linalg::util::syncStream(thread_id, stream_id);

  assert(Gamma_up_size == Gamma_up_CPU.size().first and Gamma_up_size == Gamma_up_CPU.size().second);
  assert(Gamma_dn_size == Gamma_dn_CPU.size().first and Gamma_dn_size == Gamma_dn_CPU.size().second);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::add_delayed_spin(int& delayed_index, int& Gamma_up_size,
                                                               int& Gamma_dn_size) {
  // std::cout << __FUNCTION__ << "\t|";

  assert(assert_exp_delta_V_value(HS_FIELD_DN, delayed_spins[delayed_index].random_vertex_ind,
                                  delayed_spins[delayed_index].new_HS_spin_value,
                                  delayed_spins[delayed_index].exp_delta_V_HS_field_DN));

  assert(assert_exp_delta_V_value(HS_FIELD_UP, delayed_spins[delayed_index].random_vertex_ind,
                                  delayed_spins[delayed_index].new_HS_spin_value,
                                  delayed_spins[delayed_index].exp_delta_V_HS_field_UP));

  if (delayed_spins[delayed_index].HS_current_move == CREATION)
    number_of_creations += 1;

  if (delayed_spins[delayed_index].HS_current_move == ANNIHILATION)
    number_of_annihilations += 1;

  int Gamma_index_HS_field_DN = delayed_spins[delayed_index].Gamma_index_HS_field_DN;  //-1;
  int Gamma_index_HS_field_UP = delayed_spins[delayed_index].Gamma_index_HS_field_UP;  //-1;

  Scalar exp_delta_V_HS_field_DN = delayed_spins[delayed_index].exp_delta_V_HS_field_DN;
  Scalar exp_delta_V_HS_field_UP = delayed_spins[delayed_index].exp_delta_V_HS_field_UP;

  Scalar ratio_HS_field_DN = 0.;
  Scalar ratio_HS_field_UP = 0.;

  Real tmp_up_diag_max = Gamma_up_diag_max;
  Real tmp_up_diag_min = Gamma_up_diag_min;
  Real tmp_dn_diag_max = Gamma_dn_diag_max;
  Real tmp_dn_diag_min = Gamma_dn_diag_min;

  {
    if (delayed_spins[delayed_index].e_spin_HS_field_DN == e_UP) {
      // std::cout << "\t" << "e_UP" << "\t" << Gamma_index_HS_field_DN << "\t" << Gamma_up_size <<
      // "\t|\t";

      assert(Gamma_index_HS_field_DN == Gamma_up_size);

      ratio_HS_field_DN = ctaux_tools.solve_Gamma_blocked(Gamma_index_HS_field_DN, Gamma_up_CPU,
                                                          exp_delta_V_HS_field_DN,
                                                          Gamma_up_diag_max, Gamma_up_diag_min);

      Gamma_up_size += 1;
    }
    else {
      // std::cout << "\t" << "e_DN" << "\t" << Gamma_index_HS_field_DN << "\t" << Gamma_dn_size <<
      // "\t|\t";

      assert(Gamma_index_HS_field_DN == Gamma_dn_size);

      ratio_HS_field_DN = ctaux_tools.solve_Gamma_blocked(Gamma_index_HS_field_DN, Gamma_dn_CPU,
                                                          exp_delta_V_HS_field_DN,
                                                          Gamma_dn_diag_max, Gamma_dn_diag_min);

      Gamma_dn_size += 1;
    }

    if (delayed_spins[delayed_index].e_spin_HS_field_UP == e_UP) {
      // std::cout << "\t" << "e_UP" << "\t" << Gamma_index_HS_field_UP << "\t" << Gamma_up_size <<
      // "\t|\t";

      assert(Gamma_index_HS_field_UP == Gamma_up_size);

      ratio_HS_field_UP = ctaux_tools.solve_Gamma_blocked(Gamma_index_HS_field_UP, Gamma_up_CPU,
                                                          exp_delta_V_HS_field_UP,
                                                          Gamma_up_diag_max, Gamma_up_diag_min);

      Gamma_up_size += 1;
    }
    else {
      // std::cout << "\t" << "e_DN" << "\t" << Gamma_index_HS_field_UP << "\t" << Gamma_dn_size <<
      // "\t|\t";

      assert(Gamma_index_HS_field_UP == Gamma_dn_size);

      ratio_HS_field_UP = ctaux_tools.solve_Gamma_blocked(Gamma_index_HS_field_UP, Gamma_dn_CPU,
                                                          exp_delta_V_HS_field_UP,
                                                          Gamma_dn_diag_max, Gamma_dn_diag_min);

      Gamma_dn_size += 1;
    }

    // std::cout << "\n";
  }

  const auto determinant_ratio = ratio_HS_field_UP * ratio_HS_field_DN;
  // std::cout << "ratio up: " << std::setprecision(16) << ratio_HS_field_UP
  //          << "  ratio down: " << std::setprecision(16) << ratio_HS_field_DN << '\n';
  // std::cout << "determinant_ratio: " << determinant_ratio << '\n';
  auto acceptance_ratio =
      calculate_acceptance_ratio(determinant_ratio, delayed_spins[delayed_index].HS_current_move,
                                 delayed_spins[delayed_index].QMC_factor);

  // std::cout << "determinant_ratio: " << determinant_ratio << '\n';
  // if constexpr (dca::util::IsComplex_t<decltype(acceptance_ratio)>::value)
  //   acceptance_ratio.imag(0.0);

  // std::cout << "acceptance_ratio: " << acceptance_ratio << '\n';

  if (std::abs(acceptance_ratio) >= rng()) {
    delayed_spins[delayed_index].is_accepted_move = true;

    // Update Monte Carlo weight.
    mc_log_weight_ += std::log(std::abs(determinant_ratio));
    // std::cout << "mc_log_weight: " << mc_log_weight_ << '\n';
    if (delayed_spins[delayed_index].HS_current_move == CREATION)
      mc_log_weight_ += mc_log_weight_constant_;
    else
      mc_log_weight_ -= mc_log_weight_constant_;

    // std::cout << "phase: " << phase_.getSign() << '\n';
    phase_.multiply(determinant_ratio);
    // std::cout << "phase: " << phase_.getSign() << '\n';

    assert(delayed_spins[delayed_index].delayed_spin_index == delayed_index);

    if (delayed_spins[delayed_index].HS_current_move == CREATION) {
      bennett_spins.push_back(delayed_spins[delayed_index]);

      bennett_spins.back().HS_current_move = ANNIHILATION;
    }

    if (delayed_spins[delayed_index].HS_current_move == CREATION)
      number_of_interacting_spins += 1;

    if (delayed_spins[delayed_index].HS_current_move == ANNIHILATION)
      number_of_interacting_spins -= 1;
  }
  else {
    /*
    Gamma_up_diag_max = tmp_up_diag_max<1. ? 1. : tmp_up_diag_max;
    Gamma_dn_diag_max = tmp_dn_diag_max<1. ? 1. : tmp_dn_diag_max;
    Gamma_up_diag_min = tmp_up_diag_min>1. ? 1. : tmp_up_diag_min;
    Gamma_dn_diag_min = tmp_dn_diag_min>1. ? 1. : tmp_dn_diag_min;

    delayed_spins[delayed_index].is_accepted_move  = false;

    if(delayed_spins[delayed_index].e_spin_HS_field_DN == e_UP)
      {
        CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_up_CPU, Gamma_up_size-1);
      }
    else
      {
        CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_dn_CPU, Gamma_dn_size-1);
      }

    if(delayed_spins[delayed_index].e_spin_HS_field_UP == e_UP)
      {
        CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_up_CPU, Gamma_up_size-1);
      }
    else
      {
        CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_dn_CPU, Gamma_dn_size-1);
      }
    */

    delayed_spins[delayed_index].is_accepted_move  = false;

    if (delayed_spins[delayed_index].e_spin_HS_field_DN == e_UP and
        delayed_spins[delayed_index].e_spin_HS_field_UP == e_UP) {
      Gamma_up_diag_max = tmp_up_diag_max < 1. ? 1. : tmp_up_diag_max;
      Gamma_up_diag_min = tmp_up_diag_min > 1. ? 1. : tmp_up_diag_min;

      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_up_CPU, Gamma_up_size - 2);
      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_up_CPU, Gamma_up_size - 1);
    }

    if (delayed_spins[delayed_index].e_spin_HS_field_DN == e_DN and
        delayed_spins[delayed_index].e_spin_HS_field_UP == e_UP) {
      Gamma_up_diag_max = tmp_up_diag_max < 1. ? 1. : tmp_up_diag_max;
      Gamma_dn_diag_max = tmp_dn_diag_max < 1. ? 1. : tmp_dn_diag_max;
      Gamma_up_diag_min = tmp_up_diag_min > 1. ? 1. : tmp_up_diag_min;
      Gamma_dn_diag_min = tmp_dn_diag_min > 1. ? 1. : tmp_dn_diag_min;

      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_dn_CPU, Gamma_dn_size - 1);
      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_up_CPU, Gamma_up_size - 1);
    }

    if (delayed_spins[delayed_index].e_spin_HS_field_DN == e_UP and
        delayed_spins[delayed_index].e_spin_HS_field_UP == e_DN) {
      Gamma_up_diag_max = tmp_up_diag_max < 1. ? 1. : tmp_up_diag_max;
      Gamma_dn_diag_max = tmp_dn_diag_max < 1. ? 1. : tmp_dn_diag_max;
      Gamma_up_diag_min = tmp_up_diag_min > 1. ? 1. : tmp_up_diag_min;
      Gamma_dn_diag_min = tmp_dn_diag_min > 1. ? 1. : tmp_dn_diag_min;

      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_up_CPU, Gamma_up_size - 1);
      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_dn_CPU, Gamma_dn_size - 1);
    }

    if (delayed_spins[delayed_index].e_spin_HS_field_DN == e_DN and
        delayed_spins[delayed_index].e_spin_HS_field_UP == e_DN) {
      Gamma_dn_diag_max = tmp_dn_diag_max < 1. ? 1. : tmp_dn_diag_max;
      Gamma_dn_diag_min = tmp_dn_diag_min > 1. ? 1. : tmp_dn_diag_min;

      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_dn_CPU, Gamma_dn_size - 2);
      CT_AUX_WALKER_TOOLS<dca::linalg::CPU, Scalar>::set_to_identity(Gamma_dn_CPU, Gamma_dn_size - 1);
    }
  }
}  // namespace ctaux

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::apply_bennett_on_Gamma_matrices(int& /*Gamma_up_size*/,
                                                                              int& /*Gamma_dn_size*/) {
  throw std::logic_error(__FUNCTION__);

  /*
    number_of_annihilations += 1;

    stored_Gamma_up_CPU = Gamma_up_CPU;
    stored_Gamma_dn_CPU = Gamma_dn_CPU;

    Scalar ratio_HS_field_DN = 0;
    Scalar ratio_HS_field_UP = 0;

    int bennett_index = int(rng()*Scalar(bennett_spins.size()));

    e_spin_states e_spin_HS_field_DN = bennett_spins[bennett_index].e_spin_HS_field_DN;
    e_spin_states e_spin_HS_field_UP = bennett_spins[bennett_index].e_spin_HS_field_UP;

    if(e_spin_HS_field_DN == e_UP)
    {
    int    index             = bennett_spins[bennett_index].Gamma_index_HS_field_DN;
    Scalar exp_minus_delta_V = bennett_spins[bennett_index].exp_minus_delta_V_HS_field_DN;

    ratio_HS_field_DN = ctaux_tools.apply_bennett_on_Gamma(index, Gamma_up_size, Gamma_up_CPU,
    exp_minus_delta_V);
    }
    else
    {
    int    index             = bennett_spins[bennett_index].Gamma_index_HS_field_DN;
    Scalar exp_minus_delta_V = bennett_spins[bennett_index].exp_minus_delta_V_HS_field_DN;

    ratio_HS_field_DN = ctaux_tools.apply_bennett_on_Gamma(index, Gamma_dn_size, Gamma_dn_CPU,
    exp_minus_delta_V);
    }

    if(e_spin_HS_field_UP == e_UP)
    {
    int    index             = bennett_spins[bennett_index].Gamma_index_HS_field_UP;
    Scalar exp_minus_delta_V = bennett_spins[bennett_index].exp_minus_delta_V_HS_field_UP;

    ratio_HS_field_UP = ctaux_tools.apply_bennett_on_Gamma(index, Gamma_up_size, Gamma_up_CPU,
    exp_minus_delta_V);
    }
    else
    {
    int    index             = bennett_spins[bennett_index].Gamma_index_HS_field_UP;
    Scalar exp_minus_delta_V = bennett_spins[bennett_index].exp_minus_delta_V_HS_field_UP;

    ratio_HS_field_UP = ctaux_tools.apply_bennett_on_Gamma(index, Gamma_dn_size, Gamma_dn_CPU,
    exp_minus_delta_V);
    }

    assert(bennett_spins[bennett_index].HS_current_move==ANNIHILATION);

    Scalar determinant_ratio = ratio_HS_field_UP*ratio_HS_field_DN;
    Scalar acceptance_ratio  = calculate_acceptace_ratio(determinant_ratio,
    bennett_spins[bennett_index].HS_current_move);

    if( std::abs(acceptance_ratio) >= rng() )
    {
    number_of_interacting_spins -= 1;

    if(acceptance_ratio < 0)
    phase_ *= -1;

    {// set column and row to zero
    if(bennett_spins[bennett_index].e_spin_HS_field_DN == e_UP)
    {
    int k = bennett_spins[bennett_index].Gamma_index_HS_field_DN;
    CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_up_CPU, k);
    }
    else
    {
    int k = bennett_spins[bennett_index].Gamma_index_HS_field_DN;
    CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_dn_CPU, k);
    }

    if(bennett_spins[bennett_index].e_spin_HS_field_UP == e_UP)
    {
    int k = bennett_spins[bennett_index].Gamma_index_HS_field_UP;
    CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_up_CPU, k);
    }
    else
    {
    int k = bennett_spins[bennett_index].Gamma_index_HS_field_UP;
    CT_AUX_WALKER_TOOLS<device_t, Scalar>::set_to_identity(Gamma_dn_CPU, k);
    }
    }

    {// update the delayed spins and erase the bennett-spin
    int delayed_spin_index = bennett_spins[bennett_index].delayed_spin_index;

    assert(not delayed_spins[delayed_spin_index].is_a_bennett_spin);
    delayed_spins[delayed_spin_index].is_a_bennett_spin = true;

    bennett_spins.erase(bennett_spins.begin()+bennett_index, bennett_spins.begin()+bennett_index+1);
    }
    }
    else
    {
    Gamma_up_CPU.swap(stored_Gamma_up_CPU);
    Gamma_dn_CPU.swap(stored_Gamma_dn_CPU);
    }
  */
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::update_N_matrix_with_Gamma_matrix() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  // kills Bennett-spins and puts the interacting vertices all in the left part of the configuration_
  SHRINK_TOOLS<Profiler, device_t, Scalar>::shrink_Gamma(configuration_, Gamma_up, Gamma_dn);

  N_tools_obj.rebuild_N_matrix_via_Gamma_LU(configuration_, N_up, Gamma_up, G_up, e_UP);
  N_tools_obj.rebuild_N_matrix_via_Gamma_LU(configuration_, N_dn, Gamma_dn, G_dn, e_DN);

  configuration_.commit_accepted_spins();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::clean_up_the_configuration() {
  Profiler profiler(__FUNCTION__, "CT-AUX walker", __LINE__, thread_id);

  SHRINK_tools_obj.reorganize_configuration_test(configuration_, N_up, N_dn, G0_up, G0_dn);

  assert(configuration_.assert_consistency());
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
HS_vertex_move_type CtauxWalker<device_t, Parameters, Data>::get_new_HS_move() {
  if (rng() > 0.5) {
    return CREATION;
  }

  else if (configuration_.get_number_of_interacting_HS_spins() == 0)
    return STATIC;

  else
    return ANNIHILATION;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
int CtauxWalker<device_t, Parameters, Data>::get_new_vertex_index(HS_vertex_move_type HS_current_move) {
  if (HS_current_move == CREATION)
    return configuration_.get_random_noninteracting_vertex();

  return configuration_.get_random_interacting_vertex();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
HS_spin_states_type CtauxWalker<device_t, Parameters, Data>::get_new_spin_value(
    HS_vertex_move_type HS_current_move) {
  if (HS_current_move == CREATION) {
    if (rng() > 0.5)
      return HS_UP;
    else
      return HS_DN;
  }

  return HS_ZERO;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
auto CtauxWalker<device_t, Parameters, Data>::calculate_acceptance_ratio(
    Scalar determinant_ratio, HS_vertex_move_type HS_current_move, Real QMC_factor) const {
  Real N = number_of_interacting_spins;
  Real K = parameters_.get_expansion_parameter_K();

  if (HS_current_move == CREATION) {
    return K / (N + 1.) * determinant_ratio / QMC_factor;
  }
  else {
    return (N / K) * determinant_ratio / QMC_factor;
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
bool CtauxWalker<device_t, Parameters, Data>::assert_exp_delta_V_value(
    HS_field_sign HS_field, int random_vertex_ind, HS_spin_states_type new_HS_spin_value,
    Scalar exp_delta_V) {
  switch (HS_field) {
    case HS_FIELD_DN: {
      e_spin_states e_spin_HS_field_DN = configuration_[random_vertex_ind].get_e_spins().first;
      int configuration_e_spin_index_HS_field_DN =
          configuration_[random_vertex_ind].get_configuration_e_spin_indices().first;

      vertex_singleton& v_j =
          configuration_.get(e_spin_HS_field_DN)[configuration_e_spin_index_HS_field_DN];

      if (std::abs(CV_obj.exp_delta_V(v_j, new_HS_spin_value) - exp_delta_V) > 1.e-6) {
        std::cout << HS_field << "\t" << e_spin_HS_field_DN << std::endl;
        throw std::logic_error(__FUNCTION__);
      }
    } break;

    case HS_FIELD_UP: {
      e_spin_states e_spin_HS_field_UP = configuration_[random_vertex_ind].get_e_spins().second;
      int configuration_e_spin_index_HS_field_UP =
          configuration_[random_vertex_ind].get_configuration_e_spin_indices().second;

      vertex_singleton& v_j =
          configuration_.get(e_spin_HS_field_UP)[configuration_e_spin_index_HS_field_UP];

      if (std::abs(CV_obj.exp_delta_V(v_j, new_HS_spin_value) - exp_delta_V) > 1.e-6) {
        std::cout << HS_field << "\t" << e_spin_HS_field_UP << std::endl;
        throw std::logic_error(__FUNCTION__);
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  };

  return true;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::updateShell(const int done, const int total) {
  const int milestone = std::max(total / 10, 1);
  if (concurrency_.id() == concurrency_.first() && (done % milestone == 0)) {
    std::cout.unsetf(std::ios_base::floatfield);

    std::cout << "\t\t\t" << std::setw(14)
              << static_cast<Real>(done) / static_cast<Real>(total) * 100. << " % completed"
              << "\t" << std::setw(11)
              << "<k> = " << configuration_.get_number_of_interacting_HS_spins() << "\t"
              << std::setw(11) << "N = " << configuration_.size() << "\t" << dca::util::print_time()
              << std::endl;

    std::cout << std::scientific;
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <typename AccumType>
const linalg::util::GpuEvent* CtauxWalker<device_t, Parameters, Data>::computeM(
  std::array<linalg::Matrix<AccumType, device_t>, 2>& Ms) {
  linalg::util::syncStream(thread_id, stream_id);

  for (int s = 0; s < 2; ++s) {
    const auto& config = get_configuration().get(s == 0 ? e_DN : e_UP);
    exp_v_minus_one_[s].resizeNoCopy(config.size());

    for (int i = 0; i < config.size(); ++i)
      exp_v_minus_one_[s][i] = CV_obj.exp_V(config[i]) - 1.;

    const auto& N = s == 0 ? N_dn : N_up;
    auto& M = Ms[s];

    M.resizeNoCopy(N.size());

    if (device_t == linalg::GPU) {
      exp_v_minus_one_dev_[s].setAsync(exp_v_minus_one_[s], thread_id, stream_id);
      dca::linalg::matrixop::multiplyDiagonalLeft(exp_v_minus_one_dev_[s], N, M, thread_id, stream_id);
    }
    else {
      dca::linalg::matrixop::multiplyDiagonalLeft(exp_v_minus_one_[s], N, M);
    }
  }
  m_computed_events_[0].record(linalg::util::getStream(thread_id, stream_id));
  return &m_computed_events_[0];
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::readConfig(dca::io::Buffer& buff) {
  buff >> configuration_;
  config_initialized_ = true;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
io::Buffer CtauxWalker<device_t, Parameters, Data>::dumpConfig() const {
  io::Buffer buff;
  buff << configuration_;
  return buff;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxWalker<device_t, Parameters, Data>::recomputeMCWeight() {
  mc_log_weight_ = 0.;
  phase_.reset();

  auto process_matrix = [&](auto& m) {
    if (!m.nrRows())
      return;

    const auto [log_det, phase] = linalg::matrixop::logDeterminant(m);
    mc_log_weight_ -= log_det;  // MC weight is proportional to det(N^-1)
    phase_.divide(phase);
  };

  if constexpr (decltype(N_up)::device == dca::linalg::DeviceType::GPU) {
    linalg::Matrix<typename decltype(N_up)::ValueType, dca::linalg::DeviceType::CPU> N_up_CPU{N_up};
    linalg::Matrix<typename decltype(N_dn)::ValueType, dca::linalg::DeviceType::CPU> N_dn_CPU{N_dn};
    process_matrix(N_up_CPU);
    process_matrix(N_dn_CPU);
  }
  else {
    process_matrix(N_up);
    process_matrix(N_dn);
  }

  mc_log_weight_ += mc_log_weight_constant_ * configuration_.size();
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_WALKER_HPP
