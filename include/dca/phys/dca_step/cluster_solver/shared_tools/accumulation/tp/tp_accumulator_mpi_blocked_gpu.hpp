// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Weile Wei (wwei9@lsu.edu)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// Implementation of the two particle Green's function computation on the GPU with distrubtion
// over MPI.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_MPI_BLOCKED_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_MPI_BLOCKED_GPU_HPP

#include "dca/config/dca.hpp"
#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA."
#endif

#ifndef DCA_HAVE_MPI
#error "This file requires MPI."
#endif

#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"
#include "dca/parallel/util/call_once_per_loop.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_base.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu_base.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <class Parameters>
class TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>
    : public TpAccumulatorBase<Parameters, DistType::BLOCKED>,
      public TpAccumulatorGpuBase<Parameters, DistType::BLOCKED> {
public:
  static constexpr DistType DT = DistType::BLOCKED;
  using Base = TpAccumulatorBase<Parameters, DT>;
  using BaseGpu = TpAccumulatorGpuBase<Parameters, DT>;
  using ThisType = TpAccumulator<Parameters, linalg::GPU, DT>;

private:
  // is there a smarter way to do this in c++17?
  using RDmn = typename Base::RDmn;
  using KDmn = typename Base::KDmn;
  using NuDmn = typename Base::NuDmn;
  using WDmn = typename Base::WDmn;

  using typename Base::Real;
  using typename Base::Complex;

  using Base::channels_;
  using Base::G4_;

  using Base::n_bands_;
  using Base::sign_;
  using Base::multiple_accumulators_;
  using Base::thread_id_;
  using typename Base::WTpDmn;
  using typename Base::WTpExtDmn;
  using typename Base::WTpExtPosDmn;
  using typename Base::WTpPosDmn;
  using typename BaseGpu::RMatrix;
  using typename BaseGpu::RMatrixValueType;
  using Base::beta_;

  using BaseGpu::queues_;
  using BaseGpu::ndft_objs_;
  using BaseGpu::workspaces_;
  using BaseGpu::G_;
  using BaseGpu::space_trsf_objs_;
  using BaseGpu::nr_accumulators_;
  using BaseGpu::event_;
  using BaseGpu::n_pos_frqs_;
  using BaseGpu::n_ndft_queues_;
  uint64_t start_;
  uint64_t end_;

  // Eventually distribution strategy should be pushed down into linalg::Vector but
  // I think generalization should still wait.
  using G4DevType = linalg::Vector<Complex, linalg::GPU, config::McOptions::TpAllocator<Complex>>;

  using BaseGpu::get_G0;

public:
  template <class Configuration, typename RealIn>
  float accumulate(const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M,
                   const std::array<Configuration, 2>& configs, int sign);

  // CPU input. For testing purposes.
  template <class Configuration>
  float accumulate(const std::array<linalg::Matrix<double, linalg::CPU>, 2>& M,
                   const std::array<Configuration, 2>& configs, int sign);

  // Downloads the accumulation result to the host.
  void finalize();

  TpAccumulator(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int thread_id = 0);

  // \todo resets are almost always a sign of a design problem, usually an object
  //  that should be local to the body of a loop being persisted longer.
  //  I think this is true here.
  //  1. there is an expensive allocation?
  //  2. Is interation to iteration state is being stashed in this object?
  //  3. some other reason

  // Resets the object between DCA iterations.
  void resetAccumulation(unsigned int dca_loop);

  // Returns the accumulated Green's function.
  std::vector<typename Base::TpGreensFunction>& get_G4();

  // Sums the accumulated Green's function to the accumulated Green's function of other_acc.
  void sumTo(TpAccumulator& other_acc);

private:
  static inline std::vector<G4DevType>& get_G4Dev();

  void computeGMultiband(int s);

  void computeGSingleband(int s);

  void computeG();

  template <class Configuration, typename RealIn>
  float computeM(const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M_pair,
                 const std::array<Configuration, 2>& configs);

  // The semantics of this class used to make it possible that the G4 would get resized here if
  // Base::TpDomain had changed this is no longer true as the accessible size of G4 is set in the constructor.
  void resetG4();
  // Applies pipepline ring algorithm to move G matrices around all ranks
  void ringG(float& flop);
  float updateG4(const std::size_t channel_index);

  void send(const std::array<RMatrix, 2>& data, int target, std::array<MPI_Request, 2>& request);
  void receive(std::array<RMatrix, 2>& data, int source, std::array<MPI_Request, 2>& request);

  // send buffer for pipeline ring algorithm
  std::array<RMatrix, 2> sendbuff_G_;

  std::array<MPI_Request, 2> recv_requests_{MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  std::array<MPI_Request, 2> send_requests_{MPI_REQUEST_NULL, MPI_REQUEST_NULL};

#ifndef DCA_WITH_CUDA_AWARE_MPI
  std::array<std::vector<Complex>, 2> sendbuffer_;
  std::array<std::vector<Complex>, 2> recvbuffer_;
#endif  // DCA_WITH_CUDA_AWARE_MPI
};

template <class Parameters>
TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::TpAccumulator(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int thread_id)
    : Base(G0, pars, thread_id), BaseGpu(pars, Base::get_n_pos_frqs(), thread_id_) {
  // each mpi rank only allocates memory of size 1/total_G4_size for its small portion of G4
  // static_assert(std::is_same<std::vector<TpAccumulator<DistType::BLOCKED>>, decltype(G4_)>);
  start_ = G4_[0].get_start();
  // The sense here is one past the last index held
  end_ = G4_[0].get_end() + 1;

  // possible these can both go into the parent class constructor
}

template <class Parameters>
template <class Configuration, typename RealIn>
float TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::accumulate(
    const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M,
    const std::array<Configuration, 2>& configs, const int sign) {
  // typename Base::Profiler profiler("accumulate", "tp-accumulation", __LINE__, Base::thread_id_);
  float flop = 0;

  if (!BaseGpu::initialized_)
    throw(std::logic_error("The accumulator is not ready to measure."));

  if (!(configs[0].size() + configs[0].size()))  // empty config
    return flop;

  Base::sign_ = sign;
  flop += BaseGpu::computeM(M, configs);
  computeG();

  for (std::size_t channel = 0; channel < G4_.size(); ++channel) {
    flop += updateG4(channel);
  }

  ringG(flop);

  return flop;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::computeGSingleband(const int s) {
  details::computeGSingleband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                              KDmn::dmn_size(), n_pos_frqs_, beta_, queues_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::computeGMultiband(const int s) {
  details::computeGMultiband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                             get_G0()[s].leadingDimension(), n_bands_, KDmn::dmn_size(),
                             n_pos_frqs_, beta_, queues_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::computeG() {
  if (n_ndft_queues_ == 1) {
    event_.record(queues_[0]);
    event_.block(queues_[1]);
  }
  {
    [[maybe_unused]] Profiler prf("ComputeG: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s) {
      if (n_bands_ == 1)
        computeGSingleband(s);
      else
        computeGMultiband(s);
    }
  }

  event_.record(queues_[1]);
  event_.block(queues_[0]);
}

template <class Parameters>
template <class Configuration, typename RealIn>
float TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::computeM(
    const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs) {
  auto stream_id = [&](const int s) { return n_ndft_queues_ == 1 ? 0 : s; };

  float flop = 0.;

  {
    [[maybe_unused]] Profiler prf("Frequency FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      flop += ndft_objs_[stream_id(s)].execute(configs[s], M_pair[s], G_[s]);
  }
  {
    [[maybe_unused]] Profiler prf("Space FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      flop += space_trsf_objs_[stream_id(s)].execute(G_[s]);
  }

  return flop;
}

template <class Parameters>
template <class Configuration>
float TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::accumulate(
    const std::array<linalg::Matrix<double, linalg::CPU>, 2>& M,
    const std::array<Configuration, 2>& configs, const int sign) {
  std::array<linalg::Matrix<double, linalg::GPU>, 2> M_dev;
  for (int s = 0; s < 2; ++s)
    M_dev[s].setAsync(M[s], queues_[0]);

  return accumulate(M_dev, configs, sign);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::finalize() {
  if (BaseGpu::finalized_)
    return;

  for (std::size_t channel = 0; channel < G4_.size(); ++channel) {
    // modify G4 size in G4 cpu, otherwise, copyTo() operation failed due to incomparable size
    // resize() only modifies member Nb_elements in function, does not change tp_dmn.get_size()
    G4_[channel].resize(get_G4()[channel].size());
    get_G4Dev()[channel].copyTo(G4_[channel]);
  }
  // TODO: release memory if needed by the rest of the DCA loop.
  // get_G4().clear();

  BaseGpu::finalized_ = true;
  BaseGpu::initialized_ = false;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::resetAccumulation(
    const unsigned int dca_loop) {
  static dca::util::OncePerLoopFlag flag;

  dca::util::callOncePerLoop(flag, dca_loop, [&]() {
    resetG4();
    BaseGpu::initializeG0();
    BaseGpu::synchronizeStreams();
  });

  BaseGpu::initialized_ = true;
  BaseGpu::finalized_ = false;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::resetG4() {
  // Note: this method is not thread safe by itself.
  get_G4Dev().resize(G4_.size());

  // These are the device G4's
  for (auto& G4_channel : get_G4Dev()) {
    try {
      if (!Base::multiple_accumulators_) {
        G4_channel.setStream(BaseGpu::queues_[0]);
      }

      G4_channel.resizeNoCopy(end_ - start_);
      G4_channel.setToZeroAsync(BaseGpu::queues_[0]);
    }
    catch (std::bad_alloc& err) {
      std::cerr << "Failed to allocate G4 on device.\n";
      throw(err);
    }
  }
}

template <class Parameters>
float TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::updateG4(
    const std::size_t channel_index) {
  // G4 is stored with the following band convention:
  // b1 ------------------------ b3
  //        |           |
  //        |           |
  //        |           |
  // b2 ------------------------ b4

  //  TODO: set stream only if this thread gets exclusive access to G4.
  //  get_G4().setStream(queues_[0]);

  const FourPointType channel = Base::channels_[channel_index];

  int my_rank, mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  switch (channel) {
    case FourPointType::PARTICLE_HOLE_TRANSVERSE:
      return details::updateG4<Real, FourPointType::PARTICLE_HOLE_TRANSVERSE>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start_, end_);

    case FourPointType::PARTICLE_HOLE_MAGNETIC:
      return details::updateG4<Real, FourPointType::PARTICLE_HOLE_MAGNETIC>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start_, end_);

    case FourPointType::PARTICLE_HOLE_CHARGE:
      return details::updateG4<Real, FourPointType::PARTICLE_HOLE_CHARGE>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start_, end_);

    case FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP:
      return details::updateG4<Real, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start_, end_);

    case FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN:
      return details::updateG4<Real, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start_, end_);

    case FourPointType::PARTICLE_PARTICLE_UP_DOWN:
      return details::updateG4<Real, FourPointType::PARTICLE_PARTICLE_UP_DOWN>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start_, end_);

    default:
      throw std::logic_error("Specified four point type not implemented.");
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::ringG(float& flop) {
  // get ready for send and receive

  for (int s = 0; s < 2; ++s) {
    sendbuff_G_[s] = G_[s];
  }

  int my_concurrency_id, mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_concurrency_id);

  // get rank index of left and right neighbor
  int left_neighbor = (my_concurrency_id - 1 + mpi_size) % mpi_size;
  int right_neighbor = (my_concurrency_id + 1 + mpi_size) % mpi_size;

  // Pipepline ring algorithm in the following for-loop:
  // 1) At each time step, local rank receives a new G2 from left hand neighbor,
  // makes a copy locally and uses this G2 to update G4, and
  // sends this G2 to right hand neighbor. In total, the algorithm performs (mpi_size - 1) steps.
  // 2) This algorithm currently requires parameters in input file, please refer to mci_parameters.hpp:
  //      a) #walker == #accumulator and shared-walk-and-accumulation-thread = true;
  //      b) and, local measurements are equal, and each accumulator should have same #measurement, i.e.
  //         measurements % ranks == 0 && local_measurement % threads == 0.
  for (int icount = 0; icount < (mpi_size - 1); icount++) {
    send(sendbuff_G_, right_neighbor, send_requests_);
    receive(G_, left_neighbor, recv_requests_);

    // wait for G2 to be available again
    for (int s = 0; s < 2; ++s)
      MPI_Wait(&recv_requests_[s], MPI_STATUSES_IGNORE);

    // use newly copied G2 to update G4
    for (std::size_t channel = 0; channel < G4_.size(); ++channel) {
      flop += updateG4(channel);
    }

    // wait for sendbuf_G2 to be available again
    for (int s = 0; s < 2; ++s)
      MPI_Wait(&send_requests_[s], MPI_STATUSES_IGNORE);

    // get ready for send again
    for (int s = 0; s < 2; ++s) {
      sendbuff_G_[s].swap(G_[s]);
    }
  }
}

template <class Parameters>
auto TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::get_G4Dev() -> std::vector<G4DevType>& {
  static std::vector<G4DevType> G4;
  return G4;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::send(
    const std::array<RMatrix, 2>& data, int target, std::array<MPI_Request, 2>& request) {
  using dca::parallel::MPITypeMap;
  const auto g_size = data[0].size().first * data[0].size().second;

#ifdef DCA_WITH_CUDA_AWARE_MPI
  for (int s = 0; s < 2; ++s) {
    MPI_Isend(data[s].ptr(), g_size, MPITypeMap<Complex>::value(), target, thread_id_ + 1,
              MPI_COMM_WORLD, &request[s]);
  }
#else

  for (int s = 0; s < 2; ++s) {
    sendbuffer_[s].resize(g_size);
    cudaMemcpy(sendbuffer_[s].data(), data[s].ptr(), g_size * sizeof(Complex),
               cudaMemcpyDeviceToHost);

    MPI_Isend(sendbuffer_[s].data(), g_size, MPITypeMap<Complex>::value(), target, thread_id_ + 1,
              MPI_COMM_WORLD, &request[s]);
  }
#endif  // DCA_WITH_CUDA_AWARE_MPI
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::receive(
    std::array<RMatrix, 2>& data, int source, std::array<MPI_Request, 2>& request) {
  using dca::parallel::MPITypeMap;
  const auto g_size = data[0].size().first * data[0].size().second;

#ifdef DCA_WITH_CUDA_AWARE_MPI
  for (int s = 0; s < 2; ++s) {
    MPI_Irecv(data[s].ptr(), g_size, MPITypeMap<Complex>::value(), source, thread_id_ + 1,
              MPI_COMM_WORLD, &request[s]);
  }

#else
  for (int s = 0; s < 2; ++s) {
    recvbuffer_[s].resize(g_size);
    MPI_Irecv(recvbuffer_[s].data(), g_size, MPITypeMap<Complex>::value(), source, thread_id_ + 1,
              MPI_COMM_WORLD, &request[s]);
  }
  for (int s = 0; s < 2; ++s) {
    MPI_Wait(&request[s], MPI_STATUSES_IGNORE);
    // Note: MPI can not access host memory allocated from CUDA, hence the usage of blocking communication.
    cudaMemcpy(data[s].ptr(), recvbuffer_[s].data(), g_size * sizeof(Complex),
               cudaMemcpyHostToDevice);
  }
#endif  // DCA_WITH_CUDA_AWARE_MPI
}

/** Return the G4 only as the correct type
 *
 *  the return type is quite a code smell
 */
template <class Parameters>
std::vector<typename TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::Base::TpGreensFunction>& TpAccumulator<
    Parameters, linalg::GPU, DistType::BLOCKED>::get_G4() {
  if (G4_.empty())
    throw std::logic_error("There is no G4 stored in this class.");

  return G4_;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU, DistType::BLOCKED>::sumTo(TpAccumulator& other_one) {
  BaseGpu::sumTo_(other_one);
}


}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca
#endif
