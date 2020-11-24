// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// Implementation of the two particle Green's function computation on the GPU.

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA."
#endif

#ifndef DCA_TP_ACCUMULATOR_GPU_HPP
#define DCA_TP_ACCUMULATOR_GPU_HPP

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_base.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu_base.hpp"

#include <cuda.h>
#include <mutex>
#include <vector>

#include "dca/config/mc_options.hpp"
#include "dca/config/profiler.hpp"
#include "dca/distribution/dist_types.hpp"
#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/reshapable_matrix.hpp"
#include "dca/linalg/util/allocators/managed_allocator.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_queue.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D_gpu.hpp"
#include "dca/parallel/util/call_once_per_loop.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/kernels_interface.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_gpu.hpp"

#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <class Parameters, DistType DT>
class TpAccumulator<Parameters, DT, linalg::GPU>
  : public TpAccumulatorBase<Parameters, DT>, public TpAccumulatorGpuBase<Parameters, DT>
{
public:
  static constexpr DistType dist = DT;
  using Base = TpAccumulatorBase<Parameters, DT>;
  using BaseGpu = TpAccumulatorGpuBase<Parameters, DT>;
  using ThisType = TpAccumulator<Parameters, DT, linalg::GPU>;
  using typename Base::Real;
  using RDmn = typename Base::RDmn;
  using KDmn = typename Base::KDmn;
  using NuDmn = typename Base::NuDmn;
  using WDmn = typename Base::WDmn;
  using typename Base::WTpDmn;
  using typename Base::WTpExtDmn;
  using typename Base::WTpExtPosDmn;
  using typename Base::WTpPosDmn;
  using typename Base::BDmn;
  using typename Base::SDmn;
  using typename Base::Complex;
  using typename Base::TpGreensFunction;
protected:
  using BaseGpu::queues_;
  using BaseGpu::ndft_objs_;
  using BaseGpu::workspaces_;
  using BaseGpu::G_;
  using BaseGpu::space_trsf_objs_;
  using BaseGpu::nr_accumulators_;
  using BaseGpu::event_;
  using typename BaseGpu::RMatrix;
  using typename BaseGpu::RMatrixValueType;

public:
  // Constructor:
  // In: G0: non interacting greens function.
  // In: pars: parameters object.
  // In: thread_id: thread id, only used by the profiler.
  TpAccumulator(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int thread_id = 0);

  // Resets the object between DCA iterations.
  void resetAccumulation(unsigned int dca_loop);

  // Computes the two particles Greens function from the M matrix and accumulates it internally.
  // In: M_array: stores the M matrix for each spin sector.
  // In: configs: stores the walker's configuration for each spin sector.
  // In: sign: sign of the configuration.
  // Returns: number of flop.
  template <class Configuration, typename RealIn>
  float accumulate(const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M,
                   const std::array<Configuration, 2>& configs, int sign);

  // CPU input. For testing purposes.
  template <class Configuration>
  float accumulate(const std::array<linalg::Matrix<double, linalg::CPU>, 2>& M,
                   const std::array<Configuration, 2>& configs, int sign);

  // Downloads the accumulation result to the host.
  void finalize();

  // Sums on the host the accumulated Green's function to the accumulated Green's function of
  // other_acc.
  void sumTo(ThisType& other_acc);

  const linalg::util::CudaStream* get_stream() {
    return &queues_[0].getStream();
  }

  void synchronizeCopy() {
    ndft_objs_[0].synchronizeCopy();
    ndft_objs_[1].synchronizeCopy();
  }

  void syncStreams(const linalg::util::CudaEvent& event) {
    for (const auto& stream : queues_)
      event.block(stream);
  }

  std::size_t deviceFingerprint() const {
    std::size_t res(0);
    for (const auto& work : workspaces_)
      res += work->deviceFingerprint();
    for (int s = 0; s < 2; ++s)
      res += ndft_objs_[s].deviceFingerprint() + G_[s].deviceFingerprint() +
             space_trsf_objs_[s].deviceFingerprint();
    return res;
  }

  static std::size_t staticDeviceFingerprint() {
    std::size_t res = 0;

    for (const auto& G4_channel : get_G4Dev())
      res += G4_channel.deviceFingerprint();

    res += get_G0()[0].deviceFingerprint() + get_G0()[1].deviceFingerprint();

    return res;
  }

      // Returns the accumulated Green's function.
  std::vector<typename TpAccumulator<Parameters, DT, linalg::GPU>::Base::TpGreensFunction>& get_G4();

protected:
  using typename BaseGpu::Matrix;
  using typename BaseGpu::MatrixHost;
  using Base::beta_;
  using Base::channels_;
  using Base::extension_index_offset_;
  using Base::G4_;
  using Base::multiple_accumulators_;
  using Base::n_bands_;
  using Base::n_pos_frqs_;

  using Base::non_density_density_;
  using Base::sign_;
  using Base::thread_id_;

  using BaseGpu::initialized_;
  using BaseGpu::finalized_;

  using BaseGpu::n_ndft_queues_;
  
  using G4DevType = linalg::Vector<Complex, linalg::GPU, config::McOptions::TpAllocator<Complex>>;

  using BaseGpu::get_G0;

  void initializeG0();
  void resetG4();

  void computeG();

  void computeGMultiband(int s);

  void computeGSingleband(int s);

  float updateG4(const std::size_t channel_index);

  void synchronizeStreams();

  static inline std::vector<G4DevType>& get_G4Dev();

  // For distributed G4's
  // Applies pipepline ring algorithm to move G matrices around all ranks
  void ringG(float& flop);

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

template <class Parameters, DistType DT>
TpAccumulator<Parameters, DT, linalg::GPU>::TpAccumulator(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int thread_id)
    : Base(G0, pars, thread_id),
      BaseGpu(pars, Base::get_n_pos_frqs(), thread_id) {}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::resetAccumulation(const unsigned int dca_loop) {
  static dca::util::OncePerLoopFlag flag;

  dca::util::callOncePerLoop(flag, dca_loop, [&]() {
    resetG4();
    BaseGpu::initializeG0();
    BaseGpu::synchronizeStreams();
  });

  initialized_ = true;
  finalized_ = false;
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::resetG4() {
  // Note: this method is not thread safe by itself.
  get_G4Dev().resize(G4_.size());

  for (auto& G4_channel : get_G4Dev()) {
    try {
      typename Base::TpDomain tp_dmn;
      if (!multiple_accumulators_) {
        G4_channel.setStream(queues_[0]);
      }
      G4_channel.resizeNoCopy(tp_dmn.get_size());
      G4_channel.setToZeroAsync(queues_[0]);
    }
    catch (std::bad_alloc& err) {
      std::cerr << "Failed to allocate G4 on device.\n";
      throw(err);
    }
  }
}

template <class Parameters, DistType DT>
template <class Configuration, typename RealIn>
float TpAccumulator<Parameters, DT, linalg::GPU>::accumulate(
    const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M,
    const std::array<Configuration, 2>& configs, const int sign) {
  [[maybe_unused]] Profiler profiler("accumulate", "tp-accumulation", __LINE__, thread_id_);
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

  if constexpr (dist == DistType::LINEAR || dist == DistType::BLOCKED) {
      ringG(flop);
    }
  return flop;
}

template <class Parameters, DistType DT>
template <class Configuration>
float TpAccumulator<Parameters, DT, linalg::GPU>::accumulate(
    const std::array<linalg::Matrix<double, linalg::CPU>, 2>& M,
    const std::array<Configuration, 2>& configs, const int sign) {
  std::array<linalg::Matrix<double, linalg::GPU>, 2> M_dev;
  for (int s = 0; s < 2; ++s)
    M_dev[s].setAsync(M[s], queues_[0]);

  return accumulate(M_dev, configs, sign);
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::computeG() {
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

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::computeGSingleband(const int s) {
  details::computeGSingleband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                              KDmn::dmn_size(), n_pos_frqs_, beta_, queues_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::computeGMultiband(const int s) {
  details::computeGMultiband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                             get_G0()[s].leadingDimension(), n_bands_, KDmn::dmn_size(),
                             n_pos_frqs_, beta_, queues_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters, DistType DT>
float TpAccumulator<Parameters, DT, linalg::GPU>::updateG4(const std::size_t channel_index) {
  // G4 is stored with the following band convention:
  // b1 ------------------------ b3
  //        |           |
  //        |           |
  //        |           |
  // b2 ------------------------ b4

  //  TODO: set stream only if this thread gets exclusive access to G4.
  //  get_G4().setStream(queues_[0]);
  const FourPointType channel = Base::channels_[channel_index];
  typename Base::TpDomain tp_dmn;
  uint64_t start = Base::G4_[0].get_start();
  uint64_t end = Base::G4_[0].get_end();
  switch (channel) {
    case PARTICLE_HOLE_TRANSVERSE:
      return details::updateG4<Real, PARTICLE_HOLE_TRANSVERSE>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_MAGNETIC:
      return details::updateG4<Real, PARTICLE_HOLE_MAGNETIC>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_CHARGE:
      return details::updateG4<Real, PARTICLE_HOLE_CHARGE>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_LONGITUDINAL_UP_UP:
      return details::updateG4<Real, PARTICLE_HOLE_LONGITUDINAL_UP_UP>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_LONGITUDINAL_UP_DOWN:
      return details::updateG4<Real, PARTICLE_HOLE_LONGITUDINAL_UP_DOWN>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_PARTICLE_UP_DOWN:
      return details::updateG4<Real, PARTICLE_PARTICLE_UP_DOWN>(
          get_G4Dev()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    default:
      throw std::logic_error("Specified four point type not implemented.");
  }
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::finalize() {
  if (finalized_)
    return;

  for (std::size_t channel = 0; channel < G4_.size(); ++channel) {
    get_G4Dev()[channel].copyTo(G4_[channel]);
  }
  // TODO: release memory if needed by the rest of the DCA loop.
  // get_G4().clear();

  finalized_ = true;
  initialized_ = false;
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::sumTo(ThisType& other/*other_one*/) {
  // Nothing to do: G4 on the device is shared.
  BaseGpu::sumTo_(other);
  return;
}

template <class Parameters, DistType DT>
auto TpAccumulator<Parameters, DT, linalg::GPU>::get_G4Dev() -> std::vector<G4DevType>& {
  static std::vector<G4DevType> G4;
  return G4;
}

/** Return the G4 only as the correct type
 *
 *  the return type is quite a code smell
 */
template <class Parameters, DistType DT>
std::vector<typename TpAccumulator<Parameters, DT, linalg::GPU>::Base::TpGreensFunction>& TpAccumulator<
  Parameters, DT, linalg::GPU>::get_G4() {
    
  if (G4_.empty())
    throw std::logic_error("There is no G4 stored in this class.");

  return G4_;
}


template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::ringG(float& flop) {
  // get ready for send and receive
  static_assert(dist != DistType::NONE);
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

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::send(
    const std::array<RMatrix, 2>& data, int target, std::array<MPI_Request, 2>& request) {
  static_assert(dist != DistType::NONE);
  using dca::parallel::MPITypeMap;
  const auto g_size = data[0].size().first * data[0].size().second;

#ifdef DCA_WITH_CUDA_AWARE_MPI
  for (int s = 0; s < 2; ++s) {
    MPI_Isend(data[s].ptr(), g_size, MPITypeMap<Complex>::value(), target, thread_id_ + 1,
              MPI_COMM_WORLD, &requestxo[s]);
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

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::GPU>::receive(
    std::array<RMatrix, 2>& data, int source, std::array<MPI_Request, 2>& request) {
  static_assert(dist != DistType::NONE);
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


}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACUCMULATOR_GPU_HPP
