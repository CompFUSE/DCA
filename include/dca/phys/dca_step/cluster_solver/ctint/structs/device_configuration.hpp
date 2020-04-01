// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This class defines a trivially copiable object to hold the CT-INT configuration on the device.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_DEVICE_CONFIGURATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_DEVICE_CONFIGURATION_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/util/cuda_definitions.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/sector_entry.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

struct DeviceConfiguration {
  const SectorEntry* data;

  __DEVICE__ inline unsigned short getLeftB(const int matrix_index) const;

  __DEVICE__ inline unsigned short getRightB(const int matrix_index) const;

  __DEVICE__ inline unsigned short getLeftR(const int matrix_index) const;

  __DEVICE__ inline unsigned short getRightR(const int matrix_index) const;

  __DEVICE__ inline double getTau(const int matrix_index) const;

  __DEVICE__ inline short getAuxFieldType(int matrix_index) const;
};

__DEVICE__
unsigned short DeviceConfiguration::getLeftB(const int matrix_index) const {
  return data[matrix_index].b_left_;
}

__DEVICE__
unsigned short DeviceConfiguration::getRightB(const int matrix_index) const {
  return data[matrix_index].b_right_;
}

__DEVICE__
unsigned short DeviceConfiguration::getLeftR(const int matrix_index) const {
  return data[matrix_index].r_left_;
}

__DEVICE__
unsigned short DeviceConfiguration::getRightR(const int matrix_index) const {
  return data[matrix_index].r_right_;
}

__DEVICE__
double DeviceConfiguration::getTau(const int matrix_index) const {
  return data[matrix_index].tau_;
}

__DEVICE__
short DeviceConfiguration::getAuxFieldType(int matrix_index) const {
  return data[matrix_index].aux_field_type_;
}

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_DEVICE_CONFIGURATION_HPP
