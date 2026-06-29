// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// This file tests the analytic three-band Hubbard Hamiltonian and its model parameters.

#include "dca/platform/dca_gpu.h"
#include "dca/phys/models/analytic_hamiltonians/threeband_hubbard.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/model_parameters.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"

namespace {

using PointGroup = dca::phys::domains::D4;
using Lattice = dca::phys::models::ThreebandHubbard<PointGroup>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using ModelParameters = dca::phys::params::ModelParameters<Model>;

using BandDmn = dca::func::dmn_0<dca::func::dmn<3, int>>;
using SpinDmn = dca::func::dmn_0<dca::func::dmn<2, int>>;
using BandSpinDmn = dca::func::dmn_variadic<BandDmn, SpinDmn>;
using CDA = dca::phys::ClusterDomainAliases<Lattice::DIMENSION>;
using RClusterDmn = typename CDA::RClusterDmn;

class PackingStub {
public:
  template <typename T>
  int get_buffer_size(const T&) const {
    return sizeof(T);
  }

  template <typename T>
  void pack(char* buffer, const int buffer_size, int& position, const T& value) const {
    if (position + static_cast<int>(sizeof(T)) > buffer_size)
      throw std::runtime_error("PackingStub buffer overflow.");
    std::memcpy(buffer + position, &value, sizeof(T));
    position += sizeof(T);
  }

  template <typename T>
  void unpack(char* buffer, const int buffer_size, int& position, T& value) const {
    if (position + static_cast<int>(sizeof(T)) > buffer_size)
      throw std::runtime_error("PackingStub buffer overflow.");
    std::memcpy(&value, buffer + position, sizeof(T));
    position += sizeof(T);
  }
};

void initializeRealCluster() {
  static bool initialized = false;
  if (!initialized) {
    const std::vector<std::vector<int>> cluster{{2, 0}, {0, 2}};
    dca::phys::domains::cluster_domain_initializer<RClusterDmn>::execute(
        Lattice::initializeRDCABasis(), cluster);
    initialized = true;
  }
}

template <class Container>
bool contains(const Container& container, const int value) {
  return std::find(container.begin(), container.end(), value) != container.end();
}

std::vector<int> uniqueIndices(std::initializer_list<int> indices) {
  std::vector<int> result;
  for (const int index : indices)
    if (!contains(result, index))
      result.push_back(index);
  return result;
}

int indexOfDisplacement(std::vector<double> displacement) {
  const auto& super_basis = RClusterDmn::parameter_type::get_super_basis_vectors();
  const auto& elements = RClusterDmn::parameter_type::get_elements();
  displacement = dca::phys::domains::cluster_operations::translate_inside_cluster(displacement,
                                                                                  super_basis);
  return dca::phys::domains::cluster_operations::index(displacement, elements,
                                                       dca::phys::domains::BRILLOUIN_ZONE);
}

}  // namespace

TEST(ThreebandHubbardTest, ModelParametersSetReadAndPackVpdVpp) {
  ModelParameters pars;
  EXPECT_DOUBLE_EQ(0., pars.get_V_pd());
  EXPECT_DOUBLE_EQ(0., pars.get_V_pp());

  pars.set_t_pd(1.);
  pars.set_t_pp(2.);
  pars.set_ep_d(3.);
  pars.set_ep_p(4.);
  pars.set_U_dd(5.);
  pars.set_U_pp(6.);
  pars.set_V_pd(7.);
  pars.set_V_pp(8.);
  EXPECT_DOUBLE_EQ(7., pars.get_V_pd());
  EXPECT_DOUBLE_EQ(8., pars.get_V_pp());

  PackingStub packing;
  const int buffer_size = pars.getBufferSize(packing);
  std::vector<char> buffer(buffer_size);
  int position = 0;
  pars.pack(packing, buffer.data(), buffer_size, position);
  EXPECT_EQ(buffer_size, position);

  ModelParameters unpacked;
  position = 0;
  unpacked.unpack(packing, buffer.data(), buffer_size, position);
  EXPECT_EQ(buffer_size, position);
  EXPECT_DOUBLE_EQ(1., unpacked.get_t_pd());
  EXPECT_DOUBLE_EQ(2., unpacked.get_t_pp());
  EXPECT_DOUBLE_EQ(3., unpacked.get_ep_d());
  EXPECT_DOUBLE_EQ(4., unpacked.get_ep_p());
  EXPECT_DOUBLE_EQ(5., unpacked.get_U_dd());
  EXPECT_DOUBLE_EQ(6., unpacked.get_U_pp());
  EXPECT_DOUBLE_EQ(7., unpacked.get_V_pd());
  EXPECT_DOUBLE_EQ(8., unpacked.get_V_pp());

  dca::io::JSONReader reader;
  ModelParameters read;
  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/models/analytic_hamiltonians/threeband_hubbard_input.json");
  read.readWrite(reader);
  reader.close_file();
  EXPECT_DOUBLE_EQ(1.5, read.get_t_pd());
  EXPECT_DOUBLE_EQ(0.25, read.get_t_pp());
  EXPECT_DOUBLE_EQ(-1., read.get_ep_d());
  EXPECT_DOUBLE_EQ(2., read.get_ep_p());
  EXPECT_DOUBLE_EQ(8., read.get_U_dd());
  EXPECT_DOUBLE_EQ(4., read.get_U_pp());
  EXPECT_DOUBLE_EQ(0.75, read.get_V_pd());
  EXPECT_DOUBLE_EQ(0.5, read.get_V_pp());
}

TEST(ThreebandHubbardTest, InitializeH0IncludesHartreeShifts) {
  using KDmn = dca::func::dmn_0<dca::func::dmn<2, std::vector<double>>>;
  KDmn::parameter_type::set_elements({{0., 0.}, {M_PI, M_PI}});

  dca::func::function<std::complex<double>,
                      dca::func::dmn_variadic<BandSpinDmn, BandSpinDmn, KDmn>>
      H_0;

  ModelParameters pars;
  pars.set_t_pd(1.5);
  pars.set_t_pp(0.25);
  pars.set_ep_d(-1.);
  pars.set_ep_p(2.);
  pars.set_U_dd(8.);
  pars.set_U_pp(4.);
  pars.set_V_pd(0.75);
  pars.set_V_pp(0.5);

  Lattice::initializeH0WithQ(pars, H_0, std::vector<double>{0., 0.});

  const double d_diagonal = -1. + 8. / 2. + 4. * 0.75;
  const double p_diagonal = 2. + 4. / 2. + 4. * 0.5 + 2. * 0.75;
  for (int s = 0; s < SpinDmn::dmn_size(); ++s) {
    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      EXPECT_DOUBLE_EQ(d_diagonal, H_0(0, s, 0, s, k).real());
      EXPECT_DOUBLE_EQ(p_diagonal, H_0(1, s, 1, s, k).real());
      EXPECT_DOUBLE_EQ(p_diagonal, H_0(2, s, 2, s, k).real());
    }

    EXPECT_NEAR(-2. * 1.5, H_0(0, s, 1, s, 1).imag(), 1e-14);
    EXPECT_NEAR(2. * 1.5, H_0(1, s, 0, s, 1).imag(), 1e-14);
    EXPECT_NEAR(2. * 1.5, H_0(0, s, 2, s, 1).imag(), 1e-14);
    EXPECT_NEAR(-2. * 1.5, H_0(2, s, 0, s, 1).imag(), 1e-14);
    EXPECT_NEAR(4. * 0.25, H_0(1, s, 2, s, 1).real(), 1e-14);
    EXPECT_NEAR(4. * 0.25, H_0(2, s, 1, s, 1).real(), 1e-14);
  }
}

TEST(ThreebandHubbardTest, InitializeHInteractionIncludesVpdVpp) {
  initializeRealCluster();

  dca::func::function<double, dca::func::dmn_variadic<BandSpinDmn, BandSpinDmn, RClusterDmn>>
      H_interaction;

  ModelParameters pars;
  pars.set_U_dd(8.);
  pars.set_U_pp(4.);
  pars.set_V_pd(0.75);
  pars.set_V_pp(0.5);

  Lattice::initializeHInteraction(H_interaction, pars);

  const int origin = RClusterDmn::parameter_type::origin_index();
  const auto& basis = RClusterDmn::parameter_type::get_basis_vectors();
  const int plus_x = indexOfDisplacement(basis[0]);
  const int plus_y = indexOfDisplacement(basis[1]);
  const int minus_x = RClusterDmn::parameter_type::subtract(plus_x, origin);
  const int minus_y = RClusterDmn::parameter_type::subtract(plus_y, origin);

  std::vector<double> plus_x_minus_y(Lattice::DIMENSION, 0.);
  std::vector<double> minus_x_plus_y(Lattice::DIMENSION, 0.);
  for (int d = 0; d < Lattice::DIMENSION; ++d) {
    plus_x_minus_y[d] = basis[0][d] - basis[1][d];
    minus_x_plus_y[d] = basis[1][d] - basis[0][d];
  }

  const auto d_to_px = uniqueIndices({origin, plus_x});
  const auto px_to_d = uniqueIndices({origin, minus_x});
  const auto d_to_py = uniqueIndices({origin, plus_y});
  const auto py_to_d = uniqueIndices({origin, minus_y});
  const auto px_to_py = uniqueIndices({origin, minus_x, plus_y, indexOfDisplacement(minus_x_plus_y)});
  const auto py_to_px = uniqueIndices({origin, plus_x, minus_y, indexOfDisplacement(plus_x_minus_y)});

  for (int r = 0; r < RClusterDmn::dmn_size(); ++r) {
    for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1) {
      for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2) {
        EXPECT_DOUBLE_EQ(r == origin && s1 != s2 ? 8. : 0., H_interaction(0, s1, 0, s2, r));
        EXPECT_DOUBLE_EQ(r == origin && s1 != s2 ? 4. : 0., H_interaction(1, s1, 1, s2, r));
        EXPECT_DOUBLE_EQ(r == origin && s1 != s2 ? 4. : 0., H_interaction(2, s1, 2, s2, r));

        EXPECT_DOUBLE_EQ(contains(d_to_px, r) ? 0.75 : 0., H_interaction(0, s1, 1, s2, r));
        EXPECT_DOUBLE_EQ(contains(px_to_d, r) ? 0.75 : 0., H_interaction(1, s1, 0, s2, r));
        EXPECT_DOUBLE_EQ(contains(d_to_py, r) ? 0.75 : 0., H_interaction(0, s1, 2, s2, r));
        EXPECT_DOUBLE_EQ(contains(py_to_d, r) ? 0.75 : 0., H_interaction(2, s1, 0, s2, r));

        EXPECT_DOUBLE_EQ(contains(px_to_py, r) ? 0.5 : 0., H_interaction(1, s1, 2, s2, r));
        EXPECT_DOUBLE_EQ(contains(py_to_px, r) ? 0.5 : 0., H_interaction(2, s1, 1, s2, r));
      }
    }
  }
}
