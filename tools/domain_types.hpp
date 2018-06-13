// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// The file serves as a reference for a consistent naming convention for the various 'physical'
// domains used throughout the code. In particular it should never be included!
// For more type defintions see also deprecated/python_build/prog_dca++/type_defintions.h.

#error This file should not be included!

using DCA_cluster_family_type =
    cluster_domain_family<double, lattice_type::DIMENSION, CLUSTER, BRILLOUIN_ZONE>;
using HOST_sp_cluster_family_type =
    cluster_domain_family<double, lattice_type::DIMENSION, LATTICE_SP, BRILLOUIN_ZONE>;
using HOST_tp_cluster_family_type =
    cluster_domain_family<double, lattice_type::DIMENSION, LATTICE_TP, BRILLOUIN_ZONE>;

using LDA_k_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, MOMENTUM_SPACE, PARALLELLEPIPEDUM>;
using DCA_k_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
using PCM_k_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
using host_k_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
using host_vertex_k_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, LATTICE_TP, MOMENTUM_SPACE, BRILLOUIN_ZONE>;

using LDA_r_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, REAL_SPACE, PARALLELLEPIPEDUM>;
using DCA_r_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, CLUSTER, REAL_SPACE, BRILLOUIN_ZONE>;
using PCM_r_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, CLUSTER, REAL_SPACE, BRILLOUIN_ZONE>;
using host_r_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, REAL_SPACE, BRILLOUIN_ZONE>;
using host_vertex_r_cluster_type =
    cluster_domain<double, lattice_type::DIMENSION, LATTICE_TP, REAL_SPACE, BRILLOUIN_ZONE>;

using k_LDA = func::dmn_0<LDA_k_cluster_type>;
using k_DCA = func::dmn_0<DCA_k_cluster_type>;
using k_PCM = func::dmn_0<PCM_k_cluster_type>;
using k_HOST = func::dmn_0<host_k_cluster_type>;
using k_HOST_VERTEX = func::dmn_0<host_vertex_k_cluster_type>;

using r_LDA = func::dmn_0<LDA_r_cluster_type>;
using r_DCA = func::dmn_0<DCA_r_cluster_type>;
using r_PCM = func::dmn_0<PCM_r_cluster_type>;
using r_HOST = func::dmn_0<host_r_cluster_type>;
using r_HOST_VERTEX = func::dmn_0<host_vertex_r_cluster_type>;

using crystal_harmonics_expansion = centered_cluster_domain<r_HOST_VERTEX::parameter_type>;
using crystal_harmonics_expansion_dmn_t = func::dmn_0<crystal_harmonics_expansion>;

using s = func::dmn_0<electron_spin_domain>;
using b = func::dmn_0<electron_band_domain>;

using pn__b_r_DCA_e_spin_domain_type = particle_number_domain<b, r_DCA, s>;
using pn__b_k_DCA_e_spin_domain_type = particle_number_domain<b, k_DCA, s>;
using pn__b_r_DCA_s = func::dmn_0<pn__b_r_DCA_e_spin_domain_type>;
using pn__b_k_DCA_s = func::dmn_0<pn__b_r_DCA_e_spin_domain_type>;

using sp_vertex_time_domain_type = DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN>;
using tp_vertex_time_domain_type = DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN>;
using sp_vertex_time_domain_pos_type = DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN_POSITIVE>;
using tp_vertex_time_domain_pos_type = DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN_POSITIVE>;

using compact_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::COMPACT>;
using extended_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::EXTENDED>;
using compact_positive_vertex_frequency_domain_type =
    DCA::vertex_frequency_domain<DCA::COMPACT_POSITIVE>;
using extended_positive_vertex_frequency_domain_type =
    DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>;
using bosonic_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>;

using t = func::dmn_0<time_domain>;
using w = func::dmn_0<frequency_domain>;
using w_REAL = func::dmn_0<frequency_domain_real_axis>;
using w_IMAG = func::dmn_0<frequency_domain_imag_axis>;
using sp_time_dmn_t = func::dmn_0<sp_vertex_time_domain_type>;
using tp_time_dmn_t = func::dmn_0<tp_vertex_time_domain_type>;
using sp_time_pos_dmn_t = func::dmn_0<sp_vertex_time_domain_pos_type>;
using tp_time_pos_dmn_t = func::dmn_0<tp_vertex_time_domain_pos_type>;
using w_VERTEX = func::dmn_0<compact_vertex_frequency_domain_type>;
using w_VERTEX_EXTENDED = func::dmn_0<extended_vertex_frequency_domain_type>;
using w_VERTEX_EXTENDED_POS = func::dmn_0<extended_positive_vertex_frequency_domain_type>;
using w_VERTEX_BOSONIC = func::dmn_0<bosonic_vertex_frequency_domain_type>;

using nu = func::dmn_variadic<b, s>;  // orbital-spin index
using b_r_DCA = func::dmn_variadic<b, r_DCA>;
using b_k_DCA = func::dmn_variadic<b, k_DCA>;
using nu_nu = func::dmn_variadic<nu, nu>;
using nu_nu__nu_nu = func::dmn_variadic<nu_nu, nu_nu>;
using nu_r_DCA = func::dmn_variadic<nu, r_DCA>;
using nu_k_DCA = func::dmn_variadic<nu, k_DCA>;
using nu_w_REAL = func::dmn_variadic<nu, w_REAL>;
using b_b = func::dmn_variadic<b, b>;
using b_b__b_b = func::dmn_variadic<b_b, b_b>;
using s_s = func::dmn_variadic<s, s>;
using w_k_DCA = func::dmn_variadic<w, k_DCA>;
using k_DCA_w = func::dmn_variadic<k_DCA, w>;
using k_DCA_t = func::dmn_variadic<k_DCA, t>;
using t_k_DCA = func::dmn_variadic<t, k_DCA>;
using k_DCA_w_REAL = func::dmn_variadic<k_DCA, w_REAL>;
using w_REAL_k_DCA = func::dmn_variadic<w_REAL, k_DCA>;
using k_DCA_w_VERTEX = func::dmn_variadic<k_DCA, w_VERTEX>;
using k_DCA_w_VERTEX_EXTENDED = func::dmn_variadic<k_DCA, w_VERTEX_EXTENDED>;
using nu_k_LDA = func::dmn_variadic<nu, k_LDA>;
using r_DCA_r_DCA = func::dmn_variadic<r_DCA, r_DCA>;
using r_DCA_k_DCA = func::dmn_variadic<r_DCA, k_DCA>;
using k_DCA_k_DCA = func::dmn_variadic<k_DCA, k_DCA>;
using k_DCA_r_DCA = func::dmn_variadic<k_DCA, r_DCA>;

using b_b_r_DCA = func::dmn_variadic<b, b, r_DCA>;
using b_b_k_DCA = func::dmn_variadic<b, b, k_DCA>;
using b_r_DCA_s = func::dmn_variadic<b, r_DCA, s>;
using b_k_DCA_s = func::dmn_variadic<b, k_DCA, s>;
using b_s__k_LDA = func::dmn_variadic<b, s, k_LDA>;
using b_b__s__k_LDA = func::dmn_variadic<b_b, s, k_LDA>;
using nu_nu_k_LDA = func::dmn_variadic<nu, nu, k_LDA>;
using nu_nu_r_LDA = func::dmn_variadic<nu, nu, r_LDA>;
using nu_nu_r_DCA = func::dmn_variadic<nu, nu, r_DCA>;
using nu_nu_k_DCA = func::dmn_variadic<nu, nu, k_DCA>;
using nu_nu_w = func::dmn_variadic<nu, nu, w>;
using nu_nu_w_REAL = func::dmn_variadic<nu, nu, w_REAL>;

using b_b_r_DCA_r_DCA = func::dmn_variadic<b, b, r_DCA, r_DCA>;
using nu_nu_k_DCA_w = func::dmn_variadic<nu, nu, k_DCA, w>;
using nu_nu_k_DCA_t = func::dmn_variadic<nu, nu, k_DCA, t>;
using nu_nu_r_DCA_t = func::dmn_variadic<nu, nu, r_DCA, t>;
using nu_nu_r_DCA_w = func::dmn_variadic<nu, nu, r_DCA, w>;
using nu_nu_k_HOST_w = func::dmn_variadic<nu, nu, k_HOST, w>;
using nu_nu_k_HOST_t = func::dmn_variadic<nu, nu, k_HOST, t>;
using nu_nu_r_HOST_t = func::dmn_variadic<nu, nu, r_HOST, t>;
using nu_nu_r_HOST_w = func::dmn_variadic<nu, nu, r_HOST, w>;
using t_r_DCA_nu_nu = func::dmn_variadic<t, r_DCA, nu, nu>;
using b_b_b_b = func::dmn_variadic<b, b, b, b>;
using nu_nu_k_DCA_w_REAL = func::dmn_variadic<nu, nu, k_DCA, w_REAL>;
using nu_nu_k_DCA_w_IMAG = func::dmn_variadic<nu, nu, k_DCA, w_IMAG>;

using nu_nu__nu_nu__k_DCA = func::dmn_variadic<nu, nu, nu, nu, k_DCA>;
using b_b_r_DCA_w_VERTEX_w_VERTEX = func::dmn_variadic<b, b, r_DCA, w_VERTEX, w_VERTEX>;

using nu_nu__nu_nu__k_DCA_w_VERTEX_BOSONIC =
    func::dmn_variadic<nu, nu, nu, nu, k_DCA, w_VERTEX_BOSONIC>;
using b_b__b_b__k_DCA_w_VERTEX_BOSONIC = func::dmn_variadic<b, b, b, b, k_DCA, w_VERTEX_BOSONIC>;
using b_b_k_DCA_b_b_k_DCA = func::dmn_variadic<b, b, k_DCA, b, b, k_DCA>;

// Band structure plot
using brillouin_zone_cut_domain_type = brillouin_zone_cut_domain<101>;
using k_domain_cut_dmn_type = func::dmn_0<brillouin_zone_cut_domain_type>;
using nu_k_cut = func::dmn_variadic<nu, k_domain_cut_dmn_type>;

using DCA_iteration_domain_type = func::dmn_0<DCA_iteration_domain>;

// Analysis
using b_b_r_DCA_w_VERTEX = func::dmn_variadic<b, b, r_DCA, w_VERTEX>;
using b_b_k_DCA_w_VERTEX = func::dmn_variadic<b, b, k_DCA, w_VERTEX>;
using b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX =
    func::dmn_variadic<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX>;
using b_b_k_DCA_w_VERTEX_EXTENDED = func::dmn_variadic<b, b, k_DCA, w_VERTEX_EXTENDED>;
using b_b_k_DCA_w_VERTEX_EXTENDED__b_b_k_DCA_w_VERTEX_EXTENDED =
    func::dmn_variadic<b_b_k_DCA_w_VERTEX_EXTENDED, b_b_k_DCA_w_VERTEX_EXTENDED>;
