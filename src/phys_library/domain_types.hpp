#ifndef PHYS_LIBRARY_DOMAIN_TYPES_H
#define PHYS_LIBRARY_DOMAIN_TYPES_H

#include "lattice_types.hpp"  //this file is copied to the applications folder
// INTERNAL question: does this file need to be self contained?
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/cluster/cluster_domain_family.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_real_axis.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_imag_axis.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"
#include "phys_library/domains/time_and_frequency/vertex_time_domain.h"
#include "phys_library/domains/cluster/centered_cluster_domain.h"
#include "phys_library/domains/Quantum_domain/brillouin_zone_cut_domain.h"
#include "phys_library/domains/Quantum_domain/DCA_iteration_domain.h"
// TODO  check which domains to include
#include "comp_library/function_library/domains/special_domains/dmn_0.h"
#include "comp_library/function_library/domains/product_domains/dmn_2.h"
#include "comp_library/function_library/domains/product_domains/dmn_3.h"
#include "comp_library/function_library/domains/product_domains/dmn_4.h"
#include "comp_library/function_library/domains/product_domains/dmn_5.h"
#include "comp_library/function_library/domains/product_domains/dmn_6.h"
#include "comp_library/function_library/domains/product_domains/dmn_8.h"

namespace types {

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

using compact_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::COMPACT>;
using extended_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::EXTENDED>;
using compact_positive_vertex_frequency_domain_type =
    DCA::vertex_frequency_domain<DCA::COMPACT_POSITIVE>;
using extended_positive_vertex_frequency_domain_type =
    DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>;
using bosonic_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>;

using s = dmn_0<electron_spin_domain>;
using b = dmn_0<electron_band_domain>;
using t = dmn_0<time_domain>;
using w = dmn_0<frequency_domain>;
using w_REAL = dmn_0<frequency_domain_real_axis>;
using w_IMAG = dmn_0<frequency_domain_imag_axis>;
using k_LDA = dmn_0<LDA_k_cluster_type>;
using k_DCA = dmn_0<DCA_k_cluster_type>;
using k_PCM = dmn_0<PCM_k_cluster_type>;
using k_HOST = dmn_0<host_k_cluster_type>;
using k_HOST_VERTEX = dmn_0<host_vertex_k_cluster_type>;
using r_LDA = dmn_0<LDA_r_cluster_type>;
using r_DCA = dmn_0<DCA_r_cluster_type>;
using r_PCM = dmn_0<PCM_r_cluster_type>;
using r_HOST = dmn_0<host_r_cluster_type>;
using r_HOST_VERTEX = dmn_0<host_vertex_r_cluster_type>;
using w_VERTEX = dmn_0<compact_vertex_frequency_domain_type>;
using w_VERTEX_EXTENDED = dmn_0<extended_vertex_frequency_domain_type>;
using w_VERTEX_EXTENDED_POS = dmn_0<extended_positive_vertex_frequency_domain_type>;
using w_VERTEX_BOSONIC = dmn_0<bosonic_vertex_frequency_domain_type>;

using tp_vertex_time_domain_pos_type = DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN_POSITIVE>;
using tp_time_pos_dmn_t = dmn_0<tp_vertex_time_domain_pos_type>;

using crystal_harmonics_expansion = centered_cluster_domain<r_HOST_VERTEX::parameter_type>;
using crystal_harmonics_expansion_dmn_t = dmn_0<crystal_harmonics_expansion>;

using nu = dmn_2<b, s>;  // orbital-spin index
using b_r_DCA = dmn_2<b, r_DCA>;
using b_k_DCA = dmn_2<b, k_DCA>;
using nu_nu = dmn_2<nu, nu>;
using nu_nu__nu_nu = dmn_2<nu_nu, nu_nu>;
using nu_r_DCA = dmn_2<nu, r_DCA>;
using nu_k_DCA = dmn_2<nu, k_DCA>;
using nu_w_REAL = dmn_2<nu, w_REAL>;
using b_b = dmn_2<b, b>;
using b_b__b_b = dmn_2<b_b, b_b>;
using s_s = dmn_2<s, s>;
using w_k_DCA = dmn_2<w, k_DCA>;
using k_DCA_w = dmn_2<k_DCA, w>;
using k_DCA_t = dmn_2<k_DCA, t>;
using t_k_DCA = dmn_2<t, k_DCA>;
using k_DCA_w_REAL = dmn_2<k_DCA, w_REAL>;
using w_REAL_k_DCA = dmn_2<w_REAL, k_DCA>;
using k_DCA_w_VERTEX = dmn_2<k_DCA, w_VERTEX>;
using k_DCA_w_VERTEX_EXTENDED = dmn_2<k_DCA, w_VERTEX_EXTENDED>;
using nu_k_LDA = dmn_2<nu, k_LDA>;
using r_DCA_r_DCA = dmn_2<r_DCA, r_DCA>;
using r_DCA_k_DCA = dmn_2<r_DCA, k_DCA>;
using k_DCA_k_DCA = dmn_2<k_DCA, k_DCA>;
using k_DCA_k_PCM = dmn_2<k_DCA, k_PCM>;
using k_PCM_k_DCA = dmn_2<k_PCM, k_DCA>;
using r_PCM_r_PCM = dmn_2<r_PCM, r_PCM>;
using k_PCM_k_PCM = dmn_2<k_PCM, k_PCM>;
using r_PCM_k_PCM = dmn_2<r_PCM, k_PCM>;
using r_DCA_r_PCM = dmn_2<r_DCA, r_PCM>;
using k_PCM_r_PCM = dmn_2<k_PCM, r_PCM>;
using k_DCA_r_DCA = dmn_2<k_DCA, r_DCA>;

using b_b_r_DCA = dmn_3<b, b, r_DCA>;
using b_b_k_DCA = dmn_3<b, b, k_DCA>;
using b_b_r_PCM = dmn_3<b, b, r_PCM>;
using b_b_k_PCM = dmn_3<b, b, k_PCM>;
using b_r_DCA_s = dmn_3<b, r_DCA, s>;
using b_k_DCA_s = dmn_3<b, k_DCA, s>;
using b_s__k_LDA = dmn_3<b, s, k_LDA>;
using b_b__s__k_LDA = dmn_3<b_b, s, k_LDA>;
using nu_nu_k_LDA = dmn_3<nu, nu, k_LDA>;
using nu_nu_r_LDA = dmn_3<nu, nu, r_LDA>;
using nu_nu_r_DCA = dmn_3<nu, nu, r_DCA>;
using nu_nu_k_DCA = dmn_3<nu, nu, k_DCA>;
using nu_nu_r_PCM = dmn_3<nu, nu, r_PCM>;
using nu_nu_k_PCM = dmn_3<nu, nu, k_PCM>;
using nu_nu_w = dmn_3<nu, nu, w>;
using nu_nu_w_REAL = dmn_3<nu, nu, w_REAL>;
using b_b_r_DCA_r_DCA = dmn_4<b, b, r_DCA, r_DCA>;
using nu_nu_k_DCA_w = dmn_4<nu, nu, k_DCA, w>;
using nu_nu_k_DCA_t = dmn_4<nu, nu, k_DCA, t>;
using nu_nu_r_DCA_t = dmn_4<nu, nu, r_DCA, t>;
using nu_nu_r_DCA_w = dmn_4<nu, nu, r_DCA, w>;
using nu_nu_k_PCM_w = dmn_4<nu, nu, k_PCM, w>;
using nu_nu_k_PCM_t = dmn_4<nu, nu, k_PCM, t>;
using nu_nu_r_PCM_t = dmn_4<nu, nu, r_PCM, t>;
using nu_nu_r_PCM_w = dmn_4<nu, nu, r_PCM, w>;
using nu_nu_k_HOST_w = dmn_4<nu, nu, k_HOST, w>;
using nu_nu_k_HOST_t = dmn_4<nu, nu, k_HOST, t>;
using nu_nu_r_HOST_t = dmn_4<nu, nu, r_HOST, t>;
using nu_nu_r_HOST_w = dmn_4<nu, nu, r_HOST, w>;
using t_r_DCA_nu_nu = dmn_4<t, r_DCA, nu, nu>;
using b_b_b_b = dmn_4<b, b, b, b>;
using nu_nu_k_DCA_w_REAL = dmn_4<nu, nu, k_DCA, w_REAL>;
using nu_nu_k_DCA_w_IMAG = dmn_4<nu, nu, k_DCA, w_IMAG>;
using nu_nu__nu_nu__k_DCA = dmn_5<nu, nu, nu, nu, k_DCA>;
using b_b_r_DCA_w_VERTEX_w_VERTEX = dmn_5<b, b, r_DCA, w_VERTEX, w_VERTEX>;
using b_b_r_PCM_r_PCM_w_VERTEX_w_VERTEX = dmn_6<b, b, r_PCM, r_PCM, w_VERTEX, w_VERTEX>;
using b_b_k_PCM_k_PCM_w_VERTEX_w_VERTEX = dmn_6<b, b, k_PCM, k_PCM, w_VERTEX, w_VERTEX>;
using b_b_k_DCA_b_b_k_DCA = dmn_6<b, b, k_DCA, b, b, k_DCA>;
using b_b__b_b__k_PCM_k_PCM_w_VERTEX_w_VERTEX = dmn_8<b, b, b, b, k_PCM, k_PCM, w_VERTEX, w_VERTEX>;

// Band-structure plot
using brillouin_zone_cut_domain_type = brillouin_zone_cut_domain<101>;
using k_domain_cut_dmn_type = dmn_0<brillouin_zone_cut_domain_type>;
using nu_k_cut = dmn_2<nu, k_domain_cut_dmn_type>;

using DCA_iteration_domain_type = dmn_0<DCA_iteration_domain>;
}  // types

#endif  // PHYS_LIBRARY_DOMAIN_TYPES_H
