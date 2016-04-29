//-*-C++-*-
// INTERNAL UNUSED file
#ifndef SS_CT_HYB_ACCUMULATOR_LEGENDRE_H
#define SS_CT_HYB_ACCUMULATOR_LEGENDRE_H
#include "phys_library/domain_types.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/qmci_accumulator.h"
using namespace types;

namespace DCA {
namespace QMCI {

/*!
 *  \brief   This class organizes the measurements in the single-site hybridization QMC
 *  \author  Peter Staar
 *  \author  Bart Ydens
 *  \version 1.0
 *
 * The impurity self-energy can be expressed in the following form:
 * \f{eqnarray*}{
 * \Sigma_{ab}(i \nu) = \frac12 \sum_{ij} G^{-1}_{ai}(i \nu)(U_{jb} + U_{bj}) F_{ib}^j(i \nu)
 * \f}
 * Where the impurity Green function \f$G_{ab}(\tau)\f$ is measured according to
 * \f{eqnarray*}{
 * G_{ab}(\tau) = -\frac{1}{\beta} \left< \sum_{\alpha \beta = 1}^{n} M_{\alpha \beta}
 * \delta^{-}(\tau - (\tau'_{\alpha} - \tau_{\beta}) \delta_{a,\alpha} \delta_{b,\beta} \right>_{MC}
 * \f}
 * and the correlation function \f$F_{ab}^j(\tau)\f$ is measured according to
 * \f{eqnarray*}{
 * F_{ab}^j(\tau) = -\frac{1}{\beta} \left< \sum_{\alpha \beta = 1}^{n} M_{\alpha \beta}
 * \delta^{-}(\tau - (\tau'_{\alpha} - \tau_{\beta}) n_j(\tau_{\beta}) \delta_{a,\alpha}
 * \delta_{b,\beta} \right>_{MC}
 * \f}
 * where
 * \f{eqnarray*}{
 * \delta^{\pm}(\tau) = \sum_{n} (\pm1)^n \delta(\tau - n\beta)
 * \f}
 *
 * These functions can also be measured in Legendre repesentation,
 * \f{eqnarray*}{
 * G(\tau) = \sum_{l \ge 0}\frac{\sqrt{2l+1}}{\beta} P_l(x(\tau))G_l
 * \f}
 * where \f$P_l(x)\f$ are the legendre polynomials, and \f$G_l\f$ denote the coefficients of
 * \f$G(\tau)\f$ in the Legendre basis and can be measured according to
 * \f{eqnarray*}{
 * G_{ab;l} = -\frac{\sqrt{2l+1}}{\beta} \left< \sum_{\alpha \beta = 1}^{n} M_{\alpha \beta}
 * P_l(\tau'_{\alpha} - \tau_{\beta}) \delta_{a,\alpha} \delta_{b,\beta} \right>_{MC}
 * \f}
 * and equivalent for \f$F_{ab}^j(\tau)\f$
 * \f{eqnarray*}{
 * F_{ab;l}^j = -\frac{\sqrt{2l+1}}{\beta} \left< \sum_{\alpha \beta = 1}^{n} M_{\alpha \beta}
 * P_l(\tau'_{\alpha} - \tau_{\beta})n_j(\tau_{\beta}) \delta_{a,\alpha} \delta_{b,\beta}
 * \right>_{MC}
 * \f}
 * In the segment picture, the matrix element \f$n_j(\tau_{\beta})\f$ (one or zero) of this operator
 * is simply determined by examining whether or not a segment is present in the timeline for flavor
 * j at time \f$\tau_{\beta}\f$. Hence this quantity can be measured
 * at essentially no additional computational cost. Note that this function only contributes if j is
 * different from the index b (and \f$\beta\f$) and therefore \f$n_j\f$ is never evaluated at the
 * position of the creator of flavor b at time  \f$\tau_{\beta}\f$.
 *
 * We prefer to accumulate
 * \f{eqnarray*}{
 * (G \Sigma)_{ab}(\tau) = \frac12 \sum_{j}(U_{jb}+U_{bj})F_{ab}^j(\tau)
 * \f}
 * directly instead of the individual quantities \f$F_{ab}^j(\tau)\f$
 */

template <class parameters_type, class base_cluster_type>
class MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type, base_cluster_type> {
  typedef double scalar_type;

  typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
  typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

  typedef dmn_0<r_cluster_type> r_dmn_t;
  typedef dmn_0<k_cluster_type> k_dmn_t;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::Concurrency_Type concurrency_type;

  typedef legendre_domain<time_domain, SINGLE_PARTICLE_QUANTITY> legendre_t;
  typedef dmn_0<legendre_t> legendre_dmn_t;

public:
  MC_single_particle_accumulator(parameters_type& parameters_ref);

  ~MC_single_particle_accumulator();

  void initialize(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>

                                     >& G_r_w,
                  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w);

  template <class orbital_configuration_type, class vertex_vertex_matrix_type, class H_type>
  void accumulate(FUNC_LIB::function<orbital_configuration_type, nu>& vertices,
                  FUNC_LIB::function<bool, nu>& has_full_line,
                  FUNC_LIB::function<vertex_vertex_matrix_type, nu>& M, double sign,
                  H_type& H_interactions);

  void compute(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>

                                  >& G_r_w,
               FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w);

  void finalize();

  template <class stream_type>
  void to_JSON(stream_type& ss);

private:
  template <class orbital_configuration_type, class H_type>
  int compute_U_times_n(FUNC_LIB::function<orbital_configuration_type, nu>& vertices,
                        FUNC_LIB::function<bool, nu>& has_full_line, H_type& H_interactions,
                        double t_start, int flavor);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  int N_spin_orbitals;

  std::vector<double> P;

  FUNC_LIB::function<double, dmn_4<legendre_dmn_t, nu, nu, r_dmn_t>> G_l_r;
  FUNC_LIB::function<double, dmn_4<nu, nu, r_dmn_t, legendre_dmn_t>> G_r_l;

  FUNC_LIB::function<double, dmn_4<legendre_dmn_t, nu, nu, r_dmn_t>> GS_l_r;
  FUNC_LIB::function<double, dmn_4<nu, nu, r_dmn_t, legendre_dmn_t>> GS_r_l;
};

template <class parameters_type, class base_cluster_type>
MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type,
                               base_cluster_type>::MC_single_particle_accumulator(parameters_type&
                                                                                      parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      N_spin_orbitals(b::dmn_size() * s::dmn_size()),

      P(parameters.get_nb_of_legendre_coefficients_single_particle(), 0),

      G_l_r("G_legendre_r"),
      G_r_l("G_r_legendre"),

      GS_l_r("GS_legendre_r"),
      GS_r_l("GS_r_legendre") {}

template <class parameters_type, class base_cluster_type>
MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type,
                               base_cluster_type>::~MC_single_particle_accumulator() {}

template <class parameters_type, class base_cluster_type>
void MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type, base_cluster_type>::initialize(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>

                       >& G_r_w,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w) {
  for (int i = 0; i < G_l_r.

                      size();

       i++)
    G_l_r(i) = 0;

  for (int i = 0; i < G_r_l.

                      size();

       i++)
    G_r_l(i) = 0;

  for (int i = 0; i < GS_l_r.

                      size();

       i++)
    GS_l_r(i) = 0;

  for (int i = 0; i < GS_r_l.

                      size();

       i++)
    GS_r_l(i) = 0;

  for (int i = 0; i < G_r_w.

                      size();

       i++)
    G_r_w(i) = 0;

  for (int i = 0; i < GS_r_w.

                      size();

       i++)
    GS_r_w(i) = 0;
}

template <class parameters_type, class base_cluster_type>
void MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type, base_cluster_type>::finalize() {
  concurrency.sum_and_average(G_r_l, parameters.get_number_of_measurements());
  concurrency.sum_and_average(GS_r_l, parameters.get_number_of_measurements());

  if (concurrency.id() == 0) {
    for (int l = 0; l < legendre_dmn_t::dmn_size(); l++) {
      cout << "\t" << l;
      for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
        cout << "\t" << G_r_l(0, 0, r_ind, l);
      }
      cout << "\n";
    }
  }
}

template <class parameters_type, class base_cluster_type>
template <class stream_type>
void MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type, base_cluster_type>::to_JSON(
    stream_type& ss) {
  ss << ",";
  G_r_l.to_JSON(ss);
  ss << ",";
  GS_r_l.to_JSON(ss);
}

template <class parameters_type, class base_cluster_type>
template <class orbital_configuration_type, class vertex_vertex_matrix_type, class H_type>
void MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type, base_cluster_type>::accumulate(
    FUNC_LIB::function<orbital_configuration_type, nu>& vertices,
    FUNC_LIB::function<bool, nu>& has_full_line,
    FUNC_LIB::function<vertex_vertex_matrix_type, nu>& M, double sign, H_type& H_interactions) {
  throw std::logic_error(__FUNCTION__);
}

template <class parameters_type, class base_cluster_type>
template <class orbital_configuration_type, class H_type>
int MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type, base_cluster_type>::compute_U_times_n(
    FUNC_LIB::function<orbital_configuration_type, nu>& vertices,
    FUNC_LIB::function<bool, nu>& has_full_line, H_type& H_interactions, double t_start, int flavor) {
  // typedef typename configuration_type::orbital_configuration_t orbital_configuration_type;
  int U_times_n = 0;

  for (int j = 0; j < N_spin_orbitals; j++) {
    if (j == flavor)
      continue;

    if (has_full_line(j)) {
      U_times_n += 1. / 2. * (H_interactions(j, flavor, 0) + H_interactions(flavor, j, 0));
    }
    else {
      orbital_configuration_type& vertices_j = vertices(j);
      for (int y = 0; y < int(vertices_j.size()); y++) {
        if (vertices_j[y].t_start() < vertices_j[y].t_end() && vertices_j[y].t_start() < t_start &&
            t_start < vertices_j[y].t_end()) {
          U_times_n += 1. / 2. * (H_interactions(j, flavor, 0) + H_interactions(flavor, j, 0));
        }
        else if (vertices_j[y].t_end() < vertices_j[y].t_start() &&
                 (vertices_j[y].t_start() < t_start || t_start < vertices_j[y].t_end())) {
          U_times_n += 1. / 2. * (H_interactions(j, flavor, 0) + H_interactions(flavor, j, 0));
        }
      }
    }
  }

  return U_times_n;
}

template <class parameters_type, class base_cluster_type>
void MC_single_particle_accumulator<SS_CT_HYB, LEGENDRE, parameters_type, base_cluster_type>::compute(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>

                       >& G_r_w,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w) {
  for (int l = 0; l < legendre_dmn_t::dmn_size(); l++) {
    double renorm = (2 * l + 1) / parameters.get_beta();

    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu = 0;

           nu < b::dmn_size() * s::dmn_size();

           nu++) {
        G_r_l(nu, nu, r_ind, l) = renorm * G_l_r(l, nu, nu, r_ind);
        GS_r_l(nu, nu, r_ind, l) = renorm * GS_l_r(l, nu, nu, r_ind);
      }
    }
  }

  FT<legendre_dmn_t, w>::execute(G_r_l, G_r_w);
  FT<legendre_dmn_t, w>::execute(GS_r_l, GS_r_w);

  double one_div_n_sites = 1. / double(base_cluster_type::get_cluster_size());
  G_r_w *= one_div_n_sites;
  GS_r_w *= one_div_n_sites;
}

}  // namespace QMC
}  // namespace DCA
#endif
