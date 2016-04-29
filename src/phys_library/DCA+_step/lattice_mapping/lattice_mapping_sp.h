//-*-C++-*-
//TODO  Why is there a reference to high_temperature_serie_solver ??
#ifndef DCA_LATTICE_MAP_SP_H
#define DCA_LATTICE_MAP_SP_H

#include "phys_library/DCA+_step/lattice_mapping/interpolation/interpolation_sp.h"
#include "phys_library/DCA+_step/lattice_mapping/deconvolution/deconvolution_sp.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/high_temperature_series_expansion_solver.h"

namespace DCA {
/*! \ingroup LATTICE-MAPPING
 *
 *  \author Peter Staar
 *  \brief  This class implements the lattice_map.
 *
 */
template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class lattice_map_sp {
  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

public:
  lattice_map_sp(parameters_type& parameters_ref);
  ~lattice_map_sp();

  void execute(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, source_k_dmn_t, w>>& f_source,
               FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& Sigma_interp,
               FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& Sigma_deconv,
               FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_target);

  template <typename MOMS_type, typename HTS_solver_type, typename coarsegraining_sp_type>
  void execute_with_HTS_approximation(
      MOMS_type& HTS_MOMS, HTS_solver_type& HTS_solver, coarsegraining_sp_type& cluster_mapping_obj,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, source_k_dmn_t, w>>& f_source,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& Sigma_interp,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& Sigma_deconv,
      FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_target);

private:
  void initialize();

  template <typename k_dmn_t>
  void plot_function(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_dmn_t, w>>& f);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  interpolation_sp<parameters_type, source_k_dmn_t, target_k_dmn_t> interpolation_obj;
  deconvolution_sp<parameters_type, source_k_dmn_t, target_k_dmn_t> deconvolution_obj;

  FUNC_LIB::function<double, nu> Sigma_shift;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> Sigma_cluster;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_HOST, w>> Sigma_lattice;
};

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
lattice_map_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::lattice_map_sp(
    parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      interpolation_obj(parameters),
      deconvolution_obj(parameters),

      Sigma_shift("Sigma_shift"),

      Sigma_cluster("Sigma_cluster_HTS"),
      Sigma_lattice("Sigma_lattice_HTS") {}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
lattice_map_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::~lattice_map_sp() {}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void lattice_map_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, source_k_dmn_t, w>>& f_source,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_interp,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_approx,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_target) {
  symmetrize::execute(f_source);

  // plot_function(f_source);

  interpolation_obj.execute_with_alpha_transformation(f_source, f_interp);

  // plot_function(f_interp);

  symmetrize::execute(f_interp);

  deconvolution_obj.execute(f_source, f_interp, f_approx, f_target);

  // plot_function(f_target);

  symmetrize::execute(f_target);
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename MOMS_type, typename HTS_solver_type, typename coarsegraining_sp_type>
void lattice_map_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::execute_with_HTS_approximation(
    MOMS_type& HTS_MOMS, HTS_solver_type& HTS_solver, coarsegraining_sp_type& cluster_mapping_obj,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, source_k_dmn_t, w>>& f_source,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_interp,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_approx,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, target_k_dmn_t, w>>& f_target) {
  {
    HTS_solver.initialize(0);

    HTS_solver.execute();

    HTS_solver.finalize();

    HTS_MOMS.Sigma = f_source;

    cluster_mapping_obj.compute_S_K_w(HTS_MOMS.Sigma_lattice, HTS_MOMS.Sigma_cluster);

    Sigma_cluster = HTS_MOMS.Sigma_cluster;
    Sigma_lattice = HTS_MOMS.Sigma_lattice;

    //       HTS_MOMS.write("data_HTS.json");
  }

  plot_function(Sigma_lattice);

  {
    f_source -= Sigma_cluster;

    execute(f_source, f_interp, f_approx, f_target);

    f_source += Sigma_cluster;
    f_interp += Sigma_lattice;
    f_approx += Sigma_lattice;
    f_target += Sigma_lattice;
  }

  plot_function(f_interp);
  plot_function(f_target);
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename k_dmn_t>
void lattice_map_sp<parameters_type, source_k_dmn_t, target_k_dmn_t>::plot_function(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_dmn_t, w>>& f) {
  std::vector<double> x(0);
  std::vector<double> y(0);

  std::vector<double> z_re(0);
  std::vector<double> z_im(0);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
    x.push_back(k_dmn_t::get_elements()[k_ind][0]);
    y.push_back(k_dmn_t::get_elements()[k_ind][1]);

    z_re.push_back(real(f(0, 0, k_ind, w::dmn_size() / 2)));
    z_im.push_back(imag(f(0, 0, k_ind, w::dmn_size() / 2)));
  }

  SHOW::heatmap(x, y, z_re, f.get_name());
  SHOW::heatmap(x, y, z_im, f.get_name());
}
}

#endif
