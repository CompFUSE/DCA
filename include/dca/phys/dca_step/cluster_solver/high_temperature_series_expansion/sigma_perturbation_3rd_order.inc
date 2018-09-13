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
// This class computes the third order term of the perturbation expansion of the self-energy.

template <class parameters_type, class k_dmn_t>
class sigma_perturbation<3, parameters_type, k_dmn_t> {
public:
  using w = func::dmn_0<domains::frequency_domain>;
  using w_VERTEX_BOSONIC = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED_BOSONIC>>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using ph_bubble_t = compute_bubble<ph, parameters_type, k_dmn_t, w>;
  // INTERNAL: Shouldn't the template argument be pp instead of ph?
  using pp_bubble_t = compute_bubble<ph, parameters_type, k_dmn_t, w>;

  using chi_function_type = typename ph_bubble_t::function_type;
  using phi_function_type = typename pp_bubble_t::function_type;
  using sp_function_type =
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>;
  using U_function_type = typename compute_interaction::function_type;

public:
  sigma_perturbation(parameters_type& parameters_ref, compute_interaction& interaction_obj,
                     compute_bubble<ph, parameters_type, k_dmn_t, w>& chi_obj,
                     compute_bubble<pp, parameters_type, k_dmn_t, w>& phi_obj);

  sp_function_type& get_function() {
    return Sigma;
  }

  void execute_on_cluster(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G);

  template <typename Writer>
  void write(Writer& writer);

private:
  void execute_RPA(func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G);
  void execute_VC(func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G);

  int subtract_freq_fb(int, int);

protected:
  parameters_type& parameters;

  U_function_type& U;

  chi_function_type& chi;
  phi_function_type& phi;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>> Sigma;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>> Sigma_RPA;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>> Sigma_VC;
};

template <class parameters_type, class k_dmn_t>
sigma_perturbation<3, parameters_type, k_dmn_t>::sigma_perturbation(
    parameters_type& parameters_ref, compute_interaction& interaction_obj,
    compute_bubble<ph, parameters_type, k_dmn_t, w>& chi_obj,
    compute_bubble<pp, parameters_type, k_dmn_t, w>& phi_obj)
    : parameters(parameters_ref),

      U(interaction_obj.get_function()),

      chi(chi_obj.get_function()),
      phi(phi_obj.get_function()),

      Sigma("Sigma-3rd-order"),
      Sigma_RPA("Sigma-3rd-order-RPA"),
      Sigma_VC("Sigma-3rd-order-VC") {}

template <class parameters_type, class k_dmn_t>
template <typename Writer>
void sigma_perturbation<3, parameters_type, k_dmn_t>::write(Writer& /*writer*/) {}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<3, parameters_type, k_dmn_t>::execute_on_cluster(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G) {
  std::cout << __FUNCTION__ << std::endl;

  std::cout << "\t U : " << U(0, 0, 0, 1) << std::endl;

  sigma_perturbation<3, parameters_type, k_dmn_t>::execute_RPA(G);
  sigma_perturbation<3, parameters_type, k_dmn_t>::execute_VC(G);

  Sigma = 0.;
  Sigma += Sigma_RPA;
  Sigma += Sigma_VC;
}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<3, parameters_type, k_dmn_t>::execute_RPA(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G) {
  std::cout << __FUNCTION__ << std::endl;

  double U_value = U(0, 0, 0, 1);

  Sigma_RPA = 0.;

  for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
    for (int q_ind = 0; q_ind < k_dmn_t::dmn_size(); ++q_ind) {
      int nu_c = (nu_ind - w_VERTEX_BOSONIC::dmn_size() / 2);

      for (int w_ind = std::fabs(nu_c); w_ind < w::dmn_size() - std::fabs(nu_c); ++w_ind) {
        for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
          int k_minus_q = k_dmn_t::parameters_type::subtract(q_ind, k_ind);
          int w_minus_nu = w_ind - nu_c;
          Sigma_RPA(0, 0, 0, 0, k_ind, w_ind) += G(0, 0, 0, 0, k_minus_q, w_minus_nu) *
                                                 chi(0, 0, 0, 0, q_ind, nu_ind) *
                                                 chi(0, 0, 0, 0, q_ind, nu_ind);
        }
      }
    }
  }
  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind)
      Sigma_RPA(0, 1, 0, 1, k_ind, w_ind) = Sigma_RPA(0, 0, 0, 0, k_ind, w_ind);

  double factor = 1. / (parameters.get_beta() * k_dmn_t::dmn_size()) * U_value * U_value * U_value;
  Sigma_RPA *= factor;
}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<3, parameters_type, k_dmn_t>::execute_VC(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G) {
  std::cout << __FUNCTION__ << std::endl;

  double U_value = U(0, 0, 0, 1);

  Sigma_VC = 0.;

  for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
    for (int q_ind = 0; q_ind < k_dmn_t::dmn_size(); ++q_ind) {
      for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
        int nu_minus_w = subtract_freq_fb(w_ind, nu_ind);
        if (nu_minus_w < 0 || nu_minus_w >= w::dmn_size())
          continue;

        for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
          int q_minus_k = k_dmn_t::parameters_type::subtract(k_ind, q_ind);

          Sigma_VC(0, 0, 0, 0, k_ind, w_ind) += G(0, 0, 0, 0, q_minus_k, nu_minus_w) *
                                                phi(0, 0, 0, 0, q_ind, nu_ind) *
                                                phi(0, 0, 0, 0, q_ind, nu_ind);
        }
      }
    }
  }

  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind)
      Sigma_VC(0, 1, 0, 1, k_ind, w_ind) = Sigma_VC(0, 0, 0, 0, k_ind, w_ind);

  double factor = 1. / (parameters.get_beta() * k_dmn_t::dmn_size()) * U_value * U_value * U_value;
  Sigma_VC *= factor;
}

template <class parameters_type, class k_dmn_t>
int sigma_perturbation<3, parameters_type, k_dmn_t>::subtract_freq_fb(int w1, int w2) {
  int w_f = 2 * (w1 - w::dmn_size() / 2) + 1;             // transform fermionic
  int w_b = 2 * (w2 - w_VERTEX_BOSONIC::dmn_size() / 2);  // transform bosonic
  int res = ((w_b - w_f) - 1 + w::dmn_size()) / 2;        // result is fermionic
  return res;
}
