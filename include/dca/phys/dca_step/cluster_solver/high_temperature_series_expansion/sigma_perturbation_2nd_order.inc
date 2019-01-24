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
// This class computes the second order term of the perturbation expansion of the self-energy.

template <class parameters_type, class k_dmn_t>
class sigma_perturbation<2, parameters_type, k_dmn_t> {
public:
  using concurrency_type = typename parameters_type::concurrency_type;
  using Threading = typename parameters_type::ThreadingType;
  using ThisType = sigma_perturbation<2, parameters_type, k_dmn_t>;

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
  void threaded_execute_on_cluster(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G);

  // void execute_on_cluster_2(func::function<std::complex<double>, func::dmn_variadic<nu,nu,
  // k_dmn_t,
  // w> >&
  // G);

  template <typename Writer>
  void write(Writer& writer);

private:
  void execute_2A(func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G);
  void execute_2B(func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G);

  void threaded_execute_2B(int id, int num_threads, sp_function_type& S, const sp_function_type& G,
                           const chi_function_type& chi) const;

  int subtract_freq_bf(int, int);  // fermion-boson

protected:
  parameters_type& parameters;
  concurrency_type& concurrency;

  U_function_type& U;

  chi_function_type& chi;
  phi_function_type& phi;

  func::function<std::complex<double>, func::dmn_variadic<b, b, k_dmn_t, w_VERTEX_BOSONIC>> U_chi_U;

  sp_function_type Sigma;

  sp_function_type Sigma_2A;
  sp_function_type Sigma_2B;

private:
  struct sigma_perturbation_data {
    U_function_type* U_ptr;

    sp_function_type* G_ptr;
    sp_function_type* S_ptr;

    chi_function_type* chi_ptr;
    phi_function_type* phi_ptr;

    concurrency_type* concurrency_ptr;
  };
};

template <class parameters_type, class k_dmn_t>
sigma_perturbation<2, parameters_type, k_dmn_t>::sigma_perturbation(
    parameters_type& parameters_ref, compute_interaction& interaction_obj,
    compute_bubble<ph, parameters_type, k_dmn_t, w>& chi_obj,
    compute_bubble<pp, parameters_type, k_dmn_t, w>& phi_obj)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      U(interaction_obj.get_function()),

      chi(chi_obj.get_function()),
      phi(phi_obj.get_function()),

      U_chi_U("U-chi-U"),

      Sigma("Sigma-2nd-order"),
      Sigma_2A("Sigma-2nd-order-A"),
      Sigma_2B("Sigma-2nd-order-B") {}

template <class parameters_type, class k_dmn_t>
template <typename Writer>
void sigma_perturbation<2, parameters_type, k_dmn_t>::write(Writer& /*writer*/) {}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<2, parameters_type, k_dmn_t>::execute_on_cluster(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G) {
  if (concurrency.id() == concurrency.first())
    std::cout << __FUNCTION__ << std::endl;

  // cout << "\t U : " << U(0,0,0,1) << endl;

  sigma_perturbation<2, parameters_type, k_dmn_t>::execute_2B(G);

  // Sigma_2A *= 1.;
  // Sigma_2B *= 2.;

  Sigma = 0.;
  // Sigma += Sigma_2A;
  Sigma += Sigma_2B;

  if (true) {
    std::complex<double> I(0, 1);
    for (int b_ind = 0; b_ind < 2 * b::dmn_size(); ++b_ind) {
      for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
        int wc_ind = w::dmn_size() / 8;

        double wc = w::get_elements()[wc_ind];

        std::complex<double> Sigma_wc = Sigma(b_ind, b_ind, k_ind, wc_ind);

        double alpha = real(Sigma_wc);
        double beta = imag(Sigma_wc * wc);

        for (int w_ind = 0; w_ind < wc_ind; ++w_ind) {
          Sigma(b_ind, b_ind, k_ind, w_ind) = alpha + beta * I / w::get_elements()[w_ind];
          Sigma(b_ind, b_ind, k_ind, w::dmn_size() - 1 - w_ind) =
              alpha - beta * I / w::get_elements()[w_ind];
        }
      }
    }
  }
}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<2, parameters_type, k_dmn_t>::execute_2A(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G) {
  if (concurrency.id() == concurrency.first())
    std::cout << __FUNCTION__ << std::endl;

  double U_value = U(0, 0, 0, 1);

  Sigma_2A = 0.;

  for (int nu1_ind = 0; nu1_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu1_ind) {
    for (int q1_ind = 0; q1_ind < k_dmn_t::dmn_size(); ++q1_ind) {
      for (int nu2_ind = 0; nu2_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu2_ind) {
        for (int q2_ind = 0; q2_ind < k_dmn_t::dmn_size(); ++q2_ind) {
          for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
            int w_minus_nu1 = subtract_freq_bf(nu1_ind, w_ind);
            int w_minus_nu2 = subtract_freq_bf(nu2_ind, w_ind);
            int w_minus_nu1_minus_nu2 = subtract_freq_bf(nu2_ind, w_minus_nu1);

            if (w_minus_nu1 < 0 || w_minus_nu1 >= w::dmn_size() || w_minus_nu2 < 0 ||
                w_minus_nu2 >= w::dmn_size() || w_minus_nu1_minus_nu2 < 0 ||
                w_minus_nu1_minus_nu2 >= w::dmn_size())
              continue;

            for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
              int k_minus_q1 = k_dmn_t::parameters_type::subtract(q1_ind, k_ind);
              int k_minus_q2 = k_dmn_t::parameters_type::subtract(q2_ind, k_ind);
              int k_minus_q1_minus_q2 = k_dmn_t::parameters_type::subtract(q2_ind, k_minus_q1);

              Sigma_2A(0, 0, 0, 0, k_ind, w_ind) +=
                  G(0, 0, 0, 0, k_minus_q1, w_minus_nu1) * G(0, 0, 0, 0, k_minus_q2, w_minus_nu2) *
                  G(0, 0, 0, 0, k_minus_q1_minus_q2, w_minus_nu1_minus_nu2);
            }
          }
        }
      }
    }
  }

  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind)
      Sigma_2A(0, 1, 0, 1, k_ind, w_ind) = Sigma_2A(0, 0, 0, 0, k_ind, w_ind);

  double factor = 1. / (parameters.get_beta() * parameters.get_beta() * k_dmn_t::dmn_size() *
                        k_dmn_t::dmn_size()) *
                  U_value * U_value;
  Sigma_2A *= factor;
}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<2, parameters_type, k_dmn_t>::execute_2B(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G) {
  if (concurrency.id() == concurrency.first())
    std::cout << __FUNCTION__ << std::endl;

  double U_value = U(0, 0, 0, 1);

  Sigma_2B = 0.;

  for (int q_ind = 0; q_ind < k_dmn_t::dmn_size(); ++q_ind) {
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
      int k_minus_q = k_dmn_t::parameters_type::subtract(q_ind, k_ind);

      for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
        int nu_c = (nu_ind - w_VERTEX_BOSONIC::dmn_size() / 2);

        for (int w_ind = std::fabs(nu_c); w_ind < w::dmn_size() - std::fabs(nu_c); ++w_ind) {
          int w_minus_nu = w_ind - nu_c;

          Sigma_2B(0, 0, 0, 0, k_ind, w_ind) +=
              G(0, 0, 0, 0, k_minus_q, w_minus_nu) * chi(0, 0, 0, 0, q_ind, nu_ind);
        }
      }
    }
  }

  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind)
      Sigma_2B(0, 1, 0, 1, k_ind, w_ind) = Sigma_2B(0, 0, 0, 0, k_ind, w_ind);

  double factor = 1. / (parameters.get_beta() * k_dmn_t::dmn_size()) * U_value * U_value;
  Sigma_2B *= factor;
}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<2, parameters_type, k_dmn_t>::threaded_execute_on_cluster(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_dmn_t, w>>& G) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t\t second-order Self-energy \n\n" << std::endl;

  //       std::cout << "\t U : " << U(0,0,0,1) << std::endl;

  int nr_threads = parameters.get_hts_threads();

  Threading threads;
  auto task = std::bind(&ThisType::threaded_execute_2B, this, std::placeholders::_1,
                        std::placeholders::_2, std::ref(Sigma_2B), std::ref(G), std::ref(chi));
  threads.execute(nr_threads, task);

  {
    double U_value = U(0, 0, 0, 1);
    double factor = 1. / (parameters.get_beta() * k_dmn_t::dmn_size()) * U_value * U_value;
    Sigma_2B *= factor;

    concurrency.sum(Sigma_2B);
  }

  Sigma = Sigma_2B;

  if (true) {
    std::complex<double> I(0, 1);
    for (int b_ind = 0; b_ind < 2 * b::dmn_size(); ++b_ind) {
      for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
        int wc_ind = w::dmn_size() / 8;

        double wc = w::get_elements()[wc_ind];

        std::complex<double> Sigma_wc = Sigma(b_ind, b_ind, k_ind, wc_ind);

        double alpha = real(Sigma_wc);
        double beta = imag(Sigma_wc * wc);

        for (int w_ind = 0; w_ind < wc_ind; ++w_ind) {
          Sigma(b_ind, b_ind, k_ind, w_ind) = alpha + beta * I / w::get_elements()[w_ind];
          Sigma(b_ind, b_ind, k_ind, w::dmn_size() - 1 - w_ind) =
              alpha - beta * I / w::get_elements()[w_ind];
        }
      }
    }
  }
}

template <class parameters_type, class k_dmn_t>
void sigma_perturbation<2, parameters_type, k_dmn_t>::threaded_execute_2B(
    const int id, const int nr_threads, sp_function_type& S, const sp_function_type& G,
    const chi_function_type& chi) const {
  k_dmn_t q_dmn;
  std::pair<int, int> q_bounds = concurrency.get_bounds(q_dmn);

  k_dmn_t k_dmn;
  std::pair<int, int> k_bounds = dca::parallel::util::getBounds(id, nr_threads, k_dmn);

  for (int k_ind = k_bounds.first; k_ind < k_bounds.second; k_ind++) {
    double percentage = double(k_ind - k_bounds.first) / double(k_bounds.second - k_bounds.first);

    if (concurrency.id() == concurrency.first() and id == 0 and (int(100 * percentage) % 10 == 0))
      std::cout << "\t" << int(100 * percentage) << " % finished\t" << dca::util::print_time()
                << "\n";

    for (int q_ind = q_bounds.first; q_ind < q_bounds.second; q_ind++) {
      int k_minus_q = k_dmn_t::parameter_type::subtract(q_ind, k_ind);

      for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
        int nu_c = (nu_ind - w_VERTEX_BOSONIC::dmn_size() / 2);

        for (int w_ind = std::fabs(nu_c); w_ind < w::dmn_size() - std::fabs(nu_c); ++w_ind) {
          int w_minus_nu = w_ind - nu_c;

          S(0, 0, 0, 0, k_ind, w_ind) +=
              G(0, 0, 0, 0, k_minus_q, w_minus_nu) * chi(0, 0, 0, 0, q_ind, nu_ind);  //*factor;
        }
      }
    }

    for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
      S(0, 1, 0, 1, k_ind, w_ind) = S(0, 0, 0, 0, k_ind, w_ind);
  }

  //       for(int w_ind=0; w_ind<w::dmn_size(); ++w_ind)
  // 	for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind)
  // 	  S(0,1, 0,1, k_ind, w_ind) = S(0,0, 0,0, k_ind, w_ind);

  //       Sigma_2B *= factor;
}

template <class parameters_type, class k_dmn_t>
int sigma_perturbation<2, parameters_type, k_dmn_t>::subtract_freq_bf(int w1, int w2) {
  int w_f = 2 * (w2 - w::dmn_size() / 2) + 1;             // transform fermionic
  int w_b = 2 * (w1 - w_VERTEX_BOSONIC::dmn_size() / 2);  // transform bosonic
  int res = ((w_f - w_b) - 1 + w::dmn_size()) / 2;        // result is fermionic
  return res;
}
