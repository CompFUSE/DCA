// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_HAMILTONIAN_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_HAMILTONIAN_H

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_Fock_space.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_Hilbert_spaces/Hilbert_space_phi_representation.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_Hilbert_spaces/Hilbert_space_psi_representation.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/psi_state.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/phi_operators.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_real_axis.h"

namespace DCA {
namespace ADVANCED_EXACT_DIAGONALIZATION {
// DCA::ADVANCED_EXACT_DIAGONALIZATION::

/*
 *   value c+_l c_l
 */
template <typename scalar_type>
struct V_struct {
public:
  long long index;

  scalar_type value;

public:
  void print() {
    std::cout << "\t" << index << "\t" << value << "\n";
  }
};

/*
 *   value c+_i c_j (i><j)
 */
template <typename scalar_type>
struct t_struct {
public:
  long long lhs;
  long long rhs;

  scalar_type value;

public:
  void print() {
    std::cout << "\t" << lhs << "\t" << rhs << "\t" << value << "\n";
  }
};

/*
 *   value n_i n_j
 */
template <typename scalar_type>
struct U_struct {
public:
  long long lhs;
  long long rhs;

  scalar_type value;

public:
  void print() {
    std::cout << "\t" << lhs << "\t" << rhs << "\t" << value << "\n";
  }
};

template <typename parameter_type, typename ed_options>
class fermionic_Hamiltonian {
public:
  typedef typename ed_options::b_dmn b_dmn;
  typedef typename ed_options::s_dmn s_dmn;
  typedef typename ed_options::r_dmn r_dmn;
  typedef typename ed_options::k_dmn k_dmn;

  typedef typename ed_options::profiler_t profiler_t;
  typedef typename ed_options::concurrency_type concurrency_type;

  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  typedef typename ed_options::vector_type vector_type;
  typedef typename ed_options::matrix_type matrix_type;
  typedef typename ed_options::int_matrix_type int_matrix_type;

  typedef typename ed_options::nu_dmn nu_dmn;
  typedef typename ed_options::b_s_r b_s_r_dmn_type;

  typedef typename ed_options::phi_type phi_type;

  typedef Fock_space<parameter_type, ed_options> fermionic_Fock_space_type;
  typedef Hilbert_space<parameter_type, ed_options> Hilbert_space_type;
  typedef Hilbert_space_phi_representation<parameter_type, ed_options> Hilbert_space_phi_representation_type;
  typedef psi_state<parameter_type, ed_options> psi_state_type;

  typedef dmn_0<fermionic_Fock_space_type> fermionic_Fock_dmn_type;

  typedef operators<parameter_type, ed_options> fermionic_operators_type;

  using w_REAL = dmn_0<frequency_domain_real_axis>;

public:
  fermionic_Hamiltonian(parameter_type& parameters_ref);

  void initialize(
      FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_0,
      FUNC_LIB::function<double, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_i);

  void construct_Hamiltonians(bool interacting);
  void diagonalize_Hamiltonians_st();
  // void diagonalize_Hamiltonians_mt();

  void set_spectrum(FUNC_LIB::function<double, w_REAL>& A_w);
  void print_spectrum();

  void print_Hamiltonian(const char* filename);
  void print_eigen_energies(const char* filename);
  void print_eigen_states(const char* filename);

  FUNC_LIB::function<vector_type, fermionic_Fock_dmn_type>& get_eigen_energies() {
    return eigen_energies;
  }
  FUNC_LIB::function<matrix_type, fermionic_Fock_dmn_type>& get_eigen_states() {
    return eigen_states;
  }

  double get_Z();

private:
  void initialize_t_ij_and_U_ij(
      FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_0,
      FUNC_LIB::function<double, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_i);

  void shift_the_energies();

  void add_V_to_Hamiltonian(int N, matrix_type& H, Hilbert_space_type& subspace);
  void add_T_to_Hamiltonian(int N, matrix_type& H, Hilbert_space_type& subspace);
  void add_U_to_Hamiltonian(int N, matrix_type& H, Hilbert_space_type& subspace);

  void add_V_to_Hamiltonian_old(int N, matrix_type& H, Hilbert_space_type& subspace);
  void add_T_to_Hamiltonian_old(int N, matrix_type& H, Hilbert_space_type& subspace);
  void add_U_to_Hamiltonian_old(int N, matrix_type& H, Hilbert_space_type& subspace);

  bool check_block_structure(int N, matrix_type& H, Hilbert_space_type& subspace);

private:
  parameter_type& parameters;
  concurrency_type& concurrency;

  double CUT_OFF;

  std::vector<V_struct<scalar_type>> V_i;
  std::vector<t_struct<complex_type>> t_ij;
  std::vector<U_struct<complex_type>> U_ij;

  FUNC_LIB::function<matrix_type, fermionic_Fock_dmn_type> Hamiltonians;

  FUNC_LIB::function<vector_type, fermionic_Fock_dmn_type> eigen_energies;
  FUNC_LIB::function<matrix_type, fermionic_Fock_dmn_type> eigen_states;
};

template <typename parameter_type, typename ed_options>
fermionic_Hamiltonian<parameter_type, ed_options>::fermionic_Hamiltonian(parameter_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CUT_OFF(parameters.get_eigenvalue_cut_off()),

      V_i(0),
      t_ij(0),
      U_ij(0),

      Hamiltonians("Hamiltonians"),

      eigen_energies("eigen_energies"),
      eigen_states("eigen_states") {}

template <typename parameter_type, typename ed_options>
double fermionic_Hamiltonian<parameter_type, ed_options>::get_Z() {
  assert(eigen_states.size() == fermionic_Fock_dmn_type::dmn_size());

  double beta = parameters.get_beta();

  double Z = 0;
  for (int HS_i = 0; HS_i < eigen_states.size(); ++HS_i)
    for (int n = 0; n < eigen_energies(HS_i).size(); ++n)
      Z += std::exp(-beta * eigen_energies(HS_i)[n]);

  return Z;
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::initialize(
    FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_0,
    FUNC_LIB::function<double, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_i) {
  initialize_t_ij_and_U_ij(H_0, H_i);
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::initialize_t_ij_and_U_ij(
    FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_0,
    FUNC_LIB::function<double, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn>>& H_i) {
  {
    for (int r_j = 0; r_j < r_dmn::dmn_size(); r_j++) {
      for (int s_j = 0; s_j < s_dmn::dmn_size(); s_j++) {
        for (int b_j = 0; b_j < b_dmn::dmn_size(); b_j++) {
          V_struct<scalar_type> V_obj;

          V_obj.index =
              b_j +
              b_dmn::dmn_size() * (s_j + s_dmn::dmn_size() * r_j);  // index: r -> spin -> band

          V_obj.value = -parameters.get_chemical_potential();  // ???: V = -mu = const?

          V_i.push_back(V_obj);
        }
      }
    }
  }

  for (int r_j = 0; r_j < r_dmn::dmn_size(); r_j++) {
    for (int r_i = 0; r_i < r_dmn::dmn_size(); r_i++) {
      int delta_r = r_dmn::parameter_type::subtract(r_j, r_i);

      for (int s_j = 0; s_j < s_dmn::dmn_size(); s_j++) {
        for (int b_j = 0; b_j < b_dmn::dmn_size(); b_j++) {
          for (int s_i = 0; s_i < s_dmn::dmn_size(); s_i++) {
            for (int b_i = 0; b_i < b_dmn::dmn_size(); b_i++) {
              if (abs(H_0(b_i, s_i, b_j, s_j, delta_r)) > 1.e-3) {
                t_struct<complex_type> t_obj;

                t_obj.lhs = b_i + b_dmn::dmn_size() * (s_i + s_dmn::dmn_size() * r_i);
                t_obj.rhs = b_j + b_dmn::dmn_size() * (s_j + s_dmn::dmn_size() * r_j);

                if (true)
                  t_obj.value = real(H_0(b_i, s_i, b_j, s_j, delta_r) / 2.);  // ???: factor 1/2?
                else
                  t_obj.value = real(H_0(b_i, s_i, b_j, s_j, delta_r) / 1.);  // ???: factor 1/2?

                t_ij.push_back(t_obj);
              }

              if (std::abs(H_i(b_i, s_i, b_j, s_j, delta_r)) > 1.e-3) {
                U_struct<complex_type> U_obj;

                U_obj.lhs = b_i + b_dmn::dmn_size() * (s_i + s_dmn::dmn_size() * r_i);
                U_obj.rhs = b_j + b_dmn::dmn_size() * (s_j + s_dmn::dmn_size() * r_j);

                U_obj.value = H_i(b_i, s_i, b_j, s_j, delta_r) / 2.;

                U_ij.push_back(U_obj);
              }
            }
          }
        }
      }
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::construct_Hamiltonians(bool interacting) {
  if (concurrency.id() == 0)
    std::cout << "\n\t" << __FUNCTION__ << std::endl;

  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  if (true) {
    Hamiltonians.reset();
  }

  for (int i = 0; i < Hilbert_spaces.size(); ++i) {
    if (!parameters.do_sector_check() ||
        (Hilbert_spaces[i].get_eigenvalues()[0] == parameters.get_occupation() &&
         Hilbert_spaces[i].get_eigenvalues()[1] == parameters.get_magnetization())) {
      int N = Hilbert_spaces[i].size();
      Hamiltonians(i).resize_no_copy(N);

      matrix_type& H = Hamiltonians(i);
      Hilbert_space_type& subspace = Hilbert_spaces[i];

      for (int l1 = 0; l1 < N; ++l1)
        for (int l0 = 0; l0 < N; ++l0)
          H(l0, l1) = 0;

      add_V_to_Hamiltonian(N, H, subspace);
      add_T_to_Hamiltonian(N, H, subspace);

      if (interacting)
        add_U_to_Hamiltonian(N, H, subspace);
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::add_V_to_Hamiltonian_old(
    int N, matrix_type& H, Hilbert_space_type& subspace) {
  for (int state_f = 0; state_f < N; ++state_f) {
    psi_state_type& Psi_f = subspace.get_element(state_f);
    for (int state_i = 0; state_i < N; ++state_i) {
      psi_state_type& Psi_i = subspace.get_element(state_i);

      for (int k_f = 0; k_f < Psi_f.size(); ++k_f) {
        for (int k_i = 0; k_i < Psi_i.size(); ++k_i) {
          if (Psi_f.phis[k_f] == Psi_i.phis[k_i]) {
            for (size_t l = 0; l < V_i.size(); ++l) {
              H(state_f, state_i) += V_i[l].value * conj(Psi_f.coefficients[k_f]) *
                                     Psi_i.coefficients[k_i] * (1. * Psi_i.phis[k_i][V_i[l].index]);
            }
          }
        }
      }
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::add_V_to_Hamiltonian(
    int /*N*/, matrix_type& H, Hilbert_space_type& subspace) {
  Hilbert_space_phi_representation_type& rep = subspace.get_rep();

  for (int l = 0; l < V_i.size(); ++l) {
    for (int j = 0; j < rep.size(); ++j) {
      int sign = 1;  //  no meaning, but annihilate_at needs it as argument
      phi_type phi = rep.get_phi(j);

      if (fermionic_operators_type::annihilate_at(V_i[l].index, phi, sign)) {
        std::vector<int>& column_index = rep.get_indices(j);
        std::vector<complex_type>& column_alpha = rep.get_alphas(j);

        std::vector<int>& row_index = column_index;
        std::vector<complex_type>& row_alpha = column_alpha;

        for (int c = 0; c < column_index.size(); ++c) {
          for (int r = 0; r < row_index.size(); ++r) {
            H(row_index[r], column_index[c]) += conj(row_alpha[r]) * column_alpha[c] * V_i[l].value;
          }
        }
      }
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::add_T_to_Hamiltonian_old(
    int N, matrix_type& H, Hilbert_space_type& subspace) {
  for (int state_f = 0; state_f < N; ++state_f) {
    psi_state_type& Psi_f = subspace.get_element(state_f);
    for (int state_i = 0; state_i < N; ++state_i) {
      psi_state_type& Psi_i = subspace.get_element(state_i);

      for (int k_f = 0; k_f < Psi_f.size(); ++k_f) {
        for (int k_i = 0; k_i < Psi_i.size(); ++k_i) {
          if (possible_hopping(Psi_i.phis[k_i], Psi_f.phis[k_f])) {
            for (size_t l = 0; l < t_ij.size(); l++) {
              if (Psi_i.phis[k_i][t_ij[l].rhs] == 1 and Psi_i.phis[k_i][t_ij[l].lhs] == 0 and
                  Psi_f.phis[k_f][t_ij[l].rhs] == 0 and Psi_f.phis[k_f][t_ij[l].lhs] == 1) {
                scalar_type phase = 1.;
                phi_type phi_new = Psi_i.phis[k_i];

                {
                  for (int d = 0; d < t_ij[l].rhs; d++)
                    if (phi_new[d] == 1)
                      phase *= scalar_type(-1.);

                  phi_new[t_ij[l].rhs] = 0;
                }

                {
                  for (int d = 0; d < t_ij[l].lhs; d++)
                    if (phi_new[d] == 1)
                      phase *= scalar_type(-1.);

                  phi_new[t_ij[l].lhs] = 1;
                }

                H(state_f, state_i) +=
                    phase * conj(Psi_f.coefficients[k_f]) * Psi_i.coefficients[k_i] * t_ij[l].value;
              }
            }
          }
        }
      }
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::add_T_to_Hamiltonian(
    int /*N*/, matrix_type& H, Hilbert_space_type& subspace) {
  Hilbert_space_phi_representation_type& rep = subspace.get_rep();

  for (int l = 0; l < t_ij.size(); ++l) {
    for (int j = 0; j < rep.size(); ++j) {
      int sign = 1;
      phi_type phi = rep.get_phi(j);

      if (fermionic_operators_type::annihilate_at(t_ij[l].rhs, phi, sign) &&
          fermionic_operators_type::create_at(t_ij[l].lhs, phi, sign)) {
        std::vector<int>& column_index = rep.get_indices(j);
        std::vector<complex_type>& column_alpha = rep.get_alphas(j);

        int i = rep.find(phi);

        if (i < rep.size()) {
          std::vector<int>& row_index = rep.get_indices(i);
          std::vector<complex_type>& row_alpha = rep.get_alphas(i);

          for (int c = 0; c < column_index.size(); ++c) {
            for (int r = 0; r < row_index.size(); ++r) {
              H(row_index[r], column_index[c]) +=
                  conj(row_alpha[r]) * column_alpha[c] * scalar_type(sign) * t_ij[l].value;
            }
          }
        }
      }
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::add_U_to_Hamiltonian_old(
    int N, matrix_type& H, Hilbert_space_type& subspace) {
  for (int state = 0; state < N; ++state) {
    psi_state_type& Psi = subspace.get_element(state);

    for (int k = 0; k < Psi.size(); ++k) {
      for (int l = 0; l < U_ij.size(); ++l) {
        H(state, state) +=
            U_ij[l].value / scalar_type(4.) * conj(Psi.coefficients[k]) * Psi.coefficients[k];

        if (abs(Psi.phis[k][U_ij[l].lhs] - scalar_type(1.)) < 1.e-3)
          H(state, state) -=
              U_ij[l].value / scalar_type(2.) * conj(Psi.coefficients[k]) * Psi.coefficients[k];

        if (abs(Psi.phis[k][U_ij[l].rhs] - scalar_type(1.)) < 1.e-3)
          H(state, state) -=
              U_ij[l].value / scalar_type(2.) * conj(Psi.coefficients[k]) * Psi.coefficients[k];

        if (abs(Psi.phis[k][U_ij[l].lhs] - scalar_type(1.)) < 1.e-3 and
            abs(Psi.phis[k][U_ij[l].rhs] - scalar_type(1.)) < 1.e-3)
          H(state, state) += U_ij[l].value * conj(Psi.coefficients[k]) * Psi.coefficients[k];
      }
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::add_U_to_Hamiltonian(
    int N, matrix_type& H, Hilbert_space_type& subspace) {
  Hilbert_space_phi_representation_type& rep = subspace.get_rep();

  for (int l = 0; l < U_ij.size(); ++l) {
    // n_i*n_j term
    for (int j = 0; j < rep.size(); ++j) {
      int sign = 1;  // no meaning, but annihilate_at and create_an need argument
      phi_type phi = rep.get_phi(j);

      if (fermionic_operators_type::annihilate_at(U_ij[l].rhs, phi, sign) &&
          fermionic_operators_type::create_at(U_ij[l].rhs, phi, sign) &&
          fermionic_operators_type::annihilate_at(U_ij[l].lhs, phi, sign) &&
          fermionic_operators_type::create_at(U_ij[l].lhs, phi, sign)) {
        std::vector<int>& column_index = rep.get_indices(j);
        std::vector<complex_type>& column_alpha = rep.get_alphas(j);

        int i = rep.find(phi);
        if (i < rep.size()) {
          std::vector<int>& row_index = rep.get_indices(i);
          std::vector<complex_type>& row_alpha = rep.get_alphas(i);

          for (int c = 0; c < column_index.size(); ++c) {
            for (int r = 0; r < row_index.size(); ++r) {
              H(row_index[r], column_index[c]) +=
                  conj(row_alpha[r]) * column_alpha[c] * U_ij[l].value;
            }
          }
        }
      }
    }

    // n_i term
    for (int j = 0; j < rep.size(); ++j) {
      int sign = 1;  // no meaning, but annihilate_at and create_an need argument
      phi_type phi = rep.get_phi(j);

      if (fermionic_operators_type::annihilate_at(U_ij[l].rhs, phi, sign) &&
          fermionic_operators_type::create_at(U_ij[l].rhs, phi, sign)) {
        std::vector<int>& column_index = rep.get_indices(j);
        std::vector<complex_type>& column_alpha = rep.get_alphas(j);

        int i = rep.find(phi);
        if (i < rep.size()) {
          std::vector<int>& row_index = rep.get_indices(i);
          std::vector<complex_type>& row_alpha = rep.get_alphas(i);

          for (int c = 0; c < column_index.size(); ++c) {
            for (int r = 0; r < row_index.size(); ++r) {
              H(row_index[r], column_index[c]) -=
                  conj(row_alpha[r]) * column_alpha[c] * U_ij[l].value / scalar_type(2.);
            }
          }
        }
      }
    }

    // n_j term
    for (int j = 0; j < rep.size(); ++j) {
      int sign = 1;  // no meaning, but annihilate_at and create_an need argument
      phi_type phi = rep.get_phi(j);

      if (fermionic_operators_type::annihilate_at(U_ij[l].lhs, phi, sign) &&
          fermionic_operators_type::create_at(U_ij[l].lhs, phi, sign)) {
        std::vector<int>& column_index = rep.get_indices(j);
        std::vector<complex_type>& column_alpha = rep.get_alphas(j);

        int i = rep.find(phi);
        if (i < rep.size()) {
          std::vector<int>& row_index = rep.get_indices(i);
          std::vector<complex_type>& row_alpha = rep.get_alphas(i);

          for (int c = 0; c < column_index.size(); ++c) {
            for (int r = 0; r < row_index.size(); ++r) {
              H(row_index[r], column_index[c]) -=
                  conj(row_alpha[r]) * column_alpha[c] * U_ij[l].value / scalar_type(2.);
            }
          }
        }
      }
    }
  }

  // const term
  std::complex<scalar_type> U_ij_sum = 0;

  for (int l = 0; l < U_ij.size(); ++l)
    U_ij_sum += U_ij[l].value;

  U_ij_sum /= scalar_type(4.);

  for (int state = 0; state < N; ++state)
    H(state, state) += U_ij_sum;
}

template <typename parameter_type, typename ed_options>
bool fermionic_Hamiltonian<parameter_type, ed_options>::check_block_structure(
    int N, matrix_type& H, Hilbert_space_type& subspace) {
  for (int state_f = 0; state_f < N; ++state_f) {
    psi_state_type& Psi_f = subspace.get_element(state_f);
    for (int state_i = 0; state_i < N; ++state_i) {
      psi_state_type& Psi_i = subspace.get_element(state_i);

      // cut_off?
      if (Psi_f.get_eigenvalues() != Psi_i.get_eigenvalues() &&
          (abs(H(state_f, state_i).real()) > ed_options::get_epsilon() ||
           abs(H(state_f, state_i).imag()) > ed_options::get_epsilon()))
        return false;
    }
  }

  return true;
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::diagonalize_Hamiltonians_st() {
  if (concurrency.id() == 0)
    std::cout << "\n\t" << __FUNCTION__ << "\n\n";

  int start = clock();

  eigen_energies.reset();
  eigen_states.reset();

  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i) {
    int N = Hilbert_spaces[i].size();

    std::vector<int> Hilbert_space_evals = Hilbert_spaces[i].get_eigenvalues();

    eigen_energies(i).resize(N);
    eigen_states(i).resize_no_copy(N);

    {
      if (concurrency.id() == 0)
        std::cout << "\t N_occ : " << Hilbert_space_evals[0] << ", Sz : " << Hilbert_space_evals[1]
                  << ", \t size : " << N << ", \t time : ";

      int start = clock();
      LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', Hamiltonians(i), eigen_energies(i),
                                           eigen_states(i));
      int end = clock();

      if (concurrency.id() == 0)
        std::cout << double(end - start) / double(CLOCKS_PER_SEC) << "\n";
    }
  }

  int end = clock();

  if (concurrency.id() == 0) {
    std::cout << "\n\t" << __FUNCTION__
              << "\t total time : " << double(end - start) / double(CLOCKS_PER_SEC) << "\n\n";

    print_spectrum();
  }

  shift_the_energies();
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::print_spectrum() {
  double N_BINS = 20;

  double E_MIN = 0.;
  double E_MAX = 0.;

  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
    for (int n = 0; n < Hilbert_spaces[i].size(); n++)
      E_MIN = E_MIN > eigen_energies(i)[n] ? eigen_energies(i)[n] : E_MIN;

  for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
    for (int n = 0; n < Hilbert_spaces[i].size(); n++)
      E_MAX = E_MAX < eigen_energies(i)[n] ? eigen_energies(i)[n] : E_MAX;

  if (concurrency.id() == 0)
    std::cout << "\n\t E_min : " << E_MIN << "\t E_max : " << E_MAX << "\n\n";

  if (E_MAX - E_MIN > 1.e-2) {
    size_t NUMBER_OF_LAMBDAS = 0;

    for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
      NUMBER_OF_LAMBDAS += eigen_energies(i).size();

    double delta = (E_MAX - E_MIN) / (N_BINS - 1.);
    std::vector<size_t> y(N_BINS, 0);

    for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
      for (int n = 0; n < Hilbert_spaces[i].size(); n++)
        y[int((eigen_energies(i)[n] - E_MIN) / delta)] += 1;

    if (concurrency.id() == 0) {
      std::cout << "\n\t distribution of the energies : \n\n";
      for (int l = 0; l < N_BINS; l++)
        std::cout << "\t" << E_MIN + delta / 2. + l * delta << "\t"
                  << std::string(int(double(y[l]) / double(NUMBER_OF_LAMBDAS) * 400.), '*')
                  << std::endl;
      std::cout << "\n\n";
    }
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::set_spectrum(
    FUNC_LIB::function<double, w_REAL>& A_w) {
  A_w = 0;

  std::vector<double>& w_elem = w_REAL::get_elements();

  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
    for (int n = 0; n < Hilbert_spaces[i].size(); ++n)
      for (int w_ind = 0; w_ind < w_REAL::dmn_size() - 1; ++w_ind)
        if (w_elem[w_ind] <= eigen_energies(i)[n] and eigen_energies(i)[n] < w_elem[w_ind + 1])
          A_w(w_ind) += 1.;

  double total = 0;
  for (int w_ind = 0; w_ind < w_REAL::dmn_size() - 1; w_ind++)
    total += A_w(w_ind);

  A_w /= total;
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::shift_the_energies() {
  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  double E_0 = 0.;
  for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
    for (int n = 0; n < Hilbert_spaces[i].size(); n++)
      E_0 = E_0 > eigen_energies(i)[n] ? eigen_energies(i)[n] : E_0;

  // subtract ground state energy E_0
  for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
    for (int n = 0; n < Hilbert_spaces[i].size(); n++)
      eigen_energies(i)[n] -= E_0;

  {
    double beta = parameters.get_beta();

    int number_of_eigenvalues = 0;
    for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
      for (int n = 0; n < Hilbert_spaces[i].size(); n++)
        if (std::exp(-beta * eigen_energies(i)[n]) > CUT_OFF)
          number_of_eigenvalues += 1;

    int total_of_eigenvalues = 0;
    for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
      for (int n = 0; n < Hilbert_spaces[i].size(); n++)
        total_of_eigenvalues += 1;

    if (concurrency.id() == 0)
      std::cout << "\n\n\t number of eigenvalues exp(-beta*lambda) > CUT_OFF: "
                << number_of_eigenvalues << ", " << total_of_eigenvalues << "\n\n";
  }
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::print_Hamiltonian(const char* filename) {
  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  int HS_i = 0;

  for (int i = 0; i < Hilbert_spaces.size(); ++i) {
    if (Hilbert_spaces[i].get_occupation() == parameters.get_n_0() &&
        Hilbert_spaces[i].get_magnetization() == parameters.get_Sz_0())
      HS_i = i;
  }

  std::cout << "Print Hamiltonian of Hilbert-space #" << HS_i << std::endl;

  std::ofstream data;
  data.open(filename);

  for (int m = 0; m < Hilbert_spaces[HS_i].size(); ++m) {
    for (int n = 0; n < Hilbert_spaces[HS_i].size(); ++n) {
      assert(abs(imag(Hamiltonians(HS_i)(m, n))) < ed_options::get_epsilon());
      data << real(Hamiltonians(HS_i)(m, n)) << "\t";
    }
    data << "\n";
  }
  data.close();
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::print_eigen_energies(const char* filename) {
  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  std::cout << "Print eigen-energies: " << std::endl;

  std::ofstream data;
  data.open(filename);

  for (int i = 0; i < fermionic_Fock_dmn_type::dmn_size(); ++i)
    for (int n = 0; n < Hilbert_spaces[i].size(); ++n)
      data << eigen_energies(i)[n] << "\n";

  data.close();
}

template <typename parameter_type, typename ed_options>
void fermionic_Hamiltonian<parameter_type, ed_options>::print_eigen_states(const char* filename) {
  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  int HS_i = 0;

  for (int i = 0; i < Hilbert_spaces.size(); ++i) {
    if (Hilbert_spaces[i].get_occupation() == parameters.get_n_0() &&
        Hilbert_spaces[i].get_magnetization() == parameters.get_Sz_0())
      HS_i = i;
  }

  std::cout << "Print real eigen-states of Hilbert-space #" << HS_i << std::endl;

  std::ofstream data;
  data.open(filename);

  for (int m = 0; m < Hilbert_spaces[HS_i].size(); ++m) {
    for (int n = 0; n < Hilbert_spaces[HS_i].size(); ++n) {
      assert(abs(imag(eigen_states(HS_i)(m, n))) < 1.e-10);
      data << real(eigen_states(HS_i)(m, n)) << "\t";
    }
    data << "\n";
  }
  data.close();
}

}  // ADVANCED_EXACT_DIAGONALIZATION
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_HAMILTONIAN_H
