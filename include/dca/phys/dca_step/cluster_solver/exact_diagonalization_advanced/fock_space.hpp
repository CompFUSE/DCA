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
// This class implements the fermionic Fock space.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_FOCK_SPACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_FOCK_SPACE_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/function/function.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/psi_state.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/symmetry_operation.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hilbert_spaces/hilbert_space.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/cluster_symmetry.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameters_type, typename ed_options>  // N: size of bitset sequence
class Fock_space {
public:
  typedef typename ed_options::b_dmn b_dmn;
  typedef typename ed_options::s_dmn s_dmn;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;
    
  typedef typename ed_options::profiler_t profiler_t;
  typedef typename ed_options::concurrency_type concurrency_type;

  typedef typename ed_options::int_type int_type;
  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  typedef typename ed_options::vector_type vector_type;
  typedef typename ed_options::matrix_type matrix_type;
  typedef typename ed_options::int_matrix_type int_matrix_type;

  typedef typename ed_options::nu_dmn nu_dmn;
  typedef typename ed_options::b_s_r b_s_r_dmn_type;

  typedef Hilbert_space<parameters_type, ed_options> element_type;
  typedef Fock_space<parameters_type, ed_options> this_type;

  typedef typename RClusterDmn::parameter_type r_cluster_type;

  typedef domains::cluster_symmetry<r_cluster_type> r_symmetry_type;

  typedef typename r_symmetry_type::symmetry_matrix_dmn_t r_symmetry_matrix_dmn_t;

  typedef typename r_symmetry_type::b_dmn_t b_dmn_t;
  typedef typename r_symmetry_type::c_dmn_t c_dmn_t;

  typedef typename r_symmetry_type::sym_super_cell_dmn_t sym_super_cell_dmn_t;

public:
  Fock_space(bool occupation_symmetry, bool magnetization_symmetry = false);

  void print_subspaces(bool full = false);

  static int get_size();
  static std::string get_name();
  static std::vector<element_type>& get_elements();

  void apply_rotation_symmetry(std::string symmetries, std::string ED_method);

  void apply_translation_symmetry(std::string ED_method = "default");

  bool check_orthogonality();

  void initialize_rep();

private:
  void initialize_without_symmetry();
  void initialize_with_occupation_symmetry();
  void apply_magnetization_symmetry();

  void sort_wrt_size();

  template <class symmetry_operation>
  void factorize(symmetry_operation& Op, element_type& subspace,
                 std::vector<element_type>& new_Hilbert_spaces, bool create_subspaces, int k);
};

template <typename parameter_type, typename ed_options>
bool operator>(const Hilbert_space<parameter_type, ed_options>& Hilbert_space_1,
               const Hilbert_space<parameter_type, ed_options>& Hilbert_space_2) {
  return Hilbert_space_1.size() > Hilbert_space_2.size();
}

template <typename parameter_type, typename ed_options>
Fock_space<parameter_type, ed_options>::Fock_space(bool occupation_symmetry,
                                                   bool magnetization_symmetry) {
  std::vector<element_type>& Hilbert_spaces = get_elements();

  if (occupation_symmetry) {
    initialize_with_occupation_symmetry();

    if (magnetization_symmetry)
      apply_magnetization_symmetry();
  }

  else
    initialize_without_symmetry();

  get_elements() = Hilbert_spaces;
}

template <typename parameter_type, typename ed_options>
int Fock_space<parameter_type, ed_options>::get_size() {
  return get_elements().size();
}

template <typename parameter_type, typename ed_options>
std::string Fock_space<parameter_type, ed_options>::get_name() {
  return "Fock-space-domain";
}

template <typename parameter_type, typename ed_options>
std::vector<Hilbert_space<parameter_type, ed_options>>& Fock_space<parameter_type,
                                                                   ed_options>::get_elements() {
  static std::vector<element_type> elements;
  return elements;
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::initialize_without_symmetry() {
  std::vector<element_type>& Hilbert_spaces = get_elements();

  int_type num_fermionic_states = b_dmn::dmn_size() * s_dmn::dmn_size() * RClusterDmn::dmn_size();
  int_type num_fock_states = 1ul << num_fermionic_states;

  std::vector<std::string> names;
  names.push_back("full");

  std::vector<int> eigenvalues;
  eigenvalues.push_back(0.);

  element_type Hilbert_space(names, eigenvalues);

  for (int_type s = 0; s < num_fock_states; ++s) {
    psi_state<parameter_type, ed_options> psi(s);
    Hilbert_space.insert(psi);
  }

  Hilbert_spaces.push_back(Hilbert_space);
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::initialize_with_occupation_symmetry() {
  std::vector<element_type>& Hilbert_spaces = get_elements();

  int_type N_occ_max = b_dmn::dmn_size() * s_dmn::dmn_size() * RClusterDmn::dmn_size();

  for (int_type N_occ = 0; N_occ <= N_occ_max; ++N_occ) {
    std::vector<std::string> names;
    names.push_back("occupation");

    std::vector<int> eigenvalues;
    eigenvalues.push_back(N_occ);

    Hilbert_spaces.push_back(element_type(names, eigenvalues));
  }

  int_type number_of_fock_states = 1ul << N_occ_max;

  for (int_type s = 0; s < number_of_fock_states; ++s) {
    psi_state<parameter_type, ed_options> psi(s);
    int_type N_occ = psi.occupation_number();
    Hilbert_spaces[N_occ].insert(psi);
  }
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::apply_magnetization_symmetry() {
  std::vector<element_type>& Hilbert_spaces = get_elements();

  int_type N_occ_max = b_dmn::dmn_size() * s_dmn::dmn_size() * RClusterDmn::dmn_size();

  for (int_type N_occ = 0; N_occ <= N_occ_max; ++N_occ) {
    element_type HS(*Hilbert_spaces.begin());
    Hilbert_spaces.erase(Hilbert_spaces.begin());

    int_type number_of_subspaces;

    if (N_occ <= N_occ_max / 2)
      number_of_subspaces = N_occ + 1;
    else
      number_of_subspaces = N_occ_max - N_occ + 1;

    int magnetization_of_subspace = -(number_of_subspaces - 1);

    for (int_type j = 0; j < number_of_subspaces; ++j) {
      std::vector<std::string> names;
      names.push_back("occupation");
      names.push_back("magnetization");

      std::vector<int> eigenvalues;
      eigenvalues.push_back(N_occ);
      eigenvalues.push_back(magnetization_of_subspace);

      element_type new_subspace(names, eigenvalues);

      for (typename std::vector<psi_state<parameter_type, ed_options>>::iterator it =
               HS.psi_states_begin();
           it != HS.psi_states_end(); ++it) {
        if (it->magnetization() == magnetization_of_subspace)
          new_subspace.insert(*it);
      }

      Hilbert_spaces.push_back(new_subspace);
      magnetization_of_subspace += 2;
    }
  }
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::apply_translation_symmetry(std::string ED_method) {
  bool create_subspaces = true;
  if (ED_method == "block-diagonal")
    create_subspaces = false;

  std::vector<element_type>& Hilbert_spaces = get_elements();

  int num_states = b_dmn::dmn_size() * s_dmn::dmn_size() * RClusterDmn::dmn_size();

  int DIMENSION = r_cluster_type::DIMENSION;

  std::vector<std::vector<double>>& basis = r_cluster_type::get_basis_vectors();
  std::vector<std::vector<double>>& super_basis = r_cluster_type::get_super_basis_vectors();
  std::vector<std::vector<double>>& elements = r_cluster_type::get_elements();

  std::vector<int> index(DIMENSION);

  std::vector<std::vector<int>> applied_symmetries;

  std::vector<int> identity(num_states);
  for (int i = 0; i < identity.size(); ++i) {
    identity[i] = i;
  }
  applied_symmetries.push_back(identity);

  for (int k = 0; k < DIMENSION; ++k) {
    std::vector<double> basis_vec =
        domains::cluster_operations::translate_inside_cluster(basis[k], super_basis);

    index[k] = domains::cluster_operations::index(basis_vec, elements, domains::BRILLOUIN_ZONE);

    int Nc = elements.size();

    std::vector<int> indices(elements.size());
    for (int l = 0; l < Nc; l++) {
      indices[l] = r_cluster_type::add(l, index[k]);
    }

    std::vector<int> permutation_vector(num_states);

    int idx = 0;
    for (int r = 0; r < RClusterDmn::dmn_size(); ++r) {
      for (int s = 0; s < s_dmn::dmn_size(); ++s) {
        for (int b = 0; b < b_dmn::dmn_size(); ++b) {
          permutation_vector[idx] =
              indices[r] * s_dmn::dmn_size() * b_dmn::dmn_size() + s * b_dmn::dmn_size() + b;
          ++idx;
        }
      }
    }

    if (find(applied_symmetries.begin(), applied_symmetries.end(), permutation_vector) ==
        applied_symmetries.end()) {
      applied_symmetries.push_back(permutation_vector);

      symmetry_operation<parameter_type, ed_options> Op;
      Op.initialize(permutation_vector);

      std::vector<element_type> new_Hilbert_spaces;

      for (int i = 0; i < get_size(); ++i) {
        element_type old_subspace(*(Hilbert_spaces.begin() + i));

        factorize(Op, old_subspace, new_Hilbert_spaces, create_subspaces, i);
      }

      Hilbert_spaces.swap(new_Hilbert_spaces);
    }
  }
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::print_subspaces(bool full) {
  std::vector<element_type>& Hilbert_spaces = get_elements();
  std::vector<std::string> names = Hilbert_spaces[0].get_name();

  std::cout << "Number of subspaces: " << Hilbert_spaces.size() << "\n\t#states\t#phis";
  for (int i = 0; i < names.size(); ++i) {
    std::cout << "\t" << names[i];
  }
  std::cout << "\n";
  for (int i = 0; i < Hilbert_spaces.size(); ++i) {
    std::cout << "#" << i << ":\t";
    Hilbert_spaces[i].print(full);
  }
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::sort_wrt_size() {
  std::vector<element_type>& Hilbert_spaces = get_elements();
  std::sort(Hilbert_spaces.begin(), Hilbert_spaces.end(), operator><parameter_type, ed_options>);
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::apply_rotation_symmetry(std::string /*symmetries*/,
                                                                     std::string ED_method) {
  bool create_subspaces = true;
  if (ED_method == "block-diagonal")
    create_subspaces = false;

  std::vector<element_type>& Hilbert_spaces = get_elements();

  func::function<std::pair<int, int>, r_symmetry_matrix_dmn_t>& r_symmetry_matrix =
      r_symmetry_type::get_symmetry_matrix();

  int num_states = b_dmn::dmn_size() * s_dmn::dmn_size() * RClusterDmn::dmn_size();

  std::vector<std::vector<int>> applied_symmetries;

  std::vector<int> identity(num_states);
  for (int i = 0; i < identity.size(); ++i) {
    identity[i] = i;
  }
  applied_symmetries.push_back(identity);

  for (int l = 1; l < 2; ++l) {
    std::vector<int> permutation_vector(num_states);

    int index = 0;
    for (int r = 0; r < RClusterDmn::dmn_size(); ++r) {
      for (int s = 0; s < s_dmn::dmn_size(); ++s) {
        for (int b = 0; b < b_dmn::dmn_size(); ++b) {
          permutation_vector[index] =
              r_symmetry_matrix(r, b, l).first * s_dmn::dmn_size() * b_dmn::dmn_size() +
              s * b_dmn::dmn_size() + r_symmetry_matrix(r, b, l).second;
          ++index;
        }
      }
    }

    if (find(applied_symmetries.begin(), applied_symmetries.end(), permutation_vector) ==
        applied_symmetries.end()) {
      std::cout << "apply rotation symmetry #" << l << std::endl;

      applied_symmetries.push_back(permutation_vector);

      symmetry_operation<parameter_type, ed_options> Op;
      Op.initialize(permutation_vector);

      std::vector<element_type> new_Hilbert_spaces;

      for (int i = 0; i < get_size(); ++i) {
        element_type old_subspace(*(Hilbert_spaces.begin() + i));
        factorize(Op, old_subspace, new_Hilbert_spaces, create_subspaces, i);
      }

      Hilbert_spaces.swap(new_Hilbert_spaces);
    }
  }
}

template <typename parameter_type, typename ed_options>
template <class symmetry_operation>
void Fock_space<parameter_type, ed_options>::factorize(symmetry_operation& Op, element_type& subspace,
                                                       std::vector<element_type>& new_Hilbert_spaces,
                                                       bool create_subspaces, int k) {
  std::vector<element_type> factorized_subspace;

  std::vector<std::string> names = subspace.get_name();
  names.push_back(Op.get_name());

  std::vector<int> eigenvalues = subspace.get_eigenvalues();

  for (int i = 0; i < Op.get_order(); ++i) {
    eigenvalues.push_back(i);

    element_type new_subspace(names, eigenvalues);
    factorized_subspace.push_back(new_subspace);

    eigenvalues.pop_back();
  }

  subspace.sort();

  for (int l = 0; l < subspace.size(); ++l) {
    psi_state<parameter_type, ed_options>& psi0 = subspace.get_element(l);

    if (psi0.is_marked() == false) {
      psi0.mark();

      psi_state<parameter_type, ed_options> psi_tmp = psi0;

      int order = 1;

      Op.execute(psi_tmp, true);

      while (!(psi_tmp == psi0)) {
        bool state_found = subspace.mark_state(psi_tmp);

        if (not state_found) {
          std::cout << "\n\n\t psi0 \n\n";
          psi0.print();

          std::cout << "\n\n\tstate not found !!!\n\n";

          for (int l = 0; l < subspace.size(); l++)
            subspace.get_element(l).print();

          throw std::logic_error(__FUNCTION__);
        }

        Op.execute(psi_tmp, true);

        ++order;
      }

      complex_type eval_factor = exp(complex_type(0, 2. * M_PI / order));
      complex_type eval = 1.;

      for (int i = 0; i < Op.get_order(); i += Op.get_order() / order) {
        psi_state<parameter_type, ed_options> psi_prime(eval, order, Op, psi0, k);

        if (psi_prime.size() != 0) {
          factorized_subspace[i].insert(psi_prime);
        }

        eval *= eval_factor;

        if (std::abs(real(eval)) < ed_options::get_epsilon())
          eval.real(0.);

        if (std::abs(imag(eval)) < ed_options::get_epsilon())
          eval.imag(0.);
      }
    }
  }

  if (create_subspaces)  // Create subspaces
  {
    for (int i = 0; i < Op.get_order(); ++i) {
      if (!factorized_subspace[i].empty())
        new_Hilbert_spaces.push_back(factorized_subspace[i]);
    }
  }
  else  // Create NO subspaces
  {
    int index = 0;
    while (factorized_subspace[index].empty())
      ++index;

    for (int i = index + 1; i < Op.get_order(); ++i) {
      for (int j = 0; j < factorized_subspace[i].size(); ++j) {
        factorized_subspace[index].insert(factorized_subspace[i].get_element(j));
      }
    }
    new_Hilbert_spaces.push_back(factorized_subspace[index]);
  }
}

template <typename parameter_type, typename ed_options>
bool Fock_space<parameter_type, ed_options>::check_orthogonality() {
  complex_type res = 0.;

  std::vector<element_type>& Hilbert_spaces = get_elements();

  for (int i = 0; i < get_size(); ++i) {
    element_type& subspace_1 = Hilbert_spaces[i];
    for (int j = i + 1; j < get_size(); ++j) {
      element_type& subspace_2 = Hilbert_spaces[j];

      for (int k = 0; k < subspace_1.size(); ++k) {
        psi_state<parameter_type, ed_options>& psi_1 = subspace_1.get_element(k);
        for (int l = 0; l < subspace_2.size(); ++l) {
          psi_state<parameter_type, ed_options>& psi_2 = subspace_2.get_element(l);

          res = scalar_product(psi_1, psi_2);
          if (res != complex_type(0.))
            return false;
        }
      }
    }
  }
  return true;
}

template <typename parameter_type, typename ed_options>
void Fock_space<parameter_type, ed_options>::initialize_rep() {
  std::vector<element_type>& Hilbert_spaces = get_elements();

  for (int i = 0; i < get_size(); ++i)
    Hilbert_spaces[i].initialize_rep();
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_FOCK_SPACE_HPP
