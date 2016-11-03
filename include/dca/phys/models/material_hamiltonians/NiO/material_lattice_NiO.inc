// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// Material lattice specialization for NiO.

// NiO symmetric
template <typename point_group_type>
class material_lattice<NiO_symmetric, point_group_type> {
public:
  static constexpr int DIMENSION = 3;
  static constexpr int BANDS = 8;

  typedef domains::no_symmetry<DIMENSION> LDA_point_group;
  typedef point_group_type DCA_point_group;

  static double* initialize_r_DCA_basis();
  static double* initialize_r_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(func::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  template <class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, std::vector<double> k,
                                                   int b1, int s1, int b2, int s2);

private:
  template <class parameters_type>
  static std::vector<std::vector<double>>& get_t_ij(parameters_type& parameters);

  static constexpr double a = 7.89358;  // 3.946789070;
};

template <typename point_group_type>
double* material_lattice<NiO_symmetric, point_group_type>::initialize_r_DCA_basis() {
  static std::array<double, 9> r_DCA;

  r_DCA[0] = 0;
  r_DCA[3] = a / 2.;
  r_DCA[6] = a / 2.;
  r_DCA[1] = a / 2.;
  r_DCA[4] = 0;
  r_DCA[7] = a / 2.;
  r_DCA[2] = a / 2.;
  r_DCA[5] = a / 2.;
  r_DCA[8] = 0;

  return r_DCA.data();
}

template <typename point_group_type>
double* material_lattice<NiO_symmetric, point_group_type>::initialize_r_LDA_basis() {
  static std::array<double, 9> r_LDA;

  r_LDA[0] = 0;
  r_LDA[3] = a / 2.;
  r_LDA[6] = a / 2.;
  r_LDA[1] = a / 2.;
  r_LDA[4] = 0;
  r_LDA[7] = a / 2.;
  r_LDA[2] = a / 2.;
  r_LDA[5] = a / 2.;
  r_LDA[8] = 0;

  return r_LDA.data();
}

template <typename point_group_type>
std::vector<int> material_lattice<NiO_symmetric, point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);
  flavors[0] = 0;
  flavors[1] = 0;
  flavors[2] = 1;
  flavors[3] = 0;
  flavors[4] = 1;
  flavors[5] = 2;
  flavors[6] = 2;
  flavors[7] = 2;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> material_lattice<NiO_symmetric, point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> material_lattice<
    NiO_symmetric, point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <class domain, class parameters_type>
void material_lattice<NiO_symmetric, point_group_type>::initialize_H_interaction(
    func::function<double, domain>& H_i, parameters_type& parameters) {
  std::vector<std::vector<double>> U;

  io::CSVReader::execute(parameters.get_U_ij_file_name(), U);

  for (int i = 0; i < 2 * BANDS; i++)
    for (int j = 0; j < 2 * BANDS; j++)
      H_i(i, j, 0) = U[i][j];
}

template <typename point_group_type>
template <class domain>
void material_lattice<NiO_symmetric, point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  for (int s = 0; s < 2; s++)
    for (int i = 0; i < BANDS; i++)
      H_symmetries(i, s, i, s) = get_flavors()[i];
}

template <typename point_group_type>
template <class parameters_type>
std::vector<std::vector<double>>& material_lattice<NiO_symmetric, point_group_type>::get_t_ij(
    parameters_type& parameters) {
  static std::vector<std::vector<double>> t_ij(0, std::vector<double>(0));

  static bool is_initialized = false;

  if (not is_initialized) {
    io::CSVReader::execute(parameters.get_t_ij_file_name(), t_ij);

    is_initialized = true;
  }

  return t_ij;
}

template <typename point_group_type>
template <class parameters_type>
std::complex<double> material_lattice<NiO_symmetric, point_group_type>::get_LDA_Hamiltonians(
    parameters_type& parameters, std::vector<double> k, int b1, int s1, int b2, int s2) {
  const static std::complex<double> I(0, 1);

  double* r_LDA = initialize_r_LDA_basis();

  std::vector<std::vector<double>>& t_ij = get_t_ij(parameters);

  std::complex<double> H_LDA = 0.;

  if (s1 == s2) {
    std::complex<double> val;

    for (size_t i = 0; i < t_ij.size(); i++) {
      std::vector<double>& t = t_ij[i];

      assert(t.size() == 8);

      int b_i = t[0];
      int b_j = t[1];

      if (b_i == b1 and b_j == b2) {
        double r_x = t[2] * r_LDA[0] + t[3] * r_LDA[3] + t[4] * r_LDA[6];
        double r_y = t[2] * r_LDA[1] + t[3] * r_LDA[4] + t[4] * r_LDA[7];
        double r_z = t[2] * r_LDA[2] + t[3] * r_LDA[5] + t[4] * r_LDA[8];

        val.real(t[5]);
        val.imag(t[6]);

        double multiplicity = t[7];

        H_LDA += val * std::exp(I * (k[0] * r_x + k[1] * r_y + k[2] * r_z)) / multiplicity;
      }
    }
  }

  return H_LDA;
}

// NiO unsymmetric
template <typename point_group_type>
class material_lattice<NiO_unsymmetric, point_group_type> {
public:
  static constexpr int DIMENSION = 3;
  static constexpr int BANDS = 8;

  typedef domains::no_symmetry<DIMENSION> LDA_point_group;
  typedef point_group_type DCA_point_group;

  static double* initialize_r_DCA_basis();
  static double* initialize_r_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(func::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  template <class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, std::vector<double> k,
                                                   int b1, int s1, int b2, int s2);

private:
  template <class parameters_type>
  static std::vector<std::vector<double>>& get_t_ij(parameters_type& parameters);

  static constexpr double a = 7.89358;  // 3.946789070;
};

template <typename point_group_type>
double* material_lattice<NiO_unsymmetric, point_group_type>::initialize_r_DCA_basis() {
  static std::array<double, 9> r_DCA;

  r_DCA[0] = 0;
  r_DCA[3] = a / 2.;
  r_DCA[6] = a / 2.;
  r_DCA[1] = a / 2.;
  r_DCA[4] = 0;
  r_DCA[7] = a / 2.;
  r_DCA[2] = a / 2.;
  r_DCA[5] = a / 2.;
  r_DCA[8] = 0;

  return r_DCA.data();
}

template <typename point_group_type>
double* material_lattice<NiO_unsymmetric, point_group_type>::initialize_r_LDA_basis() {
  static std::array<double, 9> r_LDA;

  r_LDA[0] = 0;
  r_LDA[3] = a / 2.;
  r_LDA[6] = a / 2.;
  r_LDA[1] = a / 2.;
  r_LDA[4] = 0;
  r_LDA[7] = a / 2.;
  r_LDA[2] = a / 2.;
  r_LDA[5] = a / 2.;
  r_LDA[8] = 0;

  return r_LDA.data();
}

template <typename point_group_type>
std::vector<int> material_lattice<NiO_unsymmetric, point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  for (int i = 0; i < BANDS; i++)
    flavors[i] = i;
  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> material_lattice<NiO_unsymmetric, point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> material_lattice<
    NiO_unsymmetric, point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <class domain, class parameters_type>
void material_lattice<NiO_unsymmetric, point_group_type>::initialize_H_interaction(
    func::function<double, domain>& H_i, parameters_type& parameters) {
  std::vector<std::vector<double>> U;

  io::CSVReader::execute(parameters.get_U_ij_file_name(), U);

  for (int i = 0; i < 2 * BANDS; i++)
    for (int j = 0; j < 2 * BANDS; j++)
      H_i(i, j, 0) = U[i][j];
}

template <typename point_group_type>
template <class domain>
void material_lattice<NiO_unsymmetric, point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  for (int s = 0; s < 2; s++)
    for (int i = 0; i < BANDS; i++)
      H_symmetries(i, s, i, s) = get_flavors()[i];
}

template <typename point_group_type>
template <class parameters_type>
std::vector<std::vector<double>>& material_lattice<NiO_unsymmetric, point_group_type>::get_t_ij(
    parameters_type& parameters) {
  static std::vector<std::vector<double>> t_ij(0, std::vector<double>(0));

  static bool is_initialized = false;

  if (not is_initialized) {
    io::CSVReader::execute(parameters.get_t_ij_file_name(), t_ij);

    is_initialized = true;
  }

  return t_ij;
}

template <typename point_group_type>
template <class parameters_type>
std::complex<double> material_lattice<NiO_unsymmetric, point_group_type>::get_LDA_Hamiltonians(
    parameters_type& parameters, std::vector<double> k, int b1, int s1, int b2, int s2) {
  const static std::complex<double> I(0, 1);

  double* r_LDA = initialize_r_LDA_basis();

  std::vector<std::vector<double>>& t_ij = get_t_ij(parameters);

  std::complex<double> H_LDA = 0.;

  if (s1 == s2) {
    std::complex<double> val;

    for (size_t i = 0; i < t_ij.size(); i++) {
      std::vector<double>& t = t_ij[i];

      assert(t.size() == 8);

      int b_i = t[0];
      int b_j = t[1];

      if (b_i == b1 and b_j == b2) {
        double r_x = t[2] * r_LDA[0] + t[3] * r_LDA[3] + t[4] * r_LDA[6];
        double r_y = t[2] * r_LDA[1] + t[3] * r_LDA[4] + t[4] * r_LDA[7];
        double r_z = t[2] * r_LDA[2] + t[3] * r_LDA[5] + t[4] * r_LDA[8];

        val.real(t[5]);
        val.imag(t[6]);

        double multiplicity = t[7];

        H_LDA += val * std::exp(I * (k[0] * r_x + k[1] * r_y + k[2] * r_z)) / multiplicity;
      }
    }
  }

  return H_LDA;
}
