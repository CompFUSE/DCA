// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// Material lattice specialization for NiO.
//
// TODO: Remove code duplication in the two NiO specializations (symmetric and unsymmetric).

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

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandSpinDmn, typename RDmn, typename parameters_type>
  static void initialize_H_interaction(
      func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are three-dimensional (access through index 0, 1, and 2).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initialize_H_0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);

  static double latticeConstant() {
    return a;
  }

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
template <typename BandSpinDmn, typename RDmn, typename parameters_type>
void material_lattice<NiO_symmetric, point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandSpinDmn::dmn_size() != 2 * BANDS)
    throw std::logic_error("Band-spin domain size must be 8*2=16.");

  const int origin = RDmn::parameter_type::origin_index();

  std::vector<std::vector<double>> U;
  io::CSVReader::execute(parameters.get_U_ij_file_name(), U);

  H_interaction = 0.;

  for (int i = 0; i < 2 * BANDS; i++)
    for (int j = 0; j < 2 * BANDS; j++)
      H_interaction(i, j, origin) = U[i][j];
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
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void material_lattice<NiO_symmetric, point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("NiO has 8 bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const static std::complex<double> I(0, 1);

  const auto& k_vecs = KDmn::get_elements();
  double* r_LDA = initialize_r_LDA_basis();

  std::vector<std::vector<double>>& t_ij = get_t_ij(parameters);

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];

    for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1) {
      for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2) {
        if (s1 == s2) {
          for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1) {
            for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2) {
              std::complex<double> val;

              for (std::size_t i = 0; i < t_ij.size(); i++) {
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

                  H_0(b1, s1, b2, s2, k_ind) +=
                      val * std::exp(I * (k[0] * r_x + k[1] * r_y + k[2] * r_z)) / multiplicity;
                }
              }
            }
          }
        }
      }
    }
  }
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

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandSpinDmn, typename RDmn, typename parameters_type>
  static void initialize_H_interaction(
      func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are three-dimensional (access through index 0, 1, and 2).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initialize_H_0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);

  static double latticeConstant() {
    return a;
  }

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
template <typename BandSpinDmn, typename RDmn, typename parameters_type>
void material_lattice<NiO_unsymmetric, point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandSpinDmn::dmn_size() != 2 * BANDS)
    throw std::logic_error("Band-spin domain size must be 8*2=16.");

  const int origin = RDmn::parameter_type::origin_index();

  std::vector<std::vector<double>> U;
  io::CSVReader::execute(parameters.get_U_ij_file_name(), U);

  H_interaction = 0.;

  for (int i = 0; i < 2 * BANDS; i++)
    for (int j = 0; j < 2 * BANDS; j++)
      H_interaction(i, j, origin) = U[i][j];
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
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void material_lattice<NiO_unsymmetric, point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("NiO has 8 bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const static std::complex<double> I(0, 1);

  const auto& k_vecs = KDmn::get_elements();
  double* r_LDA = initialize_r_LDA_basis();

  std::vector<std::vector<double>>& t_ij = get_t_ij(parameters);

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];

    for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1) {
      for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2) {
        if (s1 == s2) {
          for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1) {
            for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2) {
              std::complex<double> val;

              for (std::size_t i = 0; i < t_ij.size(); i++) {
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

                  H_0(b1, s1, b2, s2, k_ind) +=
                      val * std::exp(I * (k[0] * r_x + k[1] * r_y + k[2] * r_z)) / multiplicity;
                }
              }
            }
          }
        }
      }
    }
  }
}
