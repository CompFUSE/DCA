// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// Material lattice specialization for FeSn.
//

#include <type_traits>

template <typename point_group_type>
class material_lattice<Material::FeSn, point_group_type, Is3D<point_group_type>> {
public:
  static constexpr bool complex_g0 = true;
  static constexpr bool spin_symmetric = false;

  static constexpr int DIMENSION = 3;
  static constexpr int BANDS = 15;

  static constexpr std::string_view material_name{"FeSn"};

  static constexpr Material my_Material = Material::FeSn;
  using my_PGT = point_group_type;

  typedef domains::no_symmetry<DIMENSION> LDA_point_group;
  typedef point_group_type DCA_point_group;

  constexpr static int transformationSignOfR(int, int, int) {
    return 1;
  }
  constexpr static int transformationSignOfK(int, int, int) {
    return 1;
  }

  static std::array<double, 9> initializeRDCABasis();
  static std::array<double, 9> initializeRLDABasis();

  static std::vector<int> flavors();
  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandSpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initializeHSymmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are three-dimensional (access through index 0, 1, and 2).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);

  static double latticeConstant() {
    return a;
  }

  static constexpr double lattice_constant = 5.2969999313;

private:
  template <class parameters_type>
  static std::vector<std::vector<double>>& get_t_ij(parameters_type& parameters);

  static constexpr double a = 5.2969999313;  // 3.946789070;
};

template <typename point_group_type>
std::array<double, 9> material_lattice<Material::FeSn, point_group_type,
                                       Is3D<point_group_type>>::initializeRDCABasis() {
  std::array<double, 9> r_DCA;

  r_DCA[0] = a;
  r_DCA[1] = 0;
  r_DCA[2] = 0;
  r_DCA[3] = -a * std::cos(M_PI / 3.);
  r_DCA[4] = a * std::sin(M_PI / 3.);
  r_DCA[5] = 0;
  r_DCA[6] = 0;
  r_DCA[7] = 0;
  r_DCA[8] = 4.481;

  return r_DCA;
}

template <typename point_group_type>
std::array<double, 9> material_lattice<Material::FeSn, point_group_type,
                                       Is3D<point_group_type>>::initializeRLDABasis() {
  std::array<double, 9> r_LDA;

  r_LDA[0] = a;
  r_LDA[1] = 0;
  r_LDA[2] = 0;
  r_LDA[3] = -a * std::cos(M_PI / 3.);
  r_LDA[4] = a * std::sin(M_PI / 3.);
  r_LDA[5] = 0;
  r_LDA[6] = 0;
  r_LDA[7] = 0;
  r_LDA[8] = 4.481;

  return r_LDA;
}

template <typename point_group_type>
std::vector<int> material_lattice<Material::FeSn, point_group_type, Is3D<point_group_type>>::flavors() {
  std::vector<int> flavors(BANDS);
  for (int i = 0; i < BANDS; i++)
    flavors[i] = i;
  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> material_lattice<Material::FeSn, point_group_type,
                                                  Is3D<point_group_type>>::aVectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  a_vecs[0] = {0., 0., 0.};
  a_vecs[1] = {0., 0., 0.};
  a_vecs[2] = {0., 0., 0.};
  a_vecs[3] = {0., 0., 0.};
  a_vecs[4] = {0., 0., 0.};

  a_vecs[5] = {-0.75 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[6] = {-0.75 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[7] = {-0.75 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[8] = {-0.75 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[9] = {-0.75 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};

  a_vecs[10] = {-0.25 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[11] = {-0.25 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[12] = {-0.25 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[13] = {-0.25 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};
  a_vecs[14] = {-0.25 * a, 0.5 * std::sin(M_PI / 3.) * a, 0.0};

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> material_lattice<
    Material::FeSn, point_group_type, Is3D<point_group_type>>::orbitalPermutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandSpinDmn, typename RDmn, typename parameters_type>
void material_lattice<Material::FeSn, point_group_type, Is3D<point_group_type>>::initializeHInteraction(
    func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandSpinDmn::dmn_size() != 2 * BANDS) {
    std::ostringstream msg;
    msg << material_name << " has BandSpin domin of " << BANDS * 2 << " but BandDmn size is " << BandSpinDmn::dmn_size()
        << "!\n";
    throw std::logic_error(msg.str());
  }
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
void material_lattice<Material::FeSn, point_group_type, Is3D<point_group_type>>::initializeHSymmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  for (int s = 0; s < 2; s++)
    for (int i = 0; i < BANDS; i++)
      H_symmetries(i, s, i, s) = flavors()[i];
}

template <typename point_group_type>
template <class parameters_type>
std::vector<std::vector<double>>& material_lattice<
    Material::FeSn, point_group_type, Is3D<point_group_type>>::get_t_ij(parameters_type& parameters) {
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
void material_lattice<Material::FeSn, point_group_type, Is3D<point_group_type>>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS) {
    std::ostringstream msg;
    msg << material_name << " has " << BANDS << " but BandDmn size is " << BandDmn::dmn_size()
        << "!\n";
    throw std::logic_error(msg.str());
  }
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const std::complex<double> I(0, 1);

  const auto& k_vecs = KDmn::get_elements();
  std::array<double, 9> r_LDA = initializeRLDABasis();

  std::vector<std::vector<double>>& t_ij = get_t_ij(parameters);

  H_0 = ScalarType(0);

  if (t_ij[0].size() != 8)
    throw std::runtime_error("each row in the t_ij file must have 8 elements!");

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    std::complex<double> val;
    for (std::size_t i = 0; i < t_ij.size(); i++) {
      std::vector<double>& t = t_ij[i];

      assert(t.size() == 8);

      int b_i = t[0];
      int b_j = t[1];

      auto notInBndDmn = [](int b) -> bool { return (b >= BandDmn::dmn_size() || b < 0); };

      if (notInBndDmn(b_i) || notInBndDmn(b_j))
        throw std::runtime_error(
            "the t_ij input cannot contain bands not in the models band domain!");
      double r_x = t[2] * r_LDA[0] + t[3] * r_LDA[3] + t[4] * r_LDA[6];
      double r_y = t[2] * r_LDA[1] + t[3] * r_LDA[4] + t[4] * r_LDA[7];
      double r_z = t[2] * r_LDA[2] + t[3] * r_LDA[5] + t[4] * r_LDA[8];

      val.real(t[5]);
      val.imag(t[6]);

      double multiplicity = t[7];
      ScalarType h0_elem = val * std::exp(I * (k[0] * r_x + k[1] * r_y + k[2] * r_z)) / multiplicity;
      H_0(b_i, 0, b_j, 0, k_ind) += h0_elem;
      H_0(b_i, 1, b_j, 1, k_ind) += h0_elem;
    }
  }
}

extern template class material_lattice<Material::FeSn, dca::phys::domains::no_symmetry<3>>;
