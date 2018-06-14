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
// Material lattice specialization for CuO2.
//
// TODO: - Replace get_LDA_Hamiltonians with intialize_H_0 function (see material_lattice_NiO.inc).
//       - Use correct index of origin in initialize_H_interaction (see e.g.
//         material_lattice_NiO.inc).

template <typename point_group_type>
class material_lattice<CuO2, point_group_type> {
public:
  const static int DIMENSION = 2;
  const static int BANDS = 3;

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
};

template <typename point_group_type>
double* material_lattice<CuO2, point_group_type>::initialize_r_DCA_basis() {
  const static double a = 1.0;

  static double* r_DCA = new double[DIMENSION * DIMENSION];

  r_DCA[0] = a;
  r_DCA[2] = 0;
  r_DCA[1] = 0;
  r_DCA[3] = a;

  return r_DCA;
}

template <typename point_group_type>
double* material_lattice<CuO2, point_group_type>::initialize_r_LDA_basis() {
  const static double a = 1.0;

  static double* r_LDA = new double[DIMENSION * DIMENSION];

  r_LDA[0] = a;
  r_LDA[2] = 0;
  r_LDA[1] = 0;
  r_LDA[3] = a;

  return r_LDA;
}

template <typename point_group_type>
std::vector<int> material_lattice<CuO2, point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;
  flavors[2] = 2;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> material_lattice<CuO2, point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> material_lattice<
    CuO2, point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <class domain, class parameters_type>
void material_lattice<CuO2, point_group_type>::initialize_H_interaction(
    func::function<double, domain>& H_i, parameters_type& parameters) {
  std::vector<std::vector<double>> U;

  io::CSVReader::execute(parameters.get_U_ij_file_name(), U);

  for (int i = 0; i < 2 * BANDS; i++)
    for (int j = 0; j < 2 * BANDS; j++)
      H_i(i, j, 0) = U[i][j];

  /*
  {
    std::stringstream ss;
    ss.precision(6);
    ss<<scientific;

    ss << "\n\n U_ij \n\n";
    for(int i=0; i<2*BANDS; i++){
      for(int j=0; j<2*BANDS; j++){
        ss << "\t" << H_i(i,j,0);
      }
      ss << "\n";
    }
    ss << "\n\n";

    cout << ss.str();
  }
  */
}

template <typename point_group_type>
template <class domain>
void material_lattice<CuO2, point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  for (int s = 0; s < 2; s++)
    for (int i = 0; i < BANDS; i++)
      H_symmetries(i, s, i, s) = get_flavors()[i];
}

template <typename point_group_type>
template <class parameters_type>
std::vector<std::vector<double>>& material_lattice<CuO2, point_group_type>::get_t_ij(
    parameters_type& parameters) {
  static std::vector<std::vector<double>> t_ij(0, std::vector<double>(0));

  static bool is_initialized = false;

  if (not is_initialized) {
    io::CSVReader::execute(parameters.get_t_ij_file_name(), t_ij);

    is_initialized = true;

    /*
    cout << "\n\n";
    for(size_t i=0; i<128; i++){
      cout << "\t";
      math::util::print(t_ij[i]);
      cout << "\n";
    }
    */
  }

  return t_ij;
}

template <typename point_group_type>
template <class parameters_type>
std::complex<double> material_lattice<CuO2, point_group_type>::get_LDA_Hamiltonians(
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

      if (b_i == b1 and b_j == b2 and t[4] == 0)  // keep only the in-plane hopping
      {
        double r_x = t[2] * r_LDA[0] + t[3] * r_LDA[2];
        double r_y = t[2] * r_LDA[1] + t[3] * r_LDA[3];

        val.real(t[5]);
        val.imag(t[6]);

        H_LDA += val * std::exp(I * (k[0] * r_x + k[1] * r_y));
      }
    }
  }

  return H_LDA;
}
